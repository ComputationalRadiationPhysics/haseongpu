#include <openPMD/openPMD.hpp>
#include <openpmd/OpenPmdParser.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace io = openPMD;

namespace
{
    namespace field
    {
        constexpr char const* numberOfPoints = "number_of_points";
        constexpr char const* numberOfCells = "number_of_cells";
        constexpr char const* numberOfLevels = "number_of_levels";
        constexpr char const* thickness = "thickness";
        constexpr char const* nTot = "n_tot";
        constexpr char const* crystalTFluo = "crystal_t_fluo";
        constexpr char const* claddingNumber = "cladding_number";
        constexpr char const* claddingAbsorption = "cladding_absorption";
        constexpr char const* minRaysPerSample = "min_rays_per_sample";
        constexpr char const* maxRaysPerSample = "max_rays_per_sample";
        constexpr char const* mseThreshold = "mse_threshold";
        constexpr char const* useReflections = "use_reflections";
        constexpr char const* spectralResolution = "spectral_resolution";
        constexpr char const* monochromatic = "monochromatic";
        constexpr char const* maxSigmaAbsorption = "max_sigma_absorption";
        constexpr char const* maxSigmaEmission = "max_sigma_emission";
        constexpr char const* repetitions = "repetitions";
        constexpr char const* adaptiveSteps = "adaptive_steps";
        constexpr char const* maxGpus = "max_gpus";
        constexpr char const* backend = "backend";
        constexpr char const* parallelMode = "parallel_mode";
        constexpr char const* minSampleRange = "min_sample_range";
        constexpr char const* maxSampleRange = "max_sample_range";
        constexpr char const* rngSeed = "rng_seed";
        constexpr char const* writeVtk = "write_vtk";
        constexpr char const* devices = "devices";
    } // namespace field

    constexpr char const* OPENPMD_SST_CONFIG = R"(
{
  "backend": "adios2",
  "adios2": {
    "engine": {
      "type": "sst",
      "parameters": {
        "DataTransport": "WAN",
        "OpenTimeoutSecs": "600"
      }
    }
  }
})";

    constexpr char const* OPENPMD_DEFAULT_CONFIG = "{}";
    constexpr char const* HASE_TRANSPORT_VERSION = "0.1";

    bool hasSuffix(std::string_view value, std::string_view suffix)
    {
        return value.size() >= suffix.size() && value.substr(value.size() - suffix.size()) == suffix;
    }

    char const* seriesConfig(std::string const& stream)
    {
        return hasSuffix(stream, ".sst") ? OPENPMD_SST_CONFIG : OPENPMD_DEFAULT_CONFIG;
    }

    template<typename T>
    T attribute(io::Attributable const& obj, std::string const& name)
    {
        return obj.getAttribute(name).get<T>();
    }

    template<typename T>
    T attributeOr(io::Attributable const& obj, std::string const& name, T fallback)
    {
        if(obj.containsAttribute(name))
        {
            return attribute<T>(obj, name);
        }
        return fallback;
    }

    std::size_t elementCount(io::Extent const& extent)
    {
        return std::accumulate(extent.begin(), extent.end(), std::size_t{1}, std::multiplies<std::size_t>{});
    }

    [[noreturn]] void validationError(std::string const& field, std::string const& message)
    {
        throw std::runtime_error("openPMD validation error for '" + field + "': " + message);
    }

    std::string entityFromAxes(std::vector<std::string> const& axes);

    template<typename T>
    std::string vectorString(std::vector<T> const& values)
    {
        std::ostringstream out;
        out << "[";
        for(std::size_t i = 0; i < values.size(); ++i)
        {
            if(i != 0)
            {
                out << ", ";
            }
            out << values[i];
        }
        out << "]";
        return out.str();
    }

    std::vector<unsigned long long> extentVector(io::Extent const& extent)
    {
        return std::vector<unsigned long long>(extent.begin(), extent.end());
    }

    void validateExtent(std::string const& name, io::Extent const& extent, io::Extent const& expected)
    {
        if(extent != expected)
        {
            auto const actualVector = extentVector(extent);
            auto const expectedVector = extentVector(expected);
            validationError(
                name,
                "extent mismatch (expected " + vectorString(expectedVector) + ", got " + vectorString(actualVector)
                    + ")");
        }
    }

    template<typename T>
    void validateAttribute(
        std::string const& name,
        io::Attributable const& obj,
        std::string const& attributeName,
        T const& expected)
    {
        if(!obj.containsAttribute(attributeName))
        {
            validationError(name, "missing required attribute '" + attributeName + "'");
        }
        auto const actual = attribute<T>(obj, attributeName);
        if(actual != expected)
        {
            validationError(name, "attribute '" + attributeName + "' mismatch");
        }
    }

    std::vector<std::string> splitAxesString(std::string const& value)
    {
        std::vector<std::string> axes;
        std::stringstream stream(value);
        std::string axis;
        while(std::getline(stream, axis, ','))
        {
            if(!axis.empty())
            {
                axes.push_back(axis);
            }
        }
        return axes;
    }

    void validateAxesAttribute(
        std::string const& name,
        io::Attributable const& obj,
        std::vector<std::string> const& expected)
    {
        if(obj.containsAttribute("haseAxes"))
        {
            try
            {
                if(attribute<std::vector<std::string>>(obj, "haseAxes") == expected)
                {
                    return;
                }
            }
            catch(std::exception const&)
            {
            }
        }

        if(obj.containsAttribute("haseAxesString"))
        {
            // SST can carry Python string-list attributes as an empty scalar
            // string. Accept the scalar fallback while keeping haseAxes as the
            // canonical metadata for file-backed backends.
            if(splitAxesString(attribute<std::string>(obj, "haseAxesString")) == expected)
            {
                return;
            }
        }

        validationError(name, "attribute 'haseAxes' mismatch");
    }

    void validateTransportVersion(std::string const& name, io::Attributable const& obj)
    {
        if(obj.containsAttribute("haseTransportVersion"))
        {
            validateAttribute(name, obj, "haseTransportVersion", std::string{HASE_TRANSPORT_VERSION});
            return;
        }
        if(obj.containsAttribute("haseSchemaVersion"))
        {
            // Legacy transport files used this name before the docs separated
            // Python schemas from the openPMD transport layout.
            validateAttribute(name, obj, "haseSchemaVersion", std::string{HASE_TRANSPORT_VERSION});
            return;
        }
        validationError(name, "missing required attribute 'haseTransportVersion'");
    }

    void validateHaseMetadata(
        std::string const& name,
        io::Mesh const& record,
        std::vector<std::string> const& axes,
        std::vector<unsigned long long> const& primitiveShape,
        bool dynamic,
        bool backendRequired,
        std::string const& unit)
    {
        validateTransportVersion(name, record);
        validateAttribute(name, record, "haseEntity", entityFromAxes(axes));
        validateAxesAttribute(name, record, axes);
        validateAttribute(name, record, "haseLayoutOrder", std::string{"backendFlat"});
        validateAttribute(name, record, "hasePrimitiveShape", primitiveShape);
        validateAttribute(name, record, "haseStatic", !dynamic);
        validateAttribute(name, record, "haseDynamic", dynamic);
        validateAttribute(name, record, "haseBackendRequired", backendRequired);
        validateAttribute(name, record, "haseUnit", unit);
    }

    void validateAxisLabels(std::string const& name, io::Mesh const& record, std::vector<std::string> const& expected)
    {
        auto const actual = record.axisLabels();
        if(actual == expected)
        {
            return;
        }
        if(record.containsAttribute("haseAxisLabelsString"))
        {
            // SST can drop openPMD axisLabels on streamed mesh records. Accept
            // the scalar fallback while keeping axisLabels canonical elsewhere.
            if(splitAxesString(attribute<std::string>(record, "haseAxisLabelsString")) == expected)
            {
                return;
            }
        }
        if(actual != expected)
        {
            validationError(
                name,
                "axis labels mismatch (expected " + vectorString(expected) + ", got " + vectorString(actual) + ")");
        }
    }

    template<typename T>
    std::vector<T> loadScalar(
        io::Series& series,
        io::Iteration& iteration,
        std::string const& name,
        io::Extent const& expectedExtent,
        std::vector<std::string> const& axes,
        std::vector<unsigned long long> const& primitiveShape,
        bool dynamic,
        bool backendRequired,
        std::string const& unit = "1")
    {
        if(!iteration.meshes.contains(name))
        {
            validationError(name, "missing required mesh record");
        }
        auto record = iteration.meshes[name];
        validateHaseMetadata(name, record, axes, primitiveShape, dynamic, backendRequired, unit);
        validateAxisLabels(name, record, {"flatIndex"});
        if(!record.contains(io::MeshRecordComponent::SCALAR))
        {
            validationError(name, "missing required scalar component");
        }
        auto component = record[io::MeshRecordComponent::SCALAR];
        auto const expectedDatatype = io::determineDatatype<T>();
        auto const actualDatatype = component.getDatatype();
        if(actualDatatype != expectedDatatype)
        {
            validationError(
                name,
                "dtype mismatch (expected " + io::datatypeToString(expectedDatatype) + ", got "
                    + io::datatypeToString(actualDatatype) + ")");
        }
        auto const extent = component.getExtent();
        validateExtent(name, extent, expectedExtent);
        auto chunk = component.loadChunk<T>();
        series.flush();
        return std::vector<T>(chunk.get(), chunk.get() + elementCount(extent));
    }

    template<typename T>
    std::vector<T> loadComponent(
        io::Series& series,
        io::Iteration& iteration,
        std::string const& name,
        std::string const& componentName,
        io::Extent const& expectedExtent)
    {
        if(!iteration.meshes.contains(name))
        {
            validationError(name, "missing required mesh record");
        }
        auto record = iteration.meshes[name];
        if(!record.contains(componentName))
        {
            validationError(name + "/" + componentName, "missing required record component");
        }
        auto component = record[componentName];
        auto const expectedDatatype = io::determineDatatype<T>();
        auto const actualDatatype = component.getDatatype();
        if(actualDatatype != expectedDatatype)
        {
            validationError(
                name + "/" + componentName,
                "dtype mismatch (expected " + io::datatypeToString(expectedDatatype) + ", got "
                    + io::datatypeToString(actualDatatype) + ")");
        }
        auto const extent = component.getExtent();
        validateExtent(name + "/" + componentName, extent, expectedExtent);
        auto chunk = component.loadChunk<T>();
        series.flush();
        return std::vector<T>(chunk.get(), chunk.get() + elementCount(extent));
    }

    void validateComputeSettings(io::Iteration const& iteration)
    {
        if(iteration.containsAttribute(field::writeVtk) && attribute<bool>(iteration, field::writeVtk))
        {
            validationError(
                "compute/" + std::string{field::writeVtk},
                "unsupported compute setting; openPMD parser rejects VTK output requests");
        }
        if(iteration.containsAttribute(field::devices))
        {
            validationError(
                "compute/" + std::string{field::devices},
                "unsupported compute setting; explicit device lists are not preserved by this transport");
        }
    }

    template<typename T>
    void validateUnchangedAttribute(io::Iteration const& iteration, std::string const& name, T const& expected)
    {
        if(!iteration.containsAttribute(name))
        {
            return;
        }

        auto const actual = attribute<T>(iteration, name);
        if(actual != expected)
        {
            validationError(
                "dynamic iteration/" + name,
                "non-dynamic attribute changed after iteration 0; dynamic-only iterations may update only "
                "core_beta_volume and core_point_beta");
        }
    }

    bool isAllowedDynamicMesh(std::string const& name, std::string const& prefix)
    {
        return name == prefix + "beta_volume" || name == prefix + "point_beta";
    }

    std::string entityFromAxes(std::vector<std::string> const& axes)
    {
        std::string entity;
        for(auto const& axis : axes)
        {
            if(!entity.empty())
            {
                entity += "_";
            }
            entity += axis;
        }
        return entity;
    }

    std::vector<unsigned long long> shapeFromExtent(io::Extent const& extent)
    {
        return std::vector<unsigned long long>(extent.begin(), extent.end());
    }

    template<typename T>
    std::vector<T> concatenate(std::vector<T> first, std::vector<T> const& second)
    {
        first.insert(first.end(), second.begin(), second.end());
        return first;
    }

    struct DerivedTopology
    {
        std::vector<unsigned> trianglePointIndices;
        std::vector<double> points;
        std::vector<double> triangleCenterX;
        std::vector<double> triangleCenterY;
        std::vector<unsigned> triangleNormalPoint;
        std::vector<double> triangleNormalsX;
        std::vector<double> triangleNormalsY;
        std::vector<int> forbiddenEdge;
        std::vector<int> triangleNeighbors;
        std::vector<float> triangleSurfaces;
    };

    [[nodiscard]] unsigned triangleVertex(
        std::vector<unsigned> const& trianglePointIndices,
        unsigned numberOfCells,
        unsigned triangle,
        unsigned localVertex)
    {
        return trianglePointIndices.at(localVertex * numberOfCells + triangle);
    }

    DerivedTopology deriveBackendTopology(
        std::vector<double> points,
        std::vector<unsigned> trianglePointIndices,
        unsigned numberOfCells,
        unsigned numberOfPoints)
    {
        DerivedTopology topology;
        topology.points = std::move(points);
        topology.trianglePointIndices = std::move(trianglePointIndices);
        topology.triangleCenterX.resize(numberOfCells);
        topology.triangleCenterY.resize(numberOfCells);
        topology.triangleNormalPoint.resize(3u * numberOfCells);
        topology.triangleNormalsX.resize(3u * numberOfCells);
        topology.triangleNormalsY.resize(3u * numberOfCells);
        topology.forbiddenEdge.assign(3u * numberOfCells, -1);
        topology.triangleNeighbors.assign(3u * numberOfCells, -1);
        topology.triangleSurfaces.resize(numberOfCells);

        using EdgeOwner = std::pair<unsigned, unsigned>;
        std::map<std::pair<unsigned, unsigned>, EdgeOwner> edgeOwners;
        for(unsigned triangle = 0; triangle < numberOfCells; ++triangle)
        {
            auto const p0 = triangleVertex(topology.trianglePointIndices, numberOfCells, triangle, 0u);
            auto const p1 = triangleVertex(topology.trianglePointIndices, numberOfCells, triangle, 1u);
            auto const p2 = triangleVertex(topology.trianglePointIndices, numberOfCells, triangle, 2u);
            if(p0 >= numberOfPoints || p1 >= numberOfPoints || p2 >= numberOfPoints)
            {
                validationError(
                    "canonical topology",
                    "triangle connectivity references a point outside the base level");
            }

            auto const x0 = topology.points.at(p0);
            auto const y0 = topology.points.at(numberOfPoints + p0);
            auto const x1 = topology.points.at(p1);
            auto const y1 = topology.points.at(numberOfPoints + p1);
            auto const x2 = topology.points.at(p2);
            auto const y2 = topology.points.at(numberOfPoints + p2);
            topology.triangleCenterX.at(triangle) = (x0 + x1 + x2) / 3.0;
            topology.triangleCenterY.at(triangle) = (y0 + y1 + y2) / 3.0;
            topology.triangleSurfaces.at(triangle)
                = static_cast<float>(0.5 * std::abs((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)));

            std::array<unsigned, 3> vertices{p0, p1, p2};
            for(unsigned edge = 0; edge < 3u; ++edge)
            {
                auto const a = vertices.at(edge);
                auto const b = vertices.at((edge + 1u) % 3u);
                auto const edgeX = topology.points.at(b) - topology.points.at(a);
                auto const edgeY = topology.points.at(numberOfPoints + b) - topology.points.at(numberOfPoints + a);
                auto const norm = std::sqrt(edgeX * edgeX + edgeY * edgeY);
                if(norm <= 0.0)
                {
                    validationError("canonical topology", "triangle contains a zero-length edge");
                }
                auto const flat = edge * numberOfCells + triangle;
                topology.triangleNormalsX.at(flat) = edgeY / norm;
                topology.triangleNormalsY.at(flat) = -edgeX / norm;
                topology.triangleNormalPoint.at(flat) = a;

                auto key = std::minmax(a, b);
                auto [owner, inserted] = edgeOwners.emplace(key, EdgeOwner{triangle, edge});
                if(!inserted)
                {
                    auto const [otherTriangle, otherEdge] = owner->second;
                    topology.triangleNeighbors.at(flat) = static_cast<int>(otherTriangle);
                    topology.forbiddenEdge.at(flat) = static_cast<int>(otherEdge);
                    auto const otherFlat = otherEdge * numberOfCells + otherTriangle;
                    topology.triangleNeighbors.at(otherFlat) = static_cast<int>(triangle);
                    topology.forbiddenEdge.at(otherFlat) = static_cast<int>(edge);
                }
            }
        }
        return topology;
    }

    DerivedTopology loadLegacyTopology(
        io::Series& series,
        io::Iteration& iteration,
        std::string const& prefix,
        unsigned numberOfPoints,
        unsigned numberOfCells)
    {
        auto vertices = iteration.meshes[prefix + "vertices"];
        validateHaseMetadata(
            prefix + "vertices",
            vertices,
            {"coordinate", "point"},
            {2u, numberOfPoints},
            false,
            true,
            "m");
        validateAxisLabels(prefix + "vertices", vertices, {"point"});
        auto points = concatenate(
            loadComponent<double>(series, iteration, prefix + "vertices", "x", io::Extent{numberOfPoints}),
            loadComponent<double>(series, iteration, prefix + "vertices", "y", io::Extent{numberOfPoints}));

        auto trianglePointIndices = loadScalar<unsigned>(
            series,
            iteration,
            prefix + "connectivity",
            io::Extent{3u * numberOfCells},
            {"cell", "local_vertex"},
            {numberOfCells, 3u},
            false,
            true);
        return deriveBackendTopology(
            std::move(points),
            std::move(trianglePointIndices),
            numberOfCells,
            numberOfPoints);
    }

    DerivedTopology loadCanonicalTopology(
        io::Series& series,
        io::Iteration& iteration,
        std::string const& prefix,
        unsigned numberOfPoints,
        unsigned numberOfCells,
        unsigned numberOfLevels)
    {
        auto const numberOfPrisms = numberOfCells * (numberOfLevels - 1u);
        auto const numberOfMeshPoints = numberOfPoints * numberOfLevels;
        auto x = loadComponent<double>(series, iteration, prefix + "points", "x", io::Extent{numberOfMeshPoints});
        auto y = loadComponent<double>(series, iteration, prefix + "points", "y", io::Extent{numberOfMeshPoints});
        (void) loadComponent<double>(series, iteration, prefix + "points", "z", io::Extent{numberOfMeshPoints});
        auto connectivity = loadScalar<unsigned>(
            series,
            iteration,
            prefix + "cells_connectivity",
            io::Extent{6u * numberOfPrisms},
            {"cell", "local_vertex"},
            {numberOfPrisms, 6u},
            false,
            false);
        auto offsets = loadScalar<unsigned>(
            series,
            iteration,
            prefix + "cells_offsets",
            io::Extent{numberOfPrisms + 1u},
            {"cell_offset"},
            {numberOfPrisms + 1u},
            false,
            false);
        auto cellTypes = loadScalar<unsigned>(
            series,
            iteration,
            prefix + "cells_types",
            io::Extent{numberOfPrisms},
            {"cell"},
            {numberOfPrisms},
            false,
            false);

        for(unsigned prism = 0; prism < numberOfPrisms; ++prism)
        {
            if(offsets.at(prism) != 6u * prism || cellTypes.at(prism) != 13u)
            {
                validationError(
                    "canonical topology",
                    "only contiguous VTK_WEDGE cells are supported by the current backend");
            }
        }
        if(offsets.back() != 6u * numberOfPrisms)
        {
            validationError("canonical topology", "cell offsets do not match six-node wedge connectivity");
        }

        std::vector<double> basePoints;
        basePoints.reserve(2u * numberOfPoints);
        basePoints.insert(basePoints.end(), x.begin(), x.begin() + numberOfPoints);
        basePoints.insert(basePoints.end(), y.begin(), y.begin() + numberOfPoints);

        std::vector<unsigned> trianglePointIndices(3u * numberOfCells);
        for(unsigned local = 0; local < 3u; ++local)
        {
            for(unsigned triangle = 0; triangle < numberOfCells; ++triangle)
            {
                auto const point = connectivity.at(6u * triangle + local);
                if(point >= numberOfPoints)
                {
                    validationError(
                        "canonical topology",
                        "first-layer wedge connectivity must reference base-level points");
                }
                trianglePointIndices.at(local * numberOfCells + triangle) = point;
            }
        }
        return deriveBackendTopology(
            std::move(basePoints),
            std::move(trianglePointIndices),
            numberOfCells,
            numberOfPoints);
    }

    void initializeResultForMesh(hase::core::Result& result, hase::core::HostMesh const& mesh)
    {
        auto const numberOfSamples = mesh.numberOfPoints * mesh.numberOfLevels;
        result = hase::core::Result(
            std::vector<float>(numberOfSamples, 0.0f),
            std::vector<double>(numberOfSamples, 100000.0),
            std::vector<unsigned>(numberOfSamples, 0u),
            std::vector<double>(numberOfSamples, 0.0));
    }

    template<typename T>
    void writeScalar(
        io::Iteration& iteration,
        std::string const& name,
        std::vector<T>& values,
        io::Extent const& extent,
        std::vector<std::string> const& axisLabels,
        std::string const& unit = "1",
        double unitSI = 1.0,
        std::array<double, 7> unitDimension = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0})
    {
        auto record = iteration.meshes[name];
        record.setAttribute("geometry", "other");
        record.setAttribute("geometryParameters", "entity=point_level");
        record.setAttribute("dataOrder", "C");
        record.setAxisLabels(axisLabels);
        record.setAttribute("haseTransportVersion", std::string{HASE_TRANSPORT_VERSION});
        record.setAttribute("haseEntity", entityFromAxes(axisLabels));
        record.setAttribute("haseAxes", axisLabels);
        record.setAttribute("haseLayoutOrder", std::string{"recordC"});
        record.setAttribute("hasePrimitiveShape", shapeFromExtent(extent));
        record.setAttribute("haseStatic", false);
        record.setAttribute("haseDynamic", true);
        record.setAttribute("haseBackendRequired", false);
        record.setAttribute("haseUnit", unit);
        record.setGridSpacing(std::vector<double>(extent.size(), 1.0));
        record.setGridGlobalOffset(std::vector<double>(extent.size(), 0.0));
        record.setGridUnitSI(1.0);
        record.setUnitDimension(unitDimension);

        auto& component = record[io::MeshRecordComponent::SCALAR];
        component.setUnitSI(unitSI);
        component.setPosition(std::vector<double>(extent.size(), 0.0));
        component.resetDataset({io::determineDatatype<T>(), extent});
        component.storeChunk(values, io::Offset(extent.size(), 0u), extent);
    }
} // namespace

namespace hase::openpmd
{

    Parser::Parser(std::filesystem::path inputPath, std::filesystem::path outputPath)
        : m_inputPath(std::move(inputPath))
        , m_outputPath(std::move(outputPath))
    {
    }

#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
    Parser::Parser(std::filesystem::path inputPath, std::filesystem::path outputPath, MPI_Comm comm)
        : m_inputPath(std::move(inputPath))
        , m_outputPath(std::move(outputPath))
        , m_comm(comm)
    {
    }
#endif

    bool Parser::isHeadRank() const
    {
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
        int rank = 0;
        MPI_Comm_rank(m_comm, &rank);
        return rank == 0;
#else
        return true;
#endif
    }

    core::SimulationContext Parser::read()
    {
        auto const stream = m_inputPath.string();
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
        io::Series series(stream, io::Access::READ_LINEAR, m_comm, seriesConfig(stream));
#else
        io::Series series(stream, io::Access::READ_LINEAR, seriesConfig(stream));
#endif

        for(auto& [index, iteration] : series.snapshots())
        {
            (void) index;
            auto context = readIteration(series, iteration);
            series.close();
            return context;
        }

        throw std::runtime_error("No iteration was available in the openPMD input stream.");
    }

    bool Parser::containsStaticMeshUpdate(io::Iteration const& iteration) const
    {
        std::string const prefix = m_meshGroup + "_";
        return iteration.meshes.contains(prefix + "points") || iteration.meshes.contains(prefix + "cells_connectivity")
               || iteration.meshes.contains(prefix + "cells_offsets")
               || iteration.meshes.contains(prefix + "cells_types") || iteration.meshes.contains(prefix + "vertices")
               || iteration.meshes.contains(prefix + "connectivity");
    }

    bool Parser::hasStaticMeshUpdate(io::Iteration const& iteration) const
    {
        std::string const prefix = m_meshGroup + "_";
        bool const hasCanonical = iteration.meshes.contains(prefix + "points")
                                  || iteration.meshes.contains(prefix + "cells_connectivity")
                                  || iteration.meshes.contains(prefix + "cells_offsets")
                                  || iteration.meshes.contains(prefix + "cells_types");
        if(hasCanonical)
        {
            bool const completeCanonical = iteration.meshes.contains(prefix + "points")
                                           && iteration.meshes.contains(prefix + "cells_connectivity")
                                           && iteration.meshes.contains(prefix + "cells_offsets")
                                           && iteration.meshes.contains(prefix + "cells_types");
            if(!completeCanonical)
            {
                validationError("canonical topology", "partial static topology update");
            }
            return true;
        }

        bool const hasLegacy
            = iteration.meshes.contains(prefix + "vertices") || iteration.meshes.contains(prefix + "connectivity");
        if(hasLegacy
           && !(iteration.meshes.contains(prefix + "vertices") && iteration.meshes.contains(prefix + "connectivity")))
        {
            validationError("legacy topology", "partial static topology update");
        }
        return hasLegacy;
    }

    void Parser::validateDynamicOnlyIteration(
        io::Iteration const& iteration,
        core::SimulationContext const& simulation) const
    {
        std::string const prefix = m_meshGroup + "_";
        for(auto const& [name, record] : iteration.meshes)
        {
            (void) record;
            if(!isAllowedDynamicMesh(name, prefix))
            {
                validationError(
                    name,
                    "non-dynamic mesh record present after iteration 0; dynamic-only iterations may update only "
                    "core_beta_volume and core_point_beta");
            }
        }

        if(iteration.containsAttribute("haseStaticUpdate") && attribute<bool>(iteration, "haseStaticUpdate"))
        {
            validationError(
                "dynamic iteration/haseStaticUpdate",
                "static updates after iteration 0 are not supported; dynamic-only iterations may update only "
                "core_beta_volume and core_point_beta");
        }

        validateComputeSettings(iteration);

        validateUnchangedAttribute<unsigned>(iteration, field::numberOfPoints, simulation.mesh.numberOfPoints);
        validateUnchangedAttribute<unsigned>(iteration, field::numberOfCells, simulation.mesh.numberOfTriangles);
        validateUnchangedAttribute<unsigned>(iteration, field::numberOfLevels, simulation.mesh.numberOfLevels);
        validateUnchangedAttribute<float>(iteration, field::thickness, simulation.mesh.thickness);
        validateUnchangedAttribute<float>(iteration, field::nTot, simulation.mesh.nTot);
        validateUnchangedAttribute<float>(iteration, field::crystalTFluo, simulation.mesh.crystalTFluo);
        validateUnchangedAttribute<unsigned>(iteration, field::claddingNumber, simulation.mesh.claddingNumber);
        validateUnchangedAttribute<double>(iteration, field::claddingAbsorption, simulation.mesh.claddingAbsorption);

        validateUnchangedAttribute<unsigned>(
            iteration,
            field::minRaysPerSample,
            simulation.experiment.minRaysPerSample);
        validateUnchangedAttribute<unsigned>(
            iteration,
            field::maxRaysPerSample,
            simulation.experiment.maxRaysPerSample);
        validateUnchangedAttribute<double>(iteration, field::mseThreshold, simulation.experiment.mseThreshold);
        validateUnchangedAttribute<bool>(iteration, field::useReflections, simulation.experiment.useReflections);
        validateUnchangedAttribute<unsigned>(iteration, field::spectralResolution, simulation.experiment.spectral);
        validateUnchangedAttribute<bool>(iteration, field::monochromatic, simulation.experiment.monochromatic);
        validateUnchangedAttribute<double>(iteration, field::maxSigmaAbsorption, simulation.experiment.maxSigmaA);
        validateUnchangedAttribute<double>(iteration, field::maxSigmaEmission, simulation.experiment.maxSigmaE);

        validateUnchangedAttribute<unsigned>(iteration, field::repetitions, simulation.compute.maxRepetitions);
        validateUnchangedAttribute<unsigned>(iteration, field::adaptiveSteps, simulation.compute.adaptiveSteps);
        validateUnchangedAttribute<unsigned>(iteration, field::maxGpus, simulation.compute.numDevices);
        validateUnchangedAttribute<std::string>(iteration, field::backend, simulation.compute.backend);
        validateUnchangedAttribute<std::string>(iteration, field::parallelMode, simulation.compute.parallelMode);
        validateUnchangedAttribute<unsigned>(iteration, field::minSampleRange, simulation.compute.minSampleRange);
        validateUnchangedAttribute<unsigned>(iteration, field::maxSampleRange, simulation.compute.maxSampleRange);
        validateUnchangedAttribute<unsigned>(iteration, field::rngSeed, simulation.compute.rngSeed);
    }

    core::SimulationContext Parser::readIteration(io::Series& series, io::Iteration& iteration)
    {
        std::string const prefix = m_meshGroup + "_";

        auto const numberOfPoints = attribute<unsigned>(iteration, field::numberOfPoints);
        auto const numberOfCells = attribute<unsigned>(iteration, field::numberOfCells);
        auto const numberOfLevels = attribute<unsigned>(iteration, field::numberOfLevels);

        auto topology
            = iteration.meshes.contains(prefix + "points")
                  ? loadCanonicalTopology(series, iteration, prefix, numberOfPoints, numberOfCells, numberOfLevels)
                  : loadLegacyTopology(series, iteration, prefix, numberOfPoints, numberOfCells);
        auto betaVolume = loadScalar<double>(
            series,
            iteration,
            prefix + "beta_volume",
            io::Extent{numberOfCells * (numberOfLevels - 1u)},
            {"cell", "layer"},
            {numberOfCells, numberOfLevels - 1u},
            true,
            true);
        auto betaCells = loadScalar<double>(
            series,
            iteration,
            prefix + "point_beta",
            io::Extent{numberOfPoints * numberOfLevels},
            {"point", "level"},
            {numberOfPoints, numberOfLevels},
            true,
            true);

        core::HostMesh mesh(
            std::move(topology.trianglePointIndices),
            numberOfCells,
            numberOfLevels,
            numberOfPoints,
            attribute<float>(iteration, field::thickness),
            std::move(topology.points),
            std::move(topology.triangleCenterX),
            std::move(topology.triangleCenterY),
            std::move(topology.triangleNormalPoint),
            std::move(topology.triangleNormalsX),
            std::move(topology.triangleNormalsY),
            std::move(topology.forbiddenEdge),
            std::move(topology.triangleNeighbors),
            std::move(topology.triangleSurfaces),
            std::move(betaVolume),
            std::move(betaCells),
            loadScalar<unsigned>(
                series,
                iteration,
                prefix + "cladding_cell_type",
                io::Extent{numberOfCells},
                {"cell"},
                {numberOfCells},
                false,
                true),
            loadScalar<float>(
                series,
                iteration,
                prefix + "refractive_index",
                io::Extent{4u},
                {"interface"},
                {4u},
                false,
                true),
            loadScalar<float>(
                series,
                iteration,
                prefix + "reflectivity",
                io::Extent{2u * numberOfCells},
                {"cell", "interface"},
                {numberOfCells, 2u},
                false,
                true),
            attribute<float>(iteration, field::nTot),
            attribute<float>(iteration, field::crystalTFluo),
            attribute<unsigned>(iteration, field::claddingNumber),
            attribute<double>(iteration, field::claddingAbsorption));

        core::ExperimentParameters experiment(
            attribute<unsigned>(iteration, field::minRaysPerSample),
            attribute<unsigned>(iteration, field::maxRaysPerSample),
            loadScalar<double>(
                series,
                iteration,
                prefix + "lambda_absorption",
                io::Extent{attribute<unsigned>(iteration, field::spectralResolution)},
                {"wavelength"},
                {attribute<unsigned>(iteration, field::spectralResolution)},
                false,
                false,
                "m"),
            loadScalar<double>(
                series,
                iteration,
                prefix + "lambda_emission",
                io::Extent{attribute<unsigned>(iteration, field::spectralResolution)},
                {"wavelength"},
                {attribute<unsigned>(iteration, field::spectralResolution)},
                false,
                false,
                "m"),
            loadScalar<double>(
                series,
                iteration,
                prefix + "sigma_absorption",
                io::Extent{attribute<unsigned>(iteration, field::spectralResolution)},
                {"wavelength"},
                {attribute<unsigned>(iteration, field::spectralResolution)},
                false,
                false,
                "cm^2"),
            loadScalar<double>(
                series,
                iteration,
                prefix + "sigma_emission",
                io::Extent{attribute<unsigned>(iteration, field::spectralResolution)},
                {"wavelength"},
                {attribute<unsigned>(iteration, field::spectralResolution)},
                false,
                false,
                "cm^2"),
            0.0,
            0.0,
            attribute<double>(iteration, field::mseThreshold),
            attribute<bool>(iteration, field::useReflections),
            attribute<unsigned>(iteration, field::spectralResolution),
            attributeOr<bool>(iteration, field::monochromatic, false));

        experiment.maxSigmaA = attributeOr<double>(iteration, field::maxSigmaAbsorption, experiment.maxSigmaA);
        experiment.maxSigmaE = attributeOr<double>(iteration, field::maxSigmaEmission, experiment.maxSigmaE);

        validateComputeSettings(iteration);

        auto const numberOfSamples = numberOfPoints * numberOfLevels;
        core::ComputeParameters compute(
            attribute<unsigned>(iteration, field::repetitions),
            attribute<unsigned>(iteration, field::adaptiveSteps),
            attribute<unsigned>(iteration, field::maxGpus),
            0u,
            attribute<std::string>(iteration, field::backend),
            attributeOr<std::string>(iteration, field::parallelMode, core::ParallelMode::SINGLE),
            false,
            std::vector<unsigned>{},
            attributeOr<unsigned>(iteration, field::minSampleRange, 0u),
            attributeOr<unsigned>(iteration, field::maxSampleRange, numberOfSamples - 1u),
            attributeOr<unsigned>(iteration, field::rngSeed, core::ComputeParameters::unspecifiedRngSeed));

        mesh.calcTotalReflectionAngles();

        core::Result result;
        initializeResultForMesh(result, mesh);
        iteration.close();
        return {std::move(experiment), std::move(compute), std::move(mesh), std::move(result)};
    }

    void Parser::updateDynamicIteration(
        io::Series& series,
        io::Iteration& iteration,
        core::SimulationContext& simulation)
    {
        std::string const prefix = m_meshGroup + "_";
        auto const numberOfCells = simulation.mesh.numberOfTriangles;
        auto const numberOfPoints = simulation.mesh.numberOfPoints;
        auto const numberOfLevels = simulation.mesh.numberOfLevels;
        simulation.mesh.betaVolume = loadScalar<double>(
            series,
            iteration,
            prefix + "beta_volume",
            io::Extent{numberOfCells * (numberOfLevels - 1u)},
            {"cell", "layer"},
            {numberOfCells, numberOfLevels - 1u},
            true,
            true);
        simulation.mesh.betaCells = loadScalar<double>(
            series,
            iteration,
            prefix + "point_beta",
            io::Extent{numberOfPoints * numberOfLevels},
            {"point", "level"},
            {numberOfPoints, numberOfLevels},
            true,
            true);
        initializeResultForMesh(simulation.result, simulation.mesh);
        iteration.close();
    }

    void Parser::writeResult(core::Result const& result, core::HostMesh const& mesh)
    {
        if(!isHeadRank())
        {
            return;
        }

        auto const stream = m_outputPath.string();
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
        // Only the head rank writes result streams.  Use a self communicator so
        // openPMD backends do not wait for ranks that intentionally skipped the
        // output path.
        io::Series series(stream, io::Access::CREATE_LINEAR, MPI_COMM_SELF, seriesConfig(stream));
#else
        io::Series series(stream, io::Access::CREATE_LINEAR, seriesConfig(stream));
#endif
        series.setAttribute("haseTransportVersion", std::string{HASE_TRANSPORT_VERSION});
        writeResultIteration(series, 0u, result, mesh);
        series.close();
    }

    void Parser::writeResultIteration(
        io::Series& series,
        std::uint64_t iterationIndex,
        core::Result const& result,
        core::HostMesh const& mesh)
    {
        auto const extent = io::Extent{mesh.numberOfPoints, mesh.numberOfLevels};
        auto iterations = series.writeIterations();
        auto iteration = iterations[iterationIndex];
        iteration.setTime(0.0);
        iteration.setDt(1.0);
        iteration.setTimeUnitSI(1.0);
        iteration.setAttribute(field::numberOfPoints, mesh.numberOfPoints);
        iteration.setAttribute(field::numberOfLevels, mesh.numberOfLevels);

        std::string const prefix = m_meshGroup + "_result_";
        auto phiAse = result.phiAse;
        auto mse = result.mse;
        auto totalRays = result.totalRays;
        auto dndtAse = result.dndtAse;
        writeScalar(
            iteration,
            prefix + "phi_ase",
            phiAse,
            extent,
            {"point", "level"},
            "cm^-2 s^-1",
            1.0e4,
            {-2.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0});
        writeScalar(iteration, prefix + "mse", mse, extent, {"point", "level"});
        writeScalar(iteration, prefix + "total_rays", totalRays, extent, {"point", "level"}, "count");
        writeScalar(
            iteration,
            prefix + "dndt_ase",
            dndtAse,
            extent,
            {"point", "level"},
            "s^-1",
            1.0,
            {0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0});

        iteration.close();
    }

    void Parser::processAll(std::function<void(core::SimulationContext&)> process)
    {
        auto const inputStream = m_inputPath.string();
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
        io::Series inputSeries(inputStream, io::Access::READ_LINEAR, m_comm, seriesConfig(inputStream));
#else
        io::Series inputSeries(inputStream, io::Access::READ_LINEAR, seriesConfig(inputStream));
#endif

        std::unique_ptr<io::Series> outputSeries;
        if(isHeadRank())
        {
            auto const outputStream = m_outputPath.string();
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
            // Only rank 0 owns the output Series; keep the writer communicator
            // local to that rank instead of passing the full input communicator.
            outputSeries = std::make_unique<io::Series>(
                outputStream,
                io::Access::CREATE_LINEAR,
                MPI_COMM_SELF,
                seriesConfig(outputStream));
#else
            outputSeries
                = std::make_unique<io::Series>(outputStream, io::Access::CREATE_LINEAR, seriesConfig(outputStream));
#endif
            outputSeries->setAttribute("haseTransportVersion", std::string{HASE_TRANSPORT_VERSION});
        }

        std::optional<core::SimulationContext> simulation;
        for(auto iteration : inputSeries.readIterations())
        {
            auto const iterationIndex = iteration.iterationIndex;
            if(!simulation)
            {
                if(!hasStaticMeshUpdate(iteration))
                {
                    validationError("dynamic iteration", "arrived before any static mesh update");
                }
                simulation = readIteration(inputSeries, iteration);
            }
            else
            {
                if(containsStaticMeshUpdate(iteration))
                {
                    validationError(
                        "dynamic iteration",
                        "static topology updates after iteration 0 are not supported; dynamic-only iterations may "
                        "update only core_beta_volume and core_point_beta");
                }
                validateDynamicOnlyIteration(iteration, *simulation);
                updateDynamicIteration(inputSeries, iteration, *simulation);
            }

            process(*simulation);
            if(outputSeries)
            {
                writeResultIteration(*outputSeries, iterationIndex, simulation->result, simulation->mesh);
            }
        }

        inputSeries.close();
        if(outputSeries)
        {
            outputSeries->close();
        }
    }

} // namespace hase::openpmd
