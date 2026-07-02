#include <catch2/catch_approx.hpp>
#ifdef HASE_OPENPMD_PARSER_VALIDATION_CUSTOM_MAIN
#    include <catch2/catch_session.hpp>
#endif
#include <catch2/catch_test_macros.hpp>
#include <openPMD/openPMD.hpp>
#include <openpmd/OpenPmdParser.hpp>

#include <array>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace io = openPMD;

namespace
{
#ifndef HASE_OPENPMD_FILE_EXTENSION
#    define HASE_OPENPMD_FILE_EXTENSION "bp"
#endif
#ifndef HASE_OPENPMD_TEST_FILE_EXTENSION
#    define HASE_OPENPMD_TEST_FILE_EXTENSION HASE_OPENPMD_FILE_EXTENSION
#endif

    constexpr char const* HASE_TRANSPORT_VERSION = "0.1";

    std::filesystem::path testPath(std::string const& name)
    {
        auto path = std::filesystem::temp_directory_path()
                    / ("hase_openpmd_" + name + "." + HASE_OPENPMD_TEST_FILE_EXTENSION);
        std::filesystem::remove_all(path);
        return path;
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

    void requireNear(std::vector<double> const& actual, std::vector<double> const& expected)
    {
        REQUIRE(actual.size() == expected.size());
        for(std::size_t i = 0; i < actual.size(); ++i)
        {
            CHECK(actual[i] == Catch::Approx(expected[i]).epsilon(1e-14));
        }
    }

    void setMetadata(
        io::Mesh& record,
        std::vector<std::string> const& axes,
        std::vector<unsigned long long> const& primitiveShape,
        bool dynamic = false,
        bool backendRequired = true,
        std::string const& unit = "1",
        std::vector<std::string> const& axisLabels = {"flatIndex"})
    {
        record.setAttribute("geometry", "other");
        record.setAttribute("geometryParameters", "topology=unstructured_triangular_prism");
        record.setAttribute("dataOrder", "C");
        record.setAxisLabels(axisLabels);
        record.setAttribute("haseTransportVersion", std::string{HASE_TRANSPORT_VERSION});
        record.setAttribute("haseEntity", entityFromAxes(axes));
        record.setAttribute("haseAxes", axes);
        record.setAttribute("haseLayoutOrder", std::string{"backendFlat"});
        record.setAttribute("hasePrimitiveShape", primitiveShape);
        record.setAttribute("haseStatic", !dynamic);
        record.setAttribute("haseDynamic", dynamic);
        record.setAttribute("haseBackendRequired", backendRequired);
        record.setAttribute("haseUnit", unit);
        record.setGridSpacing(std::vector<double>(axisLabels.size(), 1.0));
        record.setGridGlobalOffset(std::vector<double>(axisLabels.size(), 0.0));
        record.setGridUnitSI(1.0);
        record.setUnitDimension(std::array<double, 7>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    }

    template<typename T>
    void writeScalar(
        io::Series& series,
        io::Iteration& iteration,
        std::string const& name,
        std::vector<T> values,
        std::vector<std::string> const& axes,
        std::vector<unsigned long long> const& primitiveShape,
        bool dynamic = false,
        bool backendRequired = true,
        std::string const& unit = "1")
    {
        auto record = iteration.meshes[name];
        setMetadata(record, axes, primitiveShape, dynamic, backendRequired, unit);
        auto& component = record[io::MeshRecordComponent::SCALAR];
        io::Extent extent{values.size()};
        component.setUnitSI(1.0);
        component.setPosition(std::vector<double>{0.0});
        component.resetDataset({io::determineDatatype<T>(), extent});
        component.storeChunk(values, io::Offset{0u}, extent);
        series.flush();
    }

    template<typename T>
    void writeComponent(
        io::Series& series,
        io::Iteration& iteration,
        std::string const& recordName,
        std::string const& componentName,
        std::vector<T> values,
        std::vector<std::string> const& axisLabels)
    {
        auto record = iteration.meshes[recordName];
        record.setAttribute("geometry", "other");
        record.setAttribute("dataOrder", "C");
        record.setAxisLabels(axisLabels);
        record.setGridSpacing(std::vector<double>(axisLabels.size(), 1.0));
        record.setGridGlobalOffset(std::vector<double>(axisLabels.size(), 0.0));
        record.setGridUnitSI(1.0);
        auto& component = record[componentName];
        io::Extent extent{values.size()};
        component.setUnitSI(1.0);
        component.setPosition(std::vector<double>{0.0});
        component.resetDataset({io::determineDatatype<T>(), extent});
        component.storeChunk(values, io::Offset{0u}, extent);
        series.flush();
    }

    std::filesystem::path writeParserInput(
        std::string const& name,
        std::function<void(io::Series&, io::Iteration&)> mutate = {},
        bool betaVolumeAsFloat = false,
        bool pointBetaTransposed = false,
        bool betaVolumeBadExtent = false)
    {
        auto path = testPath(name);
        io::Series series(path.string(), io::Access::CREATE_LINEAR, "{}");
        series.setAttribute("haseTransportVersion", std::string{HASE_TRANSPORT_VERSION});
        auto iteration = series.snapshots()[0];
        iteration.setTime(0.0);
        iteration.setDt(1.0);
        iteration.setTimeUnitSI(1.0);

        iteration.setAttribute("number_of_points", 3u);
        iteration.setAttribute("number_of_cells", 1u);
        iteration.setAttribute("number_of_levels", 2u);
        iteration.setAttribute("thickness", 0.25f);
        iteration.setAttribute("n_tot", 5.0f);
        iteration.setAttribute("crystal_t_fluo", 1.25f);
        iteration.setAttribute("cladding_number", 7u);
        iteration.setAttribute("cladding_absorption", 0.05);
        iteration.setAttribute("min_rays_per_sample", 1u);
        iteration.setAttribute("max_rays_per_sample", 2u);
        iteration.setAttribute("mse_threshold", 0.5);
        iteration.setAttribute("repetitions", 3u);
        iteration.setAttribute("adaptive_steps", 4u);
        iteration.setAttribute("max_gpus", 1u);
        iteration.setAttribute("backend", std::string{"Host_Cpu_CpuSerial"});
        iteration.setAttribute("parallel_mode", std::string{"single"});
        iteration.setAttribute("min_sample_range", 0u);
        iteration.setAttribute("max_sample_range", 5u);
        iteration.setAttribute("rng_seed", 1234u);
        iteration.setAttribute("use_reflections", true);
        iteration.setAttribute("spectral_resolution", 2u);
        iteration.setAttribute("monochromatic", false);
        iteration.setAttribute("max_sigma_absorption", 0.02);
        iteration.setAttribute("max_sigma_emission", 0.04);

        writeComponent<double>(series, iteration, "core_vertices", "x", {0.0, 1.0, 0.0}, {"point"});
        writeComponent<double>(series, iteration, "core_vertices", "y", {0.0, 0.0, 1.0}, {"point"});
        setMetadata(iteration.meshes["core_vertices"], {"coordinate", "point"}, {2u, 3u}, false, true, "m", {"point"});
        writeComponent<double>(series, iteration, "core_cell_center", "x", {1.0 / 3.0}, {"cell"});
        writeComponent<double>(series, iteration, "core_cell_center", "y", {1.0 / 3.0}, {"cell"});
        setMetadata(iteration.meshes["core_cell_center"], {"cell"}, {1u}, false, true, "m", {"cell"});

        writeScalar<unsigned>(
            series,
            iteration,
            "core_connectivity",
            {0u, 1u, 2u},
            {"cell", "local_vertex"},
            {1u, 3u});
        writeScalar<int>(series, iteration, "core_neighbors", {-1, -1, -1}, {"cell", "local_side"}, {1u, 3u});
        writeScalar<int>(series, iteration, "core_forbidden_edges", {-1, -1, -1}, {"cell", "local_side"}, {1u, 3u});
        writeScalar<unsigned>(series, iteration, "core_normal_points", {0u, 1u, 2u}, {"cell", "local_side"}, {1u, 3u});
        writeScalar<double>(
            series,
            iteration,
            "core_cell_normal_x",
            {0.0, 1.0, -1.0},
            {"cell", "local_side"},
            {1u, 3u});
        writeScalar<double>(
            series,
            iteration,
            "core_cell_normal_y",
            {-1.0, 1.0, 0.0},
            {"cell", "local_side"},
            {1u, 3u});
        writeScalar<float>(series, iteration, "core_surface", {0.5f}, {"cell"}, {1u}, false, true, "m^2");
        if(betaVolumeAsFloat)
        {
            writeScalar<float>(series, iteration, "core_beta_volume", {0.1f}, {"cell", "layer"}, {1u, 1u}, true);
        }
        else if(betaVolumeBadExtent)
        {
            writeScalar<double>(series, iteration, "core_beta_volume", {0.1, 0.2}, {"cell", "layer"}, {1u, 1u}, true);
        }
        else
        {
            writeScalar<double>(series, iteration, "core_beta_volume", {0.1}, {"cell", "layer"}, {1u, 1u}, true);
        }
        if(pointBetaTransposed)
        {
            auto record = iteration.meshes["core_point_beta"];
            setMetadata(record, {"point", "level"}, {3u, 2u}, true);
            auto& component = record[io::MeshRecordComponent::SCALAR];
            std::vector<double> values{0.1, 0.4, 0.2, 0.5, 0.3, 0.6};
            io::Extent extent{2u, 3u};
            component.setUnitSI(1.0);
            component.setPosition(std::vector<double>{0.0, 0.0});
            component.resetDataset({io::determineDatatype<double>(), extent});
            component.storeChunk(values, io::Offset{0u, 0u}, extent);
            series.flush();
        }
        else
        {
            writeScalar<double>(
                series,
                iteration,
                "core_point_beta",
                {0.1, 0.2, 0.3, 0.4, 0.5, 0.6},
                {"point", "level"},
                {3u, 2u},
                true);
        }
        writeScalar<unsigned>(series, iteration, "core_cladding_cell_type", {0u}, {"cell"}, {1u});
        writeScalar<float>(series, iteration, "core_refractive_index", {1.5f, 1.0f, 1.5f, 1.0f}, {"interface"}, {4u});
        writeScalar<float>(series, iteration, "core_reflectivity", {0.1f, 0.2f}, {"cell", "interface"}, {1u, 2u});
        writeScalar<double>(
            series,
            iteration,
            "core_lambda_absorption",
            {900.0, 910.0},
            {"wavelength"},
            {2u},
            false,
            false,
            "m");
        writeScalar<double>(
            series,
            iteration,
            "core_lambda_emission",
            {1000.0, 1010.0},
            {"wavelength"},
            {2u},
            false,
            false,
            "m");
        writeScalar<double>(
            series,
            iteration,
            "core_sigma_absorption",
            {0.01, 0.02},
            {"wavelength"},
            {2u},
            false,
            false,
            "cm^2");
        writeScalar<double>(
            series,
            iteration,
            "core_sigma_emission",
            {0.03, 0.04},
            {"wavelength"},
            {2u},
            false,
            false,
            "cm^2");

        if(mutate)
        {
            mutate(series, iteration);
        }

        iteration.close();
        series.close();
        return path;
    }

    std::filesystem::path writeCanonicalStaticDynamicInput(
        std::string const& name,
        std::function<void(io::Series&, io::Iteration&)> mutateSecond = {})
    {
        auto path = testPath(name);
        io::Series series(path.string(), io::Access::CREATE_LINEAR, "{}");
        series.setAttribute("haseTransportVersion", std::string{HASE_TRANSPORT_VERSION});

        auto setCommonAttributes = [](io::Iteration& iteration)
        {
            iteration.setTime(0.0);
            iteration.setDt(1.0);
            iteration.setTimeUnitSI(1.0);
            iteration.setAttribute("number_of_points", 3u);
            iteration.setAttribute("number_of_cells", 1u);
            iteration.setAttribute("number_of_levels", 2u);
            iteration.setAttribute("thickness", 0.25f);
            iteration.setAttribute("n_tot", 5.0f);
            iteration.setAttribute("crystal_t_fluo", 1.25f);
            iteration.setAttribute("cladding_number", 7u);
            iteration.setAttribute("cladding_absorption", 0.05);
            iteration.setAttribute("min_rays_per_sample", 1u);
            iteration.setAttribute("max_rays_per_sample", 2u);
            iteration.setAttribute("mse_threshold", 0.5);
            iteration.setAttribute("repetitions", 3u);
            iteration.setAttribute("adaptive_steps", 4u);
            iteration.setAttribute("max_gpus", 1u);
            iteration.setAttribute("backend", std::string{"Host_Cpu_CpuSerial"});
            iteration.setAttribute("parallel_mode", std::string{"single"});
            iteration.setAttribute("min_sample_range", 0u);
            iteration.setAttribute("max_sample_range", 5u);
            iteration.setAttribute("rng_seed", 1234u);
            iteration.setAttribute("use_reflections", true);
            iteration.setAttribute("spectral_resolution", 2u);
            iteration.setAttribute("monochromatic", false);
            iteration.setAttribute("max_sigma_absorption", 0.02);
            iteration.setAttribute("max_sigma_emission", 0.04);
        };

        auto first = series.snapshots()[0];
        setCommonAttributes(first);
        first.setAttribute("haseStaticUpdate", true);
        writeComponent<double>(series, first, "core_points", "x", {0.0, 1.0, 0.0, 0.0, 1.0, 0.0}, {"mesh_point"});
        writeComponent<double>(series, first, "core_points", "y", {0.0, 0.0, 1.0, 0.0, 0.0, 1.0}, {"mesh_point"});
        writeComponent<double>(series, first, "core_points", "z", {0.0, 0.0, 0.0, 0.25, 0.25, 0.25}, {"mesh_point"});
        setMetadata(
            first.meshes["core_points"],
            {"coordinate", "mesh_point"},
            {3u, 6u},
            false,
            false,
            "m",
            {"mesh_point"});
        writeScalar<unsigned>(
            series,
            first,
            "core_cells_connectivity",
            {0u, 1u, 2u, 3u, 4u, 5u},
            {"cell", "local_vertex"},
            {1u, 6u},
            false,
            false);
        writeScalar<unsigned>(series, first, "core_cells_offsets", {0u, 6u}, {"cell_offset"}, {2u}, false, false);
        writeScalar<unsigned>(series, first, "core_cells_types", {13u}, {"cell"}, {1u}, false, false);
        writeScalar<double>(series, first, "core_beta_volume", {0.1}, {"cell", "layer"}, {1u, 1u}, true);
        writeScalar<double>(
            series,
            first,
            "core_point_beta",
            {0.1, 0.2, 0.3, 0.4, 0.5, 0.6},
            {"point", "level"},
            {3u, 2u},
            true);
        writeScalar<unsigned>(series, first, "core_cladding_cell_type", {0u}, {"cell"}, {1u});
        writeScalar<float>(series, first, "core_refractive_index", {1.5f, 1.0f, 1.5f, 1.0f}, {"interface"}, {4u});
        writeScalar<float>(series, first, "core_reflectivity", {0.1f, 0.2f}, {"cell", "interface"}, {1u, 2u});
        writeScalar<
            double>(series, first, "core_lambda_absorption", {900.0, 910.0}, {"wavelength"}, {2u}, false, false, "m");
        writeScalar<
            double>(series, first, "core_lambda_emission", {1000.0, 1010.0}, {"wavelength"}, {2u}, false, false, "m");
        writeScalar<
            double>(series, first, "core_sigma_absorption", {0.01, 0.02}, {"wavelength"}, {2u}, false, false, "cm^2");
        writeScalar<
            double>(series, first, "core_sigma_emission", {0.03, 0.04}, {"wavelength"}, {2u}, false, false, "cm^2");
        first.close();

        auto second = series.snapshots()[1];
        second.setTime(1.0);
        second.setDt(1.0);
        second.setTimeUnitSI(1.0);
        second.setAttribute("haseStaticUpdate", false);
        writeScalar<double>(series, second, "core_beta_volume", {0.9}, {"cell", "layer"}, {1u, 1u}, true);
        writeScalar<double>(
            series,
            second,
            "core_point_beta",
            {1.1, 1.2, 1.3, 1.4, 1.5, 1.6},
            {"point", "level"},
            {3u, 2u},
            true);
        if(mutateSecond)
        {
            mutateSecond(series, second);
        }
        second.close();
        series.close();
        return path;
    }

    std::filesystem::path writePythonContractEquivalentInput(std::string const& name)
    {
        auto path = testPath(name);
        io::Series series(path.string(), io::Access::CREATE_LINEAR, "{}");
        series.setAttribute("haseTransportVersion", std::string{HASE_TRANSPORT_VERSION});
        auto iteration = series.snapshots()[0];
        iteration.setTime(0.0);
        iteration.setDt(1.0);
        iteration.setTimeUnitSI(1.0);

        constexpr unsigned numberOfPoints = 5u;
        constexpr unsigned numberOfCells = 3u;
        constexpr unsigned numberOfLevels = 6u;
        constexpr unsigned numberOfPrisms = numberOfCells * (numberOfLevels - 1u);
        std::array<double, numberOfPoints> const x{0.0, 1.5, 0.25, 2.25, -0.75};
        std::array<double, numberOfPoints> const y{0.0, -0.5, 1.25, 2.5, 3.75};
        std::array<std::array<unsigned, 3>, numberOfCells> const triangles{{{0u, 1u, 2u}, {2u, 3u, 4u}, {0u, 2u, 4u}}};

        iteration.setAttribute("number_of_points", numberOfPoints);
        iteration.setAttribute("number_of_cells", numberOfCells);
        iteration.setAttribute("number_of_levels", numberOfLevels);
        iteration.setAttribute("thickness", 0.375f);
        iteration.setAttribute("n_tot", 7.5f);
        iteration.setAttribute("crystal_t_fluo", 1.75f);
        iteration.setAttribute("cladding_number", 3u);
        iteration.setAttribute("cladding_absorption", 0.075);
        iteration.setAttribute("min_rays_per_sample", 1u);
        iteration.setAttribute("max_rays_per_sample", 1u);
        iteration.setAttribute("mse_threshold", 0.25);
        iteration.setAttribute("repetitions", 1u);
        iteration.setAttribute("adaptive_steps", 1u);
        iteration.setAttribute("max_gpus", 1u);
        iteration.setAttribute("backend", std::string{"Host_Cpu_CpuSerial"});
        iteration.setAttribute("parallel_mode", std::string{"single"});
        iteration.setAttribute("min_sample_range", 0u);
        iteration.setAttribute("max_sample_range", 0u);
        iteration.setAttribute("rng_seed", 1234u);
        iteration.setAttribute("use_reflections", true);
        iteration.setAttribute("spectral_resolution", 3u);
        iteration.setAttribute("monochromatic", false);
        iteration.setAttribute("max_sigma_absorption", 0.040);
        iteration.setAttribute("max_sigma_emission", 0.050);

        std::vector<double> pointsX;
        std::vector<double> pointsY;
        std::vector<double> pointsZ;
        pointsX.reserve(numberOfPoints * numberOfLevels);
        pointsY.reserve(numberOfPoints * numberOfLevels);
        pointsZ.reserve(numberOfPoints * numberOfLevels);
        for(unsigned level = 0u; level < numberOfLevels; ++level)
        {
            for(unsigned point = 0u; point < numberOfPoints; ++point)
            {
                pointsX.push_back(x.at(point));
                pointsY.push_back(y.at(point));
                pointsZ.push_back(0.375 * static_cast<double>(level));
            }
        }
        writeComponent<double>(series, iteration, "core_points", "x", pointsX, {"mesh_point"});
        writeComponent<double>(series, iteration, "core_points", "y", pointsY, {"mesh_point"});
        writeComponent<double>(series, iteration, "core_points", "z", pointsZ, {"mesh_point"});
        setMetadata(
            iteration.meshes["core_points"],
            {"coordinate", "mesh_point"},
            {3u, numberOfPoints * numberOfLevels},
            false,
            false,
            "m",
            {"mesh_point"});

        std::vector<unsigned> connectivity;
        connectivity.reserve(6u * numberOfPrisms);
        for(unsigned level = 0u; level + 1u < numberOfLevels; ++level)
        {
            auto const lower = level * numberOfPoints;
            auto const upper = (level + 1u) * numberOfPoints;
            for(auto const& triangle : triangles)
            {
                connectivity.push_back(triangle.at(0u) + lower);
                connectivity.push_back(triangle.at(1u) + lower);
                connectivity.push_back(triangle.at(2u) + lower);
                connectivity.push_back(triangle.at(0u) + upper);
                connectivity.push_back(triangle.at(1u) + upper);
                connectivity.push_back(triangle.at(2u) + upper);
            }
        }
        std::vector<unsigned> offsets(numberOfPrisms + 1u);
        for(unsigned prism = 0u; prism <= numberOfPrisms; ++prism)
        {
            offsets.at(prism) = 6u * prism;
        }
        writeScalar<unsigned>(
            series,
            iteration,
            "core_cells_connectivity",
            connectivity,
            {"cell", "local_vertex"},
            {numberOfPrisms, 6u},
            false,
            false);
        writeScalar<unsigned>(
            series,
            iteration,
            "core_cells_offsets",
            offsets,
            {"cell_offset"},
            {numberOfPrisms + 1u},
            false,
            false);
        writeScalar<unsigned>(
            series,
            iteration,
            "core_cells_types",
            std::vector<unsigned>(numberOfPrisms, 13u),
            {"cell"},
            {numberOfPrisms},
            false,
            false);

        writeScalar<double>(
            series,
            iteration,
            "core_beta_volume",
            {0.11, 0.21, 0.31, 0.12, 0.22, 0.32, 0.13, 0.23, 0.33, 0.14, 0.24, 0.34, 0.15, 0.25, 0.35},
            {"cell", "layer"},
            {numberOfCells, numberOfLevels - 1u},
            true);
        writeScalar<double>(
            series,
            iteration,
            "core_point_beta",
            {100.0, 110.0, 120.0, 130.0, 140.0, 101.0, 111.0, 121.0, 131.0, 141.0, 102.0, 112.0, 122.0, 132.0, 142.0,
             103.0, 113.0, 123.0, 133.0, 143.0, 104.0, 114.0, 124.0, 134.0, 144.0, 105.0, 115.0, 125.0, 135.0, 145.0},
            {"point", "level"},
            {numberOfPoints, numberOfLevels},
            true);
        writeScalar<unsigned>(series, iteration, "core_cladding_cell_type", {0u, 2u, 1u}, {"cell"}, {numberOfCells});
        writeScalar<float>(
            series,
            iteration,
            "core_refractive_index",
            {1.80f, 1.20f, 1.65f, 1.05f},
            {"interface"},
            {4u});
        writeScalar<float>(
            series,
            iteration,
            "core_reflectivity",
            {0.01f, 0.03f, 0.05f, 0.02f, 0.04f, 0.06f},
            {"cell", "interface"},
            {numberOfCells, 2u});
        writeScalar<double>(
            series,
            iteration,
            "core_lambda_absorption",
            {900e-9, 910e-9, 930e-9},
            {"wavelength"},
            {3u},
            false,
            false,
            "m");
        writeScalar<double>(
            series,
            iteration,
            "core_lambda_emission",
            {1000e-9, 1015e-9, 1040e-9},
            {"wavelength"},
            {3u},
            false,
            false,
            "m");
        writeScalar<double>(
            series,
            iteration,
            "core_sigma_absorption",
            {0.010, 0.025, 0.040},
            {"wavelength"},
            {3u},
            false,
            false,
            "cm^2");
        writeScalar<double>(
            series,
            iteration,
            "core_sigma_emission",
            {0.050, 0.035, 0.020},
            {"wavelength"},
            {3u},
            false,
            false,
            "cm^2");

        iteration.close();
        series.close();
        return path;
    }

    std::string parserError(std::filesystem::path const& path)
    {
        hase::openpmd::Parser parser{path, testPath("unused-output")};
        try
        {
            (void) parser.read();
        }
        catch(std::runtime_error const& err)
        {
            return err.what();
        }
        return {};
    }
} // namespace

TEST_CASE("openPMD parser reads a transport-valid openPMD record", "[openpmd][parser]")
{
    auto const path = writeParserInput("valid");
    hase::openpmd::Parser parser{path, testPath("valid-output")};
    auto context = parser.read();

    REQUIRE(context.mesh.numberOfPoints == 3u);
    REQUIRE(context.mesh.numberOfTriangles == 1u);
    REQUIRE(context.mesh.numberOfLevels == 2u);
    REQUIRE(context.mesh.trianglePointIndices == std::vector<unsigned>{0u, 1u, 2u});
    REQUIRE(context.mesh.betaCells.size() == 6u);
    REQUIRE(context.experiment.spectral == 2u);
    REQUIRE(context.compute.maxRepetitions == 3u);
    REQUIRE(context.compute.writeVtk == false);
    REQUIRE(context.compute.devices.empty());
    REQUIRE(context.compute.rngSeed == 1234u);
    REQUIRE(context.run.enableAse == true);
    REQUIRE(context.run.timeIntegration.method == "explicit-euler");
}

TEST_CASE("openPMD parser reads compiled simulation run-control attributes", "[openpmd][parser]")
{
    auto const path = writeParserInput(
        "run_control",
        [](io::Series& series, io::Iteration& iteration)
        {
            (void) series;
            iteration.setAttribute("time_step", 2.0e-5);
            iteration.setAttribute("number_of_steps", 100u);
            iteration.setAttribute("enable_ase", false);
            iteration.setAttribute("pump_steps", 50u);
            iteration.setAttribute("time_integrator", std::string{"frozen-phi-ase-runge-kutta-4"});
        });
    hase::openpmd::Parser parser{path, testPath("run-control-output")};
    auto context = parser.read();

    REQUIRE(context.run.timeStep == Catch::Approx(2.0e-5));
    REQUIRE(context.run.numberOfSteps == 100u);
    REQUIRE(context.run.enableAse == false);
    REQUIRE(context.run.pumpSteps == 50u);
    REQUIRE(context.run.timeIntegration.method == "frozen-phi-ase-runge-kutta-4");
}

TEST_CASE("openPMD parser rejects malformed fields before HostMesh construction", "[openpmd][parser]")
{
    SECTION("extent")
    {
        auto const path = writeParserInput("bad_extent", {}, false, false, true);
        auto const error = parserError(path);
        REQUIRE(error.find("openPMD validation error for 'core_beta_volume'") != std::string::npos);
        REQUIRE(error.find("extent mismatch") != std::string::npos);
    }

    SECTION("accidental transpose extent with matching element count")
    {
        auto const path = writeParserInput("bad_transpose_extent", {}, false, true);
        auto const error = parserError(path);
        REQUIRE(error.find("openPMD validation error for 'core_point_beta'") != std::string::npos);
        REQUIRE(error.find("extent mismatch") != std::string::npos);
    }

    SECTION("dtype")
    {
        auto const path = writeParserInput("bad_dtype", {}, true);
        auto const error = parserError(path);
        REQUIRE(error.find("openPMD validation error for 'core_beta_volume'") != std::string::npos);
        REQUIRE(error.find("dtype mismatch") != std::string::npos);
    }

    SECTION("metadata")
    {
        auto const path = writeParserInput(
            "bad_metadata",
            [](io::Series& series, io::Iteration& iteration)
            {
                (void) series;
                iteration.meshes["core_connectivity"].setAttribute("haseTransportVersion", std::string{"999"});
            });
        auto const error = parserError(path);
        REQUIRE(error.find("openPMD validation error for 'core_connectivity'") != std::string::npos);
        REQUIRE(error.find("haseTransportVersion") != std::string::npos);
    }

    SECTION("role")
    {
        auto const path = writeParserInput(
            "bad_role",
            [](io::Series& series, io::Iteration& iteration)
            {
                (void) series;
                iteration.meshes["core_beta_volume"].setAttribute("haseDynamic", false);
            });
        auto const error = parserError(path);
        REQUIRE(error.find("openPMD validation error for 'core_beta_volume'") != std::string::npos);
        REQUIRE(error.find("haseDynamic") != std::string::npos);
    }
}

TEST_CASE("openPMD parser rejects unsupported compute settings explicitly", "[openpmd][parser]")
{
    auto const path = writeParserInput(
        "bad_compute",
        [](io::Series& series, io::Iteration& iteration)
        {
            (void) series;
            iteration.setAttribute("write_vtk", true);
        });
    auto const error = parserError(path);
    REQUIRE(error.find("openPMD validation error for 'compute/write_vtk'") != std::string::npos);
    REQUIRE(error.find("unsupported compute setting") != std::string::npos);
}

TEST_CASE("openPMD parser processAll consumes stream until producer close", "[openpmd][parser]")
{
    auto const input = writeParserInput("process_all");
    auto const output = testPath("process_all_output");

    hase::openpmd::Parser parser{input, output};
    unsigned calls = 0u;
    parser.processAll(
        [&calls](hase::core::SimulationContext& context)
        {
            ++calls;
            context.result.phiAse = std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
        });

    REQUIRE(calls == 1u);

    io::Series series(output.string(), io::Access::READ_LINEAR, "{}");
    unsigned outputIterations = 0u;
    for(auto iteration : series.readIterations())
    {
        ++outputIterations;
        REQUIRE(iteration.iterationIndex == 0u);
        REQUIRE(iteration.getAttribute("number_of_points").get<unsigned>() == 3u);
        REQUIRE(iteration.getAttribute("number_of_levels").get<unsigned>() == 2u);

        auto component = iteration.meshes["core_result_phi_ase"][io::MeshRecordComponent::SCALAR];
        auto chunk = component.loadChunk<float>();
        series.flush();
        std::vector<float> phiAse(chunk.get(), chunk.get() + 6u);
        REQUIRE(phiAse == std::vector<float>{1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f});
        iteration.close();
    }
    series.close();

    REQUIRE(outputIterations == 1u);
}

TEST_CASE("openPMD parser processAll reuses cached topology for dynamic-only iterations", "[openpmd][parser]")
{
    auto const input = writeCanonicalStaticDynamicInput("process_all_dynamic");
    auto const output = testPath("process_all_dynamic_output");

    hase::openpmd::Parser parser{input, output};
    unsigned calls = 0u;
    parser.processAll(
        [&calls](hase::core::SimulationContext& context)
        {
            ++calls;
            REQUIRE(context.mesh.trianglePointIndices == std::vector<unsigned>{0u, 1u, 2u});
            if(calls == 1u)
            {
                REQUIRE(context.mesh.betaVolume == std::vector<double>{0.1});
                REQUIRE(context.mesh.betaCells == std::vector<double>{0.1, 0.2, 0.3, 0.4, 0.5, 0.6});
            }
            else
            {
                REQUIRE(context.mesh.betaVolume == std::vector<double>{0.9});
                REQUIRE(context.mesh.betaCells == std::vector<double>{1.1, 1.2, 1.3, 1.4, 1.5, 1.6});
            }
            context.result.phiAse = std::vector<float>(6u, static_cast<float>(calls));
        });

    REQUIRE(calls == 2u);

    io::Series series(output.string(), io::Access::READ_LINEAR, "{}");
    unsigned outputIterations = 0u;
    for(auto iteration : series.readIterations())
    {
        ++outputIterations;
        auto const iterationIndex = iteration.iterationIndex;
        REQUIRE(iterationIndex < 2u);

        auto component = iteration.meshes["core_result_phi_ase"][io::MeshRecordComponent::SCALAR];
        auto chunk = component.loadChunk<float>();
        series.flush();
        std::vector<float> phiAse(chunk.get(), chunk.get() + 6u);
        REQUIRE(phiAse == std::vector<float>(6u, static_cast<float>(iterationIndex + 1u)));
        iteration.close();
    }
    series.close();

    REQUIRE(outputIterations == 2u);
}

TEST_CASE("openPMD parser rejects non-dynamic changes after iteration 0", "[openpmd][parser]")
{
    auto const expectDynamicError = [](std::filesystem::path const& input)
    {
        hase::openpmd::Parser parser{input, testPath("dynamic_contract_output")};
        try
        {
            parser.processAll([](hase::core::SimulationContext&) {});
        }
        catch(std::runtime_error const& err)
        {
            return std::string{err.what()};
        }
        return std::string{};
    };

    SECTION("extra non-dynamic mesh record")
    {
        auto const input = writeCanonicalStaticDynamicInput(
            "dynamic_contract_static_record",
            [](io::Series& series, io::Iteration& iteration)
            { writeScalar<unsigned>(series, iteration, "core_cladding_cell_type", {0u}, {"cell"}, {1u}); });
        auto const error = expectDynamicError(input);
        REQUIRE(error.find("openPMD validation error for 'core_cladding_cell_type'") != std::string::npos);
        REQUIRE(error.find("non-dynamic mesh record") != std::string::npos);
    }

    SECTION("changed static attribute")
    {
        auto const input = writeCanonicalStaticDynamicInput(
            "dynamic_contract_static_attr",
            [](io::Series& series, io::Iteration& iteration)
            {
                (void) series;
                iteration.setAttribute("number_of_cells", 2u);
            });
        auto const error = expectDynamicError(input);
        REQUIRE(error.find("openPMD validation error for 'dynamic iteration/number_of_cells'") != std::string::npos);
        REQUIRE(error.find("non-dynamic attribute changed after iteration 0") != std::string::npos);
    }

    SECTION("changed spectral attribute")
    {
        auto const input = writeCanonicalStaticDynamicInput(
            "dynamic_contract_spectral_attr",
            [](io::Series& series, io::Iteration& iteration)
            {
                (void) series;
                iteration.setAttribute("spectral_resolution", 3u);
            });
        auto const error = expectDynamicError(input);
        REQUIRE(
            error.find("openPMD validation error for 'dynamic iteration/spectral_resolution'") != std::string::npos);
        REQUIRE(error.find("non-dynamic attribute changed after iteration 0") != std::string::npos);
    }

    SECTION("changed compute attribute")
    {
        auto const input = writeCanonicalStaticDynamicInput(
            "dynamic_contract_compute_attr",
            [](io::Series& series, io::Iteration& iteration)
            {
                (void) series;
                iteration.setAttribute("repetitions", 9u);
            });
        auto const error = expectDynamicError(input);
        REQUIRE(error.find("openPMD validation error for 'dynamic iteration/repetitions'") != std::string::npos);
        REQUIRE(error.find("non-dynamic attribute changed after iteration 0") != std::string::npos);
    }

    SECTION("declared static update")
    {
        auto const input = writeCanonicalStaticDynamicInput(
            "dynamic_contract_static_update",
            [](io::Series& series, io::Iteration& iteration)
            {
                (void) series;
                iteration.setAttribute("haseStaticUpdate", true);
            });
        auto const error = expectDynamicError(input);
        REQUIRE(error.find("openPMD validation error for 'dynamic iteration/haseStaticUpdate'") != std::string::npos);
        REQUIRE(error.find("static updates after iteration 0 are not supported") != std::string::npos);
    }
}

TEST_CASE("openPMD parser round-trips a Python writer contract input", "[openpmd][parser][python]")
{
    char const* inputEnv = std::getenv("HASE_OPENPMD_PYTHON_CONTRACT_INPUT");
    auto const input = inputEnv == nullptr ? writePythonContractEquivalentInput("python_contract_input")
                                           : std::filesystem::path{inputEnv};

    char const* outputEnv = std::getenv("HASE_OPENPMD_PYTHON_CONTRACT_OUTPUT");
    auto const output = outputEnv == nullptr ? testPath("python_contract_result") : std::filesystem::path{outputEnv};
    std::filesystem::remove_all(output);

    hase::openpmd::Parser parser{input, output};
    auto context = parser.read();

    REQUIRE(context.mesh.numberOfTriangles == 3u);
    REQUIRE(context.mesh.numberOfLevels == 6u);
    REQUIRE(context.mesh.numberOfPoints == 5u);
    REQUIRE(context.mesh.thickness == 0.375f);
    REQUIRE(context.mesh.trianglePointIndices == std::vector<unsigned>{0u, 2u, 0u, 1u, 3u, 2u, 2u, 4u, 4u});
    REQUIRE(context.mesh.triangleNeighbors == std::vector<int>{-1, -1, 0, -1, -1, 1, 2, 2, -1});
    REQUIRE(context.mesh.forbiddenEdge == std::vector<int>{-1, -1, 2, -1, -1, 2, 0, 1, -1});
    REQUIRE(context.mesh.triangleNormalPoint == std::vector<unsigned>{0u, 2u, 0u, 1u, 3u, 2u, 2u, 4u, 4u});
    REQUIRE(context.mesh.points == std::vector<double>{0.0, 1.5, 0.25, 2.25, -0.75, 0.0, -0.5, 1.25, 2.5, 3.75});
    requireNear(
        context.mesh.triangleCenterX,
        std::vector<double>{0.58333333333333337, 0.58333333333333337, -0.16666666666666666});
    requireNear(context.mesh.triangleCenterY, std::vector<double>{0.25, 2.5, 1.66666666666666674});
    REQUIRE(
        context.mesh.betaVolume
        == std::vector<
            double>{0.11, 0.21, 0.31, 0.12, 0.22, 0.32, 0.13, 0.23, 0.33, 0.14, 0.24, 0.34, 0.15, 0.25, 0.35});
    REQUIRE(context.mesh.betaCells == std::vector<double>{100.0, 110.0, 120.0, 130.0, 140.0, 101.0, 111.0, 121.0,
                                                          131.0, 141.0, 102.0, 112.0, 122.0, 132.0, 142.0, 103.0,
                                                          113.0, 123.0, 133.0, 143.0, 104.0, 114.0, 124.0, 134.0,
                                                          144.0, 105.0, 115.0, 125.0, 135.0, 145.0});
    REQUIRE(context.mesh.claddingCellTypes == std::vector<unsigned>{0u, 2u, 1u});
    REQUIRE(context.mesh.refractiveIndices == std::vector<float>{1.80f, 1.20f, 1.65f, 1.05f});
    REQUIRE(context.mesh.reflectivities == std::vector<float>{0.01f, 0.03f, 0.05f, 0.02f, 0.04f, 0.06f});
    REQUIRE(context.experiment.spectral == 3u);
    REQUIRE(context.experiment.lambdaA == std::vector<double>{900e-9, 910e-9, 930e-9});
    REQUIRE(context.experiment.lambdaE == std::vector<double>{1000e-9, 1015e-9, 1040e-9});
    REQUIRE(context.experiment.sigmaA == std::vector<double>{0.010, 0.025, 0.040});
    REQUIRE(context.experiment.sigmaE == std::vector<double>{0.050, 0.035, 0.020});
    REQUIRE(context.compute.maxRepetitions == 1u);
    REQUIRE(context.compute.minSampleRange == 0u);
    REQUIRE(context.compute.maxSampleRange == 0u);
    REQUIRE(context.compute.rngSeed == 1234u);

    std::vector<float> phiAse(30u);
    std::vector<double> mse(30u);
    std::vector<unsigned> totalRays(30u);
    std::vector<double> dndtAse(30u);
    for(unsigned i = 0; i < 30u; ++i)
    {
        phiAse[i] = 0.5f + static_cast<float>(i);
        mse[i] = 1000.0 + static_cast<double>(i);
        totalRays[i] = 200u + i;
        dndtAse[i] = -10.0 - static_cast<double>(i);
    }

    parser.writeResult(hase::core::Result{phiAse, mse, totalRays, dndtAse}, context.mesh);
}

#ifdef HASE_OPENPMD_PARSER_VALIDATION_CUSTOM_MAIN
int main(int argc, char* argv[])
{
#    if defined(MPI_FOUND) && !defined(DISABLE_MPI)
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if(!mpiInitialized)
    {
        MPI_Init(&argc, &argv);
    }

    int const result = Catch::Session().run(argc, argv);

    int mpiFinalized = 0;
    MPI_Finalized(&mpiFinalized);
    if(!mpiFinalized)
    {
        MPI_Finalize();
    }
    return result;
#    else
    return Catch::Session().run(argc, argv);
#    endif
}
#endif
