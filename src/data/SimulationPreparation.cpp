#include <data/SimulationPreparation.hpp>
#include <data/TraceData.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace hase::internal::simulationPreparation
{
    using namespace hase::data;
    using hase::data::ExcitationLayout;

    constexpr std::size_t coordinates = 3u;
    constexpr std::size_t verticesPerCell = 4u;
    constexpr std::size_t facesPerCell = 4u;
    constexpr std::size_t verticesPerFace = 3u;

    [[noreturn]] void invalid(std::string const& message)
    {
        throw std::runtime_error("simulation preparation: " + message);
    }

    template<typename T>
    T narrow(std::uint64_t value, std::string const& name)
    {
        if(value > static_cast<std::uint64_t>(std::numeric_limits<T>::max()))
            invalid(name + " is outside the runtime range");
        return static_cast<T>(value);
    }

    template<typename T>
    void requireSize(std::vector<T> const& values, std::size_t expected, std::string const& name)
    {
        if(values.size() != expected)
            invalid(name + " has " + std::to_string(values.size()) + " values; expected " + std::to_string(expected));
    }

    struct Shard
    {
        std::shared_ptr<VolumeTopology> topology;
        std::vector<std::uint8_t> mask;
    };

    std::vector<Shard> shards(Domain const& domain)
    {
        if(domain.masks.offsets.size() != domain.topologies.size() + 1u || domain.masks.offsets.front() != 0u
           || domain.masks.offsets.back() != domain.masks.values.size())
            invalid("domain mask offsets do not match its topology references");
        std::vector<Shard> result;
        result.reserve(domain.topologies.size());
        for(std::size_t index = 0u; index < domain.topologies.size(); ++index)
        {
            auto const begin = domain.masks.offsets[index];
            auto const end = domain.masks.offsets[index + 1u];
            if(begin > end || end > domain.masks.values.size())
                invalid("domain mask offsets are not monotonic");
            result.push_back(
                {domain.topologies[index],
                 {domain.masks.values.begin() + static_cast<std::ptrdiff_t>(begin),
                  domain.masks.values.begin() + static_cast<std::ptrdiff_t>(end)}});
        }
        return result;
    }

    std::size_t cellCount(VolumeTopology const& topology)
    {
        if(topology.cellPointIndices.values.size() % verticesPerCell != 0u)
            invalid("volume topology connectivity is not Tet4");
        auto const count = topology.cellPointIndices.values.size() / verticesPerCell;
        requireSize(topology.cellTypes.values, count, "volume topology cellTypes");
        requireSize(topology.cellDomains.values, count, "volume topology cellDomains");
        requireSize(
            topology.facePointIndices.values,
            count * facesPerCell * verticesPerFace,
            "volume topology facePointIndices");
        requireSize(topology.neighborCells.values, count * facesPerCell, "volume topology neighborCells");
        requireSize(topology.neighborLocalFaces.values, count * facesPerCell, "volume topology neighborLocalFaces");
        requireSize(topology.faceBoundaries.values, count * facesPerCell, "volume topology faceBoundaries");
        requireSize(topology.cellCenters.values, count * coordinates, "volume topology cellCenters");
        requireSize(topology.cellVolumes.values, count, "volume topology cellVolumes");
        return count;
    }

    std::size_t pointCount(VolumeTopology const& topology)
    {
        if(topology.points.values.size() % coordinates != 0u)
            invalid("volume topology points are not three-dimensional");
        return topology.points.values.size() / coordinates;
    }

    template<typename T>
    T const& soa2(std::vector<T> const& values, std::size_t items, std::size_t component, std::size_t item)
    {
        return values.at(component * items + item);
    }

    template<typename T>
    T const& soa3(
        std::vector<T> const& values,
        std::size_t middle,
        std::size_t items,
        std::size_t component,
        std::size_t middleIndex,
        std::size_t item)
    {
        return values.at((component * middle + middleIndex) * items + item);
    }

    struct TopologyMap
    {
        std::shared_ptr<VolumeTopology> topology;
        std::vector<std::uint8_t> selected;
        std::vector<int> cells;
        std::vector<int> points;
    };

    struct Assembly
    {
        std::vector<TopologyMap> topologies;
        std::unordered_map<VolumeTopology const*, std::size_t> lookup;
        std::vector<unsigned> cellPointIndices;
        std::vector<unsigned> cellTypes;
        std::vector<int> cellFaces;
        std::vector<int> neighborCells;
        std::vector<int> neighborLocalFaces;
        std::vector<int> faceBoundaries;
        std::vector<double> cellVolumes;
        std::vector<double> points;
        std::vector<double> samplePoints;
        std::vector<double> cellCenters;
    };

    TopologyMap& topologyMap(Assembly& assembly, std::shared_ptr<VolumeTopology> const& topology)
    {
        if(!topology)
            invalid("domain references a null topology");
        if(auto const found = assembly.lookup.find(topology.get()); found != assembly.lookup.end())
            return assembly.topologies[found->second];
        auto const count = cellCount(*topology);
        auto const index = assembly.topologies.size();
        assembly.lookup.emplace(topology.get(), index);
        assembly.topologies.push_back(
            {topology, std::vector<std::uint8_t>(count, 0u), std::vector<int>(count, -1), {}});
        return assembly.topologies.back();
    }

    void selectDomain(Assembly& assembly, Domain const& domain, bool rejectOverlap)
    {
        if(domain.entityKind != "volume")
            invalid("optical components must reference volume domains");
        for(auto& shard : shards(domain))
        {
            auto& mapping = topologyMap(assembly, shard.topology);
            requireSize(shard.mask, mapping.selected.size(), "volume domain mask");
            for(std::size_t cell = 0u; cell < shard.mask.size(); ++cell)
            {
                if(!shard.mask[cell])
                    continue;
                if(rejectOverlap && mapping.selected[cell])
                    invalid("optical component volume domains overlap");
                mapping.selected[cell] = 1u;
            }
        }
    }

    std::vector<double> structureOfArrays(std::vector<double> const& arrayOfStructures)
    {
        if(arrayOfStructures.size() % coordinates != 0u)
            invalid("coordinate array is not three-dimensional");
        auto const count = arrayOfStructures.size() / coordinates;
        std::vector<double> result(arrayOfStructures.size());
        for(std::size_t item = 0u; item < count; ++item)
            for(std::size_t coordinate = 0u; coordinate < coordinates; ++coordinate)
                result[coordinate * count + item] = arrayOfStructures[item * coordinates + coordinate];
        return result;
    }

    void buildAssemblyArrays(Assembly& assembly)
    {
        std::vector<double> pointsAoS;
        std::vector<double> centersAoS;
        std::size_t nextCell = 0u;
        std::size_t nextPoint = 0u;

        for(auto& mapping : assembly.topologies)
        {
            auto const& topology = *mapping.topology;
            auto const cells = cellCount(topology);
            auto const points = pointCount(topology);
            mapping.points.assign(points, -1);

            std::vector<std::uint8_t> used(points, 0u);
            bool const complete = std::all_of(
                mapping.selected.begin(),
                mapping.selected.end(),
                [](auto value) { return value != 0u; });
            if(complete)
                std::fill(used.begin(), used.end(), 1u);
            else
            {
                for(std::size_t cell = 0u; cell < cells; ++cell)
                    if(mapping.selected[cell])
                        for(std::size_t vertex = 0u; vertex < verticesPerCell; ++vertex)
                        {
                            auto const point = soa2(topology.cellPointIndices.values, cells, vertex, cell);
                            if(point >= points)
                                invalid("volume topology connectivity references an invalid point");
                            used[point] = 1u;
                        }
            }
            for(std::size_t point = 0u; point < points; ++point)
            {
                if(!used[point])
                    continue;
                mapping.points[point] = static_cast<int>(nextPoint++);
                pointsAoS.insert(
                    pointsAoS.end(),
                    {soa2(topology.points.values, points, 0u, point),
                     soa2(topology.points.values, points, 1u, point),
                     soa2(topology.points.values, points, 2u, point)});
            }
            for(std::size_t cell = 0u; cell < cells; ++cell)
                if(mapping.selected[cell])
                    mapping.cells[cell] = static_cast<int>(nextCell++);
        }

        for(auto const& mapping : assembly.topologies)
        {
            auto const& topology = *mapping.topology;
            auto const cells = mapping.selected.size();
            for(std::size_t cell = 0u; cell < cells; ++cell)
            {
                if(!mapping.selected[cell])
                    continue;
                for(std::size_t vertex = 0u; vertex < verticesPerCell; ++vertex)
                {
                    auto const source = soa2(topology.cellPointIndices.values, cells, vertex, cell);
                    auto const target = mapping.points.at(source);
                    if(target < 0)
                        invalid("selected cell references an unselected point");
                    assembly.cellPointIndices.push_back(static_cast<unsigned>(target));
                }
                assembly.cellTypes.push_back(topology.cellTypes.values[cell]);
                for(std::size_t face = 0u; face < facesPerCell; ++face)
                {
                    for(std::size_t vertex = 0u; vertex < verticesPerFace; ++vertex)
                    {
                        auto const source
                            = soa3(topology.facePointIndices.values, facesPerCell, cells, vertex, face, cell);
                        if(source < 0 || static_cast<std::size_t>(source) >= mapping.points.size()
                           || mapping.points[source] < 0)
                            invalid("selected cell face references an invalid point");
                        assembly.cellFaces.push_back(mapping.points[source]);
                    }
                    auto const sourceNeighbor = soa2(topology.neighborCells.values, cells, face, cell);
                    int targetNeighbor = -1;
                    if(sourceNeighbor >= 0 && static_cast<std::size_t>(sourceNeighbor) < mapping.cells.size())
                        targetNeighbor = mapping.cells[sourceNeighbor];
                    assembly.neighborCells.push_back(targetNeighbor);
                    assembly.neighborLocalFaces.push_back(
                        targetNeighbor < 0 ? -1 : soa2(topology.neighborLocalFaces.values, cells, face, cell));
                    auto boundary = soa2(topology.faceBoundaries.values, cells, face, cell);
                    if(targetNeighbor < 0 && boundary == 0)
                        boundary = -1;
                    assembly.faceBoundaries.push_back(boundary);
                }
                assembly.cellVolumes.push_back(topology.cellVolumes.values[cell]);
                centersAoS.insert(
                    centersAoS.end(),
                    {soa2(topology.cellCenters.values, cells, 0u, cell),
                     soa2(topology.cellCenters.values, cells, 1u, cell),
                     soa2(topology.cellCenters.values, cells, 2u, cell)});
            }
        }
        if(assembly.cellTypes.empty())
            invalid("optical assembly is empty");
        assembly.points = structureOfArrays(pointsAoS);
        assembly.cellCenters = structureOfArrays(centersAoS);
        assembly.samplePoints = assembly.cellCenters;
    }

    std::vector<unsigned> domainCells(Assembly const& assembly, Domain const& domain)
    {
        if(domain.entityKind != "volume")
            invalid("cell assignment references a non-volume domain");
        std::vector<unsigned> result;
        for(auto const& shard : shards(domain))
        {
            auto const found = assembly.lookup.find(shard.topology.get());
            if(found == assembly.lookup.end())
                invalid("domain contains cells outside the optical assembly");
            auto const& mapping = assembly.topologies[found->second];
            requireSize(shard.mask, mapping.cells.size(), "volume domain mask");
            for(std::size_t cell = 0u; cell < shard.mask.size(); ++cell)
            {
                if(!shard.mask[cell])
                    continue;
                if(mapping.cells[cell] < 0)
                    invalid("domain contains cells outside the optical assembly");
                result.push_back(static_cast<unsigned>(mapping.cells[cell]));
            }
        }
        return result;
    }

    bool sameDomain(Domain const& left, Domain const& right)
    {
        if(left.entityKind != right.entityKind || left.topologies.size() != right.topologies.size()
           || left.masks.offsets != right.masks.offsets || left.masks.values != right.masks.values)
            return false;
        for(std::size_t index = 0u; index < left.topologies.size(); ++index)
            if(left.topologies[index].get() != right.topologies[index].get())
                return false;
        return true;
    }

    class SurfaceDomains
    {
    public:
        explicit SurfaceDomains(Assembly& assembly) : m_assembly(assembly)
        {
            for(auto boundary : assembly.faceBoundaries)
                if(boundary > 0)
                    m_next = std::max(m_next, boundary + 1);
        }

        int assign(std::shared_ptr<Domain> const& domain)
        {
            if(!domain || domain->entityKind != "surface")
                invalid("surface assignment references a non-surface domain");
            for(auto const& [known, identifier] : m_domains)
                if(sameDomain(*known, *domain))
                    return identifier;
            int const identifier = m_next++;
            bool any = false;
            for(auto const& shard : shards(*domain))
            {
                auto const found = m_assembly.lookup.find(shard.topology.get());
                if(found == m_assembly.lookup.end())
                    invalid("surface domain is outside the optical assembly");
                auto const& mapping = m_assembly.topologies[found->second];
                requireSize(shard.mask, mapping.cells.size() * facesPerCell, "surface domain mask");
                for(std::size_t index = 0u; index < shard.mask.size(); ++index)
                {
                    if(!shard.mask[index])
                        continue;
                    auto const cell = index / facesPerCell;
                    auto const face = index % facesPerCell;
                    if(mapping.cells[cell] < 0)
                        invalid("surface domain is outside the selected optical assembly");
                    m_assembly.faceBoundaries[static_cast<std::size_t>(mapping.cells[cell]) * facesPerCell + face]
                        = identifier;
                    any = true;
                }
            }
            if(!any)
                invalid("surface domain is empty");
            m_domains.emplace_back(domain, identifier);
            return identifier;
        }

    private:
        Assembly& m_assembly;
        int m_next = 1;
        std::vector<std::pair<std::shared_ptr<Domain>, int>> m_domains;
    };

    bool domainsOverlap(Domain const& left, Domain const& right)
    {
        for(auto const& leftShard : shards(left))
            for(auto const& rightShard : shards(right))
                if(leftShard.topology.get() == rightShard.topology.get())
                {
                    requireSize(rightShard.mask, leftShard.mask.size(), "surface domain mask");
                    for(std::size_t index = 0u; index < leftShard.mask.size(); ++index)
                        if(leftShard.mask[index] && rightShard.mask[index])
                            return true;
                }
        return false;
    }

    ExcitationLayout excitationUpdatePlan(Simulation const& simulation, Assembly const& assembly)
    {
        if(!simulation.excitationState || !simulation.gainMedium)
            invalid("simulation has no excitation state or gain medium");
        std::vector<std::uint8_t> active(assembly.cellTypes.size(), 0u);
        for(auto const& component : simulation.gainMedium->components)
            for(auto cell : domainCells(assembly, *component->domain))
                active[cell] = 1u;
        std::vector<std::vector<unsigned>> domainCellIndices;
        domainCellIndices.reserve(simulation.excitationState->domains.size());
        for(auto const& domain : simulation.excitationState->domains)
            domainCellIndices.push_back(domainCells(assembly, *domain));
        return ExcitationLayout{std::move(active), std::move(domainCellIndices)};
    }

    void validateCrossSections(CrossSectionTable const& table, std::string const& materialName)
    {
        auto const& wavelengths = table.wavelengths.values;
        if(wavelengths.empty() || wavelengths.size() != table.absorption.values.size()
           || wavelengths.size() != table.emission.values.size())
            invalid("material '" + materialName + "' has an empty or inconsistent cross-section table");
        for(std::size_t sample = 0u; sample < wavelengths.size(); ++sample)
        {
            if(!std::isfinite(wavelengths[sample]) || wavelengths[sample] <= 0.0
               || !std::isfinite(table.absorption.values[sample]) || !std::isfinite(table.emission.values[sample])
               || table.absorption.values[sample] < 0.0 || table.emission.values[sample] < 0.0)
                invalid("material '" + materialName + "' has invalid cross-section samples");
        }
        bool const monochromatic
            = std::ranges::all_of(wavelengths, [&](double wavelength) { return wavelength == wavelengths.front(); });
        if(!monochromatic)
            for(std::size_t sample = 1u; sample < wavelengths.size(); ++sample)
                if(wavelengths[sample] <= wavelengths[sample - 1u])
                    invalid("material '" + materialName + "' wavelengths must be strictly increasing");
    }

    struct MaterialAssembly
    {
        std::vector<unsigned> cellMaterialIds;
        std::vector<std::uint8_t> active;
        std::vector<double> refractiveIndices;
        std::vector<double> activeIonDensities;
        std::vector<double> fluorescenceLifetimes;
        std::vector<double> bulkAttenuations;
        std::vector<double> peakAbsorption;
        std::vector<double> peakEmission;
        std::vector<unsigned> crossSectionOffsets{0u};
        std::vector<double> wavelengths;
        std::vector<double> absorption;
        std::vector<double> emission;
    };

    MaterialAssembly assembleMaterials(Simulation const& simulation, Assembly const& assembly)
    {
        MaterialAssembly result;
        result.cellMaterialIds.assign(assembly.cellTypes.size(), std::numeric_limits<unsigned>::max());
        std::unordered_map<Material const*, unsigned> materialIds;

        for(auto const& component : simulation.opticalComponents)
        {
            auto const& material = component->material;
            auto [position, inserted] = materialIds.emplace(material.get(), static_cast<unsigned>(materialIds.size()));
            unsigned const materialId = position->second;
            if(inserted)
            {
                if(!std::isfinite(material->refractiveIndex) || material->refractiveIndex <= 0.0)
                    invalid("material '" + material->materialName + "' has an invalid refractive index");
                if(!std::isfinite(material->activeIonDensity) || material->activeIonDensity < 0.0)
                    invalid("material '" + material->materialName + "' has an invalid active-ion density");
                double const lifetime = material->fluorescenceLifetime.value_or(0.0);
                if(material->active
                   && (material->activeIonDensity <= 0.0 || lifetime <= 0.0 || !material->crossSections))
                    invalid(
                        "active material '" + material->materialName
                        + "' requires positive active-ion density, fluorescence lifetime, and cross sections");
                double const attenuation = material->bulkAttenuation.value_or(0.0);
                if(!std::isfinite(attenuation) || attenuation < 0.0)
                    invalid("material '" + material->materialName + "' has invalid bulk attenuation");

                result.active.push_back(material->active ? 1u : 0u);
                result.refractiveIndices.push_back(material->refractiveIndex);
                result.activeIonDensities.push_back(material->activeIonDensity);
                result.fluorescenceLifetimes.push_back(lifetime);
                result.bulkAttenuations.push_back(attenuation);

                double peakAbsorption = 0.0;
                double peakEmission = 0.0;
                if(material->crossSections)
                {
                    validateCrossSections(*material->crossSections, material->materialName);
                    auto const& table = *material->crossSections;
                    auto const peak = std::max_element(table.emission.values.begin(), table.emission.values.end());
                    auto const peakIndex = static_cast<std::size_t>(peak - table.emission.values.begin());
                    peakAbsorption = table.absorption.values[peakIndex];
                    peakEmission = table.emission.values[peakIndex];
                    result.wavelengths.insert(
                        result.wavelengths.end(),
                        table.wavelengths.values.begin(),
                        table.wavelengths.values.end());
                    result.absorption.insert(
                        result.absorption.end(),
                        table.absorption.values.begin(),
                        table.absorption.values.end());
                    result.emission.insert(
                        result.emission.end(),
                        table.emission.values.begin(),
                        table.emission.values.end());
                }
                result.peakAbsorption.push_back(peakAbsorption);
                result.peakEmission.push_back(peakEmission);
                result.crossSectionOffsets.push_back(static_cast<unsigned>(result.wavelengths.size()));
            }

            for(auto const cell : domainCells(assembly, *component->domain))
            {
                if(result.cellMaterialIds[cell] != std::numeric_limits<unsigned>::max())
                    invalid("optical component domains overlap");
                result.cellMaterialIds[cell] = materialId;
            }
        }

        if(std::ranges::find(result.cellMaterialIds, std::numeric_limits<unsigned>::max())
           != result.cellMaterialIds.end())
            invalid("optical components do not assign a material to every trace cell");
        return result;
    }

    core::PumpProfileParameters pumpProfile(PumpProfile const& profile)
    {
        core::PumpProfileParameters result;
        std::visit(
            [&result](auto const& value)
            {
                using T = std::decay_t<decltype(value)>;
                if constexpr(std::is_same_v<T, UniformPumpProfile>)
                    result.kind = 0u;
                else
                {
                    result.kind = 1u;
                    result.radiusU = value.radiusU;
                    result.radiusV = value.radiusV;
                    result.exponent = value.exponent;
                    requireSize(value.center.values, coordinates, "pump profile center");
                    requireSize(value.axisU.values, coordinates, "pump profile axisU");
                    requireSize(value.axisV.values, coordinates, "pump profile axisV");
                    for(std::size_t coordinate = 0u; coordinate < coordinates; ++coordinate)
                    {
                        result.center[coordinate] = value.center.values[coordinate];
                        result.axisU[coordinate] = value.axisU.values[coordinate];
                        result.axisV[coordinate] = value.axisV.values[coordinate];
                    }
                }
            },
            profile.value);
        return result;
    }
} // namespace hase::internal::simulationPreparation

namespace hase::data
{
    using namespace hase::internal::simulationPreparation;

    ExcitationLayout::ExcitationLayout(
        std::vector<std::uint8_t> active,
        std::vector<std::vector<unsigned>> domainCells)
        : m_active(std::move(active))
        , m_domainCells(std::move(domainCells))
    {
    }

    std::vector<double> ExcitationLayout::values(ExcitationState const& state) const
    {
        if(state.values.offsets.size() != m_domainCells.size() + 1u || state.values.offsets.empty()
           || state.values.offsets.front() != 0u || state.values.offsets.back() != state.values.values.size())
            invalid("excitation value offsets do not match its domains");
        std::vector<double> result(m_active.size(), 0.0);
        std::vector<std::uint8_t> covered(m_active.size(), 0u);
        for(std::size_t domainIndex = 0u; domainIndex < m_domainCells.size(); ++domainIndex)
        {
            auto const& cells = m_domainCells[domainIndex];
            auto const begin = state.values.offsets[domainIndex];
            auto const end = state.values.offsets[domainIndex + 1u];
            auto const count = end - begin;
            if(count != 1u && count != cells.size())
                invalid("excitation values must be scalar or match their domain size");
            for(std::size_t index = 0u; index < cells.size(); ++index)
            {
                auto const cell = cells[index];
                if(!m_active[cell])
                    invalid("excitation domain extends outside the gain medium");
                if(covered[cell])
                    invalid("excitation domains overlap");
                auto const value = state.values.values[begin + (count == 1u ? 0u : index)];
                if(!std::isfinite(value) || value < 0.0 || value > 1.0)
                    invalid("excitation values must be finite and within [0, 1]");
                result[cell] = value;
                covered[cell] = 1u;
            }
        }
        for(std::size_t cell = 0u; cell < m_active.size(); ++cell)
            if(m_active[cell] && !covered[cell])
                invalid("excitation domains do not cover the gain medium");
        return result;
    }

    void ExcitationLayout::apply(ExcitationState const& state, std::vector<double>& destination) const
    {
        destination = values(state);
    }

    SimulationState prepareSimulation(Simulation const& simulation)
    {
        return std::move(prepareSimulationWithUpdates(simulation).state);
    }

    SimulationPreparation prepareSimulationWithUpdates(Simulation const& simulation)
    {
        if(simulation.opticalComponents.empty() || !simulation.phiAse || !simulation.timeIntegrationSolver)
            invalid("simulation is missing required primitive references");

        Assembly assembly;
        for(auto const& component : simulation.opticalComponents)
        {
            if(!component || !component->domain || !component->material)
                invalid("optical component is incomplete");
            selectDomain(assembly, *component->domain, true);
        }
        buildAssemblyArrays(assembly);
        auto excitationPlan = excitationUpdatePlan(simulation, assembly);
        auto betaVolume = excitationPlan.values(*simulation.excitationState);
        auto materials = assembleMaterials(simulation, assembly);

        SurfaceDomains surfaceDomains(assembly);

        struct OpticsEntry
        {
            int identifier;
            std::shared_ptr<Domain> domain;
            std::shared_ptr<SurfaceOptics> optics;
        };

        std::vector<OpticsEntry> optics;
        for(auto const& component : simulation.opticalComponents)
            for(auto const& assignment : component->surfaceOptics)
            {
                if(!assignment || !assignment->domain || !assignment->optics)
                    invalid("surface-optics assignment is incomplete");
                auto const identifier = surfaceDomains.assign(assignment->domain);
                auto const existing = std::find_if(
                    optics.begin(),
                    optics.end(),
                    [identifier](auto const& value) { return value.identifier == identifier; });
                auto const overlap = std::find_if(
                    optics.begin(),
                    optics.end(),
                    [&](auto const& value) { return domainsOverlap(*value.domain, *assignment->domain); });
                if(overlap != optics.end() && overlap->identifier != identifier)
                    invalid("surface-optics assignments overlap on one face");
                if(existing != optics.end()
                   && (existing->optics->reflectivity != assignment->optics->reflectivity
                       || existing->optics->nInside != assignment->optics->nInside
                       || existing->optics->nOutside != assignment->optics->nOutside))
                    invalid("one surface domain has conflicting optical properties");
                if(existing == optics.end())
                    optics.push_back({identifier, assignment->domain, assignment->optics});
            }

        core::SimulationControls run;
        run.timeStep = simulation.timeStep;
        if(!simulation.simulationSteps)
            invalid("simulation has no requested step count");
        run.numberOfSteps = narrow<unsigned>(*simulation.simulationSteps, "simulationSteps");
        run.firstSimulationStep = narrow<unsigned>(simulation.currentStep, "currentStep");
        run.aseSteps = narrow<unsigned>(simulation.phiAse->aseSteps.value_or(0u), "aseSteps");
        run.prePump = simulation.prePump;
        run.executionMode = simulation.executionMode;
        if(simulation.outputSteps)
            for(auto step : simulation.outputSteps->values)
                run.outputSteps.push_back(narrow<unsigned>(step, "outputSteps"));
        run.outputFields = simulation.outputFields;
        run.controlFields = simulation.controlFields;
        run.timeIntegration.method = simulation.timeIntegrationSolver->name;
        run.timeIntegration.implicitIterations
            = narrow<unsigned>(simulation.timeIntegrationSolver->iterations.value_or(8u), "timeIntegrator.iterations");
        run.timeIntegration.implicitTolerance = simulation.timeIntegrationSolver->tolerance.value_or(1.0e-10);

        for(auto const& registration : simulation.pumpRegistrations)
        {
            if(!registration || !registration->pump || !registration->injectionMethod)
                invalid("pump registration is incomplete");
            auto const& pump = *registration->pump;
            if(!pump.spectrum || !pump.profile || !pump.angularDistribution)
                invalid("pump is missing its spectrum, profile, or angular distribution");
            core::PumpSourceParameters source;
            source.totalPower = pump.totalPower;
            source.rayCount = narrow<unsigned>(pump.rayCount, "pump.rayCount");
            source.pumpSteps = narrow<unsigned>(pump.pumpSteps.value_or(0u), "pump.pumpSteps");
            source.rngSeed = narrow<std::uint32_t>(pump.rngSeed, "pump.rngSeed");
            for(auto const& domain : registration->injectionMethod->surfaceDomains)
                source.surfaces.push_back(surfaceDomains.assign(domain));
            source.wavelengths = pump.spectrum->wavelengths.values;
            source.spectralWeights = pump.spectrum->weights.values;
            source.polarAngles = pump.angularDistribution->polarAngles.values;
            source.azimuthalAngles = pump.angularDistribution->azimuthalAngles.values;
            source.angularWeights = pump.angularDistribution->weights.values;
            source.profile = pumpProfile(*pump.profile);
            for(auto const& relayValue : registration->relays)
            {
                core::PumpRelayParameters relay;
                for(auto const& domain : relayValue->exitDomains)
                    relay.exitSurfaces.push_back(surfaceDomains.assign(domain));
                for(auto const& domain : relayValue->entryDomains)
                    relay.entrySurfaces.push_back(surfaceDomains.assign(domain));
                relay.flipU = relayValue->flipU;
                relay.flipV = relayValue->flipV;
                relay.rotation = relayValue->rotation;
                requireSize(relayValue->offset.values, 2u, "pump relay offset");
                requireSize(relayValue->tilt.values, 2u, "pump relay tilt");
                relay.offset[0] = relayValue->offset.values[0];
                relay.offset[1] = relayValue->offset.values[1];
                relay.tilt[0] = relayValue->tilt.values[0];
                relay.tilt[1] = relayValue->tilt.values[1];
                relay.magnification = relayValue->magnification;
                relay.transmission = relayValue->transmission;
                source.relays.push_back(std::move(relay));
            }
            run.pump.sources.push_back(std::move(source));
        }

        int maximumSurface = 0;
        for(auto boundary : assembly.faceBoundaries)
            maximumSurface = std::max(maximumSurface, boundary);
        auto const surfaceCount = maximumSurface == 0 ? 0u : static_cast<std::size_t>(maximumSurface + 1);
        std::vector<float> surfaceReflectivity(surfaceCount, 0.0f);
        std::vector<float> surfaceInside(surfaceCount, 1.0f);
        std::vector<float> surfaceOutside(surfaceCount, 1.0f);
        for(auto const& entry : optics)
        {
            surfaceReflectivity.at(entry.identifier) = static_cast<float>(entry.optics->reflectivity);
            surfaceInside.at(entry.identifier) = static_cast<float>(entry.optics->nInside);
            surfaceOutside.at(entry.identifier) = static_cast<float>(entry.optics->nOutside);
        }

        TraceData trace(
            std::move(assembly.cellPointIndices),
            std::move(assembly.cellTypes),
            std::move(assembly.cellFaces),
            std::move(assembly.neighborCells),
            std::move(assembly.neighborLocalFaces),
            std::move(assembly.faceBoundaries),
            std::move(assembly.cellVolumes),
            std::move(assembly.points),
            std::move(assembly.samplePoints),
            std::move(assembly.cellCenters),
            std::move(betaVolume),
            std::move(materials.cellMaterialIds),
            std::move(materials.active),
            std::move(materials.refractiveIndices),
            std::move(materials.activeIonDensities),
            std::move(materials.fluorescenceLifetimes),
            std::move(materials.bulkAttenuations),
            std::move(materials.peakAbsorption),
            std::move(materials.peakEmission),
            std::move(materials.crossSectionOffsets),
            std::move(materials.wavelengths),
            std::move(materials.absorption),
            std::move(materials.emission),
            std::move(surfaceReflectivity),
            std::move(surfaceInside),
            std::move(surfaceOutside));

        core::AseTraceControls ase;
        ase.minRays = narrow<unsigned>(simulation.phiAse->minRays, "phiAse.minRays");
        ase.maxRays = narrow<unsigned>(simulation.phiAse->maxRays, "phiAse.maxRays");
        ase.relativeStandardErrorThreshold = simulation.phiAse->relativeStandardErrorThreshold;
        ase.trackRayVisits = simulation.phiAse->trackRayVisits;
        ase.useReflections = simulation.phiAse->useReflections;
        ase.reflectionMode = simulation.phiAse->reflectionMode;
        ase.surfaceReservoirSize
            = narrow<unsigned>(simulation.phiAse->surfaceReservoirSize, "phiAse.surfaceReservoirSize");
        ase.srmPositionMode = simulation.phiAse->srmPositionMode;
        if(ase.reflectionMode != "direct" && ase.reflectionMode != "srm")
            invalid("phiAse.reflectionMode must be 'direct' or 'srm'");
        if(ase.surfaceReservoirSize == 0u)
            invalid("phiAse.surfaceReservoirSize must be positive");
        if(ase.srmPositionMode != "exact" && ase.srmPositionMode != "centroid")
            invalid("phiAse.srmPositionMode must be 'exact' or 'centroid'");
        ase.monochromatic = simulation.phiAse->monochromatic;
        ase.propagationMode = simulation.phiAse->propagationMode;
        ase.forwardRayCount
            = narrow<unsigned>(simulation.phiAse->forwardRayCount.value_or(0u), "phiAse.forwardRayCount");
        ase.reflectionMaxIterations
            = narrow<unsigned>(simulation.phiAse->reflectionMaxIterations, "phiAse.reflectionMaxIterations");
        ase.reflectionTolerance = simulation.phiAse->reflectionTolerance;

        if(!simulation.phiAse->backend)
            invalid("phiAse backend was not resolved before transport");
        core::ExecutionPolicy compute(
            narrow<unsigned>(simulation.phiAse->repetitions, "phiAse.repetitions"),
            narrow<unsigned>(simulation.phiAse->adaptiveSteps, "phiAse.adaptiveSteps"),
            narrow<unsigned>(simulation.phiAse->numDevices, "phiAse.numDevices"),
            0u,
            *simulation.phiAse->backend,
            simulation.phiAse->parallelMode,
            simulation.phiAse->writeVtk,
            [&simulation]
            {
                std::vector<unsigned> devices;
                if(!simulation.phiAse->devices)
                    return devices;
                devices.reserve(simulation.phiAse->devices->values.size());
                for(auto const device : simulation.phiAse->devices->values)
                    devices.push_back(narrow<unsigned>(device, "phiAse.devices"));
                return devices;
            }(),
            narrow<unsigned>(simulation.phiAse->minSampleRange.value_or(0u), "phiAse.minSampleRange"),
            narrow<unsigned>(
                simulation.phiAse->maxSampleRange.value_or(trace.numberOfCells - 1u),
                "phiAse.maxSampleRange"),
            narrow<unsigned>(
                simulation.phiAse->rngSeed.value_or(core::ExecutionPolicy::unspecifiedRngSeed),
                "phiAse.rngSeed"));

        PhiAseResult result(
            std::vector<float>(trace.numberOfCells, 0.0f),
            std::vector<double>(trace.numberOfCells, 0.0),
            std::vector<double>(trace.numberOfCells, 0.0),
            std::vector<unsigned>(trace.numberOfCells, 0u),
            std::vector<double>(trace.numberOfCells, 0.0));
        SimulationState state{std::move(ase), std::move(compute), std::move(trace), std::move(result), std::move(run)};
        return {std::move(state), std::move(excitationPlan)};
    }
} // namespace hase::data
