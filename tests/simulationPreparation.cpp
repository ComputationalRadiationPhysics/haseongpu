#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <data/SimulationPreparation.hpp>

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

TEST_CASE("cross-section replacement is explicit and may resize", "[transport][material][update]")
{
    hase::data::CrossSectionTable table;
    table.replaceSamples({900e-9}, {1.0e-25}, {2.0e-25});
    table.replaceSamples({900e-9, 950e-9, 1000e-9}, {1.0e-25, 2.0e-25, 3.0e-25}, {4.0e-25, 5.0e-25, 6.0e-25});

    CHECK(table.wavelengths.values.size() == 3u);
    CHECK(table.absorption.values.back() == Catch::Approx(3.0e-25));
    CHECK(table.emission.values.back() == Catch::Approx(6.0e-25));
    CHECK_THROWS_AS(table.replaceSamples({900e-9, 850e-9}, {1.0, 1.0}, {1.0, 1.0}), std::invalid_argument);
    CHECK(table.wavelengths.values == std::vector<double>{900e-9, 950e-9, 1000e-9});
    CHECK(table.absorption.values == std::vector<double>{1.0e-25, 2.0e-25, 3.0e-25});
    CHECK(table.emission.values == std::vector<double>{4.0e-25, 5.0e-25, 6.0e-25});
}

TEST_CASE("simulation preparation builds material-local trace data", "[transport][preparation]")
{
    using namespace hase::data;

    auto topology = std::make_shared<VolumeTopology>();
    topology->points.values = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    topology->cellPointIndices.values = {0u, 1u, 2u, 3u};
    topology->cellTypes.values = {10u};
    topology->cellDomains.values = {1};
    topology->facePointIndices.values = {0, 0, 1, 2, 2, 1, 2, 0, 1, 3, 3, 3};
    topology->neighborCells.values = {-1, -1, -1, -1};
    topology->neighborLocalFaces.values = {-1, -1, -1, -1};
    topology->faceBoundaries.values = {-1, -1, -1, -1};
    topology->cellCenters.values = {0.25, 0.25, 0.25};
    topology->cellVolumes.values = {1.0 / 6.0};

    auto volume = std::make_shared<Domain>();
    volume->entityKind = "volume";
    volume->masks = {{1u}, {0u, 1u}};
    volume->topologies = {topology};

    auto table = std::make_shared<CrossSectionTable>();
    table->wavelengths.values = {940.0e-9};
    table->absorption.values = {1.2e-25};
    table->emission.values = {2.1e-24};

    auto material = std::make_shared<Material>();
    material->materialName = "test gain material";
    material->active = true;
    material->refractiveIndex = 1.8;
    material->fluorescenceLifetime = 9.5e-4;
    material->activeIonDensity = 2.76e26;
    material->crossSections = table;

    auto component = std::make_shared<OpticalComponent>();
    component->domain = volume;
    component->material = material;

    auto gainMedium = std::make_shared<GainMedium>();
    gainMedium->components = {component};

    auto excitation = std::make_shared<ExcitationState>();
    excitation->domains = {volume};
    excitation->values = {{0.25}, {0u, 1u}};

    auto phiAse = std::make_shared<PhiAse>();
    phiAse->propagationMode = "forward";
    phiAse->minRays = 10u;
    phiAse->maxRays = 20u;
    phiAse->relativeStandardErrorThreshold = 0.1;
    phiAse->trackRayVisits = true;
    phiAse->repetitions = 2u;
    phiAse->adaptiveSteps = 3u;
    phiAse->reflectionMode = "srm";
    phiAse->surfaceReservoirSize = 256u;
    phiAse->srmPositionMode = "centroid";
    phiAse->reflectionMaxIterations = 40u;
    phiAse->reflectionTolerance = 1.0e-4;
    phiAse->backend = "test-backend";
    phiAse->parallelMode = "single";
    phiAse->numDevices = 1u;
    phiAse->devices = hase::transport::Array<std::uint64_t>{};

    auto solver = std::make_shared<TimeIntegrationSolver>();
    solver->name = "explicit-euler";

    Simulation simulation;
    simulation.timeStep = 1.0e-6;
    simulation.simulationSteps = 2u;
    simulation.executionMode = "autonomous";
    simulation.outputFields = {"beta_volume", "phi_ase"};
    simulation.opticalComponents = {component};
    simulation.gainMedium = gainMedium;
    simulation.exteriorSurface = std::make_shared<Domain>();
    simulation.excitationState = excitation;
    simulation.phiAse = phiAse;
    simulation.timeIntegrationSolver = solver;

    auto preparation = prepareSimulationWithUpdates(simulation);
    auto& state = preparation.state;

    CHECK(state.trace.numberOfCells == 1u);
    CHECK(state.trace.numberOfMeshPoints == 4u);
    CHECK(state.trace.betaVolume == std::vector<double>{0.25});
    CHECK(state.trace.numberOfMaterials == 1u);
    CHECK(state.trace.cellMaterialIds == std::vector<unsigned>{0u});
    CHECK(state.trace.materialActiveIonDensities.front() == Catch::Approx(2.76e26));
    CHECK(state.trace.materialFluorescenceLifetimes.front() == Catch::Approx(9.5e-4));
    CHECK(state.trace.crossSectionAbsorption.front() == Catch::Approx(1.2e-25));
    CHECK(state.trace.crossSectionEmission.front() == Catch::Approx(2.1e-24));
    CHECK(state.controls.numberOfSteps == 2u);
    CHECK(state.controls.timeIntegration.method == "explicit-euler");
    CHECK(state.ase.reflectionMode == "srm");
    CHECK(state.ase.surfaceReservoirSize == 256u);
    CHECK(state.ase.srmPositionMode == "centroid");
    CHECK(state.ase.trackRayVisits);

    excitation->values.values = {0.75};
    preparation.excitation.apply(*excitation, state.trace.betaVolume);
    CHECK(state.trace.betaVolume == std::vector<double>{0.75});

    table->replaceSamples({940.0e-9, 1030.0e-9}, {1.3e-25, 1.4e-25}, {2.2e-24, 2.3e-24});
    auto updated = prepareSimulation(simulation);
    CHECK_FALSE(state.trace.hasSameMaterialData(updated.trace));
    state.trace.replaceMaterialData(updated.trace);
    CHECK(state.trace.cellMaterialIds == std::vector<unsigned>{0u});
    CHECK(state.trace.materialCrossSectionOffsets == std::vector<unsigned>{0u, 2u});
    CHECK(state.trace.crossSectionWavelengths == std::vector<double>{940.0e-9, 1030.0e-9});
    CHECK(state.trace.crossSectionAbsorption == std::vector<double>{1.3e-25, 1.4e-25});
    CHECK(state.trace.crossSectionEmission == std::vector<double>{2.2e-24, 2.3e-24});
    CHECK(state.trace.hasSameMaterialData(updated.trace));
}

TEST_CASE(
    "simulation preparation assembles independent domain-local topology shards",
    "[transport][preparation][partition]")
{
    using namespace hase::data;

    auto makeTopology = [](double offset)
    {
        auto topology = std::make_shared<VolumeTopology>();
        topology->points.values = {offset, offset + 1.0, offset, offset, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
        topology->cellPointIndices.values = {0u, 1u, 2u, 3u};
        topology->cellTypes.values = {10u};
        topology->cellDomains.values = {1};
        topology->facePointIndices.values = {0, 0, 1, 2, 2, 1, 2, 0, 1, 3, 3, 3};
        topology->neighborCells.values = {-1, -1, -1, -1};
        topology->neighborLocalFaces.values = {-1, -1, -1, -1};
        topology->faceBoundaries.values = {-1, -1, -1, -1};
        topology->cellCenters.values = {offset + 0.25, 0.25, 0.25};
        topology->cellVolumes.values = {1.0 / 6.0};
        return topology;
    };
    auto makeDomain = [](std::shared_ptr<VolumeTopology> const& topology)
    {
        auto domain = std::make_shared<Domain>();
        domain->entityKind = "volume";
        domain->masks = {{1u}, {0u, 1u}};
        domain->topologies = {topology};
        return domain;
    };

    auto makeMaterial =
        [](std::string name, double wavelength, double absorption, double emission, double density, double attenuation)
    {
        auto table = std::make_shared<CrossSectionTable>();
        table->wavelengths.values = {wavelength};
        table->absorption.values = {absorption};
        table->emission.values = {emission};
        auto material = std::make_shared<Material>();
        material->materialName = std::move(name);
        material->active = true;
        material->refractiveIndex = 1.8;
        material->fluorescenceLifetime = 1.0e-3;
        material->activeIonDensity = density;
        material->bulkAttenuation = attenuation;
        material->crossSections = std::move(table);
        return material;
    };
    std::array materials{
        makeMaterial("first gain material", 1030.0e-9, 1.0e-25, 2.0e-24, 2.76e26, 0.5),
        makeMaterial("second gain material", 1064.0e-9, 3.0e-25, 4.0e-24, 3.10e26, 1.5)};

    std::vector<std::shared_ptr<OpticalComponent>> components;
    std::vector<std::shared_ptr<Domain>> domains;
    for(std::size_t index = 0u; index < materials.size(); ++index)
    {
        auto domain = makeDomain(makeTopology(index == 0u ? 0.0 : 10.0));
        auto component = std::make_shared<OpticalComponent>();
        component->domain = domain;
        component->material = materials[index];
        domains.push_back(domain);
        components.push_back(component);
    }

    auto gainMedium = std::make_shared<GainMedium>();
    gainMedium->components = components;
    auto excitation = std::make_shared<ExcitationState>();
    excitation->domains = domains;
    excitation->values = {{0.2, 0.4}, {0u, 1u, 2u}};
    auto phiAse = std::make_shared<PhiAse>();
    phiAse->propagationMode = "forward";
    phiAse->minRays = 10u;
    phiAse->maxRays = 10u;
    phiAse->relativeStandardErrorThreshold = 0.1;
    phiAse->repetitions = 1u;
    phiAse->adaptiveSteps = 0u;
    phiAse->reflectionMode = "direct";
    phiAse->surfaceReservoirSize = 64u;
    phiAse->srmPositionMode = "exact";
    phiAse->reflectionMaxIterations = 40u;
    phiAse->reflectionTolerance = 1.0e-4;
    phiAse->backend = "test-backend";
    phiAse->parallelMode = "single";
    phiAse->numDevices = 1u;
    phiAse->devices = hase::transport::Array<std::uint64_t>{};
    auto solver = std::make_shared<TimeIntegrationSolver>();
    solver->name = "explicit-euler";

    Simulation simulation;
    simulation.timeStep = 1.0e-6;
    simulation.simulationSteps = 1u;
    simulation.executionMode = "autonomous";
    simulation.outputFields = {"beta_volume"};
    simulation.opticalComponents = components;
    simulation.gainMedium = gainMedium;
    simulation.exteriorSurface = std::make_shared<Domain>();
    simulation.excitationState = excitation;
    simulation.phiAse = phiAse;
    simulation.timeIntegrationSolver = solver;

    auto state = prepareSimulation(simulation);
    CHECK(state.trace.numberOfCells == 2u);
    CHECK(state.trace.numberOfMeshPoints == 8u);
    CHECK(state.trace.betaVolume == std::vector<double>{0.2, 0.4});
    CHECK(state.trace.cellPointIndices == std::vector<unsigned>{0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u});
    CHECK(state.trace.cellCenters == std::vector<double>{0.25, 10.25, 0.25, 0.25, 0.25, 0.25});
    CHECK(state.trace.numberOfMaterials == 2u);
    CHECK(state.trace.cellMaterialIds == std::vector<unsigned>{0u, 1u});
    CHECK(state.trace.materialCrossSectionOffsets == std::vector<unsigned>{0u, 1u, 2u});
    CHECK(state.trace.crossSectionWavelengths == std::vector<double>{1030.0e-9, 1064.0e-9});
    CHECK(state.trace.crossSectionAbsorption == std::vector<double>{1.0e-25, 3.0e-25});
    CHECK(state.trace.crossSectionEmission == std::vector<double>{2.0e-24, 4.0e-24});
    CHECK(state.trace.materialActiveIonDensities == std::vector<double>{2.76e26, 3.10e26});
    CHECK(state.trace.materialBulkAttenuations == std::vector<double>{0.5, 1.5});
}
