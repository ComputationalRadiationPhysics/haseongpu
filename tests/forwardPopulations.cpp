#include "include/populationCube.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <core/forwardPhiAseEvaluator.hpp>

#include <algorithm>
#include <cmath>

namespace
{
    using TestBackends = std::decay_t<
        decltype(alpaka::onHost::allBackends(alpaka::onHost::enabledDeviceSpecs, alpaka::exec::enabledExecutors))>;

    hase::data::AseDomainGraph populationGraph(hase::data::TraceData const& mesh)
    {
        auto graph = hase::test::populationSources(mesh);
        graph.domains.resize(1u);
        graph.domains[0u].id = 0u;
        graph.domains[0u].trace = mesh;
        graph.domains[0u].localToGlobalCells = graph.domainGlobalCells;
        graph.globalCellDomains.assign(mesh.numberOfCells, 0u);
        auto const faces = mesh.numberOfCells * mesh.numberOfFacesPerCell;
        graph.globalBoundaryTargetDomains.assign(faces, hase::data::invalidDomainId);
        graph.globalBoundaryTargetCells.assign(faces, 0u);
        graph.globalBoundaryTargetFaces.assign(faces, 0u);
        graph.globalBoundaryReflectivities.assign(faces, 0.0f);
        graph.globalBoundarySourceRefractiveIndices.assign(faces, 1.0f);
        graph.globalBoundaryTargetRefractiveIndices.assign(faces, 1.0f);
        return graph;
    }

    // Three slabs with two emitting domains and one unpumped transport domain.
    hase::data::AseDomainGraph slabPopulationGraph(hase::data::TraceData& mesh)
    {
        auto graph = populationGraph(mesh);
        graph.domains.resize(3u);
        graph.domainCellOffsets = {0u};
        graph.domainGlobalCells.clear();
        for(std::uint32_t domain = 0u; domain < 3u; ++domain)
        {
            auto& entry = graph.domains[domain];
            entry.id = domain;
            entry.localToGlobalCells.clear();
            for(std::uint32_t cell = 0u; cell < mesh.numberOfCells; ++cell)
                if(static_cast<std::uint32_t>(3.0 * mesh.cellCenters[cell]) == domain)
                {
                    entry.localToGlobalCells.push_back(cell);
                    graph.domainGlobalCells.push_back(cell);
                    graph.globalCellDomains[cell] = domain;
                    mesh.betaVolume[cell] = domain == 2u ? 0.0 : (domain == 0u ? 0.5 : 0.25);
                }
            entry.trace = mesh;
            graph.domainCellOffsets.push_back(static_cast<std::uint32_t>(graph.domainGlobalCells.size()));
        }
        mesh.rebuildStaticPrefixes();
        graph.domainSourceStrengthTotals.assign(3u, 0.0);
        graph.domainSourceStrengthPrefix.clear();
        for(auto const cell : graph.domainGlobalCells)
        {
            double const strength = mesh.betaVolume[cell] * mesh.cellVolumes[cell];
            auto& total = graph.domainSourceStrengthTotals[graph.globalCellDomains[cell]];
            graph.domainSourceStrengthPrefix.push_back(total += strength);
        }
        for(std::uint32_t cell = 0u; cell < mesh.numberOfCells; ++cell)
            for(std::uint32_t face = 0u; face < 4u; ++face)
            {
                auto const index = 4u * cell + face;
                auto const neighbor = mesh.cellNeighborCells[index];
                if(neighbor < 0 || graph.globalCellDomains[cell] == graph.globalCellDomains[neighbor])
                    continue;
                graph.globalBoundaryTargetDomains[index] = graph.globalCellDomains[neighbor];
                graph.globalBoundaryTargetCells[index] = static_cast<std::uint32_t>(neighbor);
                graph.globalBoundaryTargetFaces[index]
                    = static_cast<std::uint32_t>(mesh.cellNeighborLocalFaces[index]);
                mesh.cellNeighborCells[index] = -1;
                mesh.cellFaceBoundaries[index] = 1;
            }
        return graph;
    }
} // namespace

TEMPLATE_LIST_TEST_CASE(
    "ray population estimates are independent of worker partitioning",
    "[forward][populations][correctness]",
    TestBackends)
{
    auto selector = alpaka::onHost::makeDeviceSelector(TestType::makeDict());
    if(!selector.isAvailable())
        SKIP("Requested test backend has no available device");
    auto device = selector.makeDevice(0u);
    auto const executor = alpaka::getExecutor(TestType::makeDict());
    for(auto const reflectivity : {0.0f, 0.25f})
        for(auto const populations : {3u, 8u, 12u})
        {
            auto mesh = hase::test::populationCube(2u, reflectivity);
            auto graph = populationGraph(mesh);
            hase::core::AseTraceControls controls;
            controls.numIndependentRayPopulations = populations;
            controls.forwardRayCount = 8195u;
            controls.minRays = controls.maxRays = controls.forwardRayCount;
            controls.useReflections = reflectivity > 0.0f;
            controls.boundaryMaxPasses = 40u;
            controls.reflectionTolerance = 1.0e-4;
            hase::data::PhiAseResult reference;
            for(auto const workers : {1u, 2u, 5u, 13u})
            {
                CAPTURE(reflectivity, populations, workers);
                hase::core::ForwardPhiAseContext
                    context(std::vector(workers, device), executor, controls, mesh, graph);
                hase::core::ExecutionPolicy compute(
                    1u,
                    0u,
                    workers,
                    0u,
                    "population-test",
                    hase::core::ParallelMode::SINGLE,
                    false,
                    {0u},
                    0u,
                    mesh.numberOfCells,
                    137u);
                hase::data::PhiAseResult result;
                auto const evaluation = context.evaluate(controls, compute, mesh, context.primaryBetaVolume(), result);
                REQUIRE(evaluation.rayCount == controls.forwardRayCount);
                REQUIRE(evaluation.adaptiveLaunches == 1u);
                REQUIRE(
                    result.boundaryStatus
                    == (reflectivity == 0.0f ? hase::data::BoundaryStatus::disabled
                                             : hase::data::BoundaryStatus::converged));
                REQUIRE(result.boundaryMaxPasses == (reflectivity == 0.0f ? 0u : controls.boundaryMaxPasses));
                if(workers == 1u)
                    reference = result;
                REQUIRE(result.phiAse.size() == mesh.numberOfCells);
                for(std::uint32_t cell = 0u; cell < mesh.numberOfCells; ++cell)
                {
                    CHECK(result.phiAse[cell] == Catch::Approx(reference.phiAse[cell]).epsilon(1.0e-6));
                    CHECK(result.standardError[cell] == Catch::Approx(reference.standardError[cell]).epsilon(1.0e-10));
                    CHECK(
                        result.relativeStandardError[cell]
                        == Catch::Approx(reference.relativeStandardError[cell]).epsilon(1.0e-10));
                }
                CHECK(result.boundaryPasses == reference.boundaryPasses);
            }
        }
}

TEST_CASE("population plan covers every ray independently of execution chunk size", "[forward][populations]")
{
    using namespace hase::core;
    std::vector<DomainQuota> quotas{{0u, 101u, 1.0, 0.0}, {1u, 53u, 2.0, 0.0}, {2u, 0u, 0.0, 0.0}};
    std::vector<std::uint32_t> counts{101u, 53u, 0u};
    for(auto const populations : {1u, 3u, 8u})
        for(auto const chunk : {1u, 7u, 64u})
        {
            auto const plan = makeForwardPopulationBatches(counts, quotas, 3.0, populations, chunk);
            std::vector<std::vector<std::uint32_t>> seen(populations);
            std::uint32_t total = 0u;
            for(auto const& work : plan)
            {
                CHECK(work.domainId != 2u);
                CHECK(work.rayCount <= chunk);
                for(std::uint32_t i = 0u; i < work.rayCount; ++i)
                    seen.at(work.rayPopulationId).push_back(work.populationOffset + work.rayOffset + i);
                total += work.rayCount;
            }
            CHECK(total == 154u);
            for(auto const& rays : seen)
                for(std::size_t i = 0u; i < rays.size(); ++i)
                    CHECK(rays[i] == i);
        }
    CHECK_THROWS_AS(makeForwardPopulationBatches(counts, quotas, 3.0, 54u, 16u), std::invalid_argument);
}

TEMPLATE_LIST_TEST_CASE(
    "logical SRM batches preserve statistics across worker ownership",
    "[forward][populations][srm]",
    TestBackends)
{
    auto selector = alpaka::onHost::makeDeviceSelector(TestType::makeDict());
    if(!selector.isAvailable())
        SKIP("Requested test backend has no available device");
    auto device = selector.makeDevice(0u);
    auto const executor = alpaka::getExecutor(TestType::makeDict());
    for(auto const* position : {"exact", "centroid"})
        for(bool const multipleDomains : {false, true})
            DYNAMIC_SECTION("SRM " << position << ", multiple domains: " << multipleDomains)
            {
                auto mesh = hase::test::populationCube(multipleDomains ? 3u : 1u, 0.25f);
                auto graph = multipleDomains ? slabPopulationGraph(mesh) : populationGraph(mesh);
                hase::core::AseTraceControls controls;
                controls.reflectionMode = "srm";
                controls.srmPositionMode = position;
                controls.surfaceReservoirSize = 4u;
                controls.domainCount = graph.domains.size();
                controls.numIndependentRayPopulations = multipleDomains ? 8u : 2u;
                // Two logical batches per population in the single-domain case, including a small tail.
                controls.forwardRayCount = multipleDomains ? 8195u : 2u * hase::core::maxLogicalSrmBatchRays + 3u;
                controls.minRays = controls.maxRays = controls.forwardRayCount;
                controls.useReflections = true;
                controls.boundaryMaxPasses = 40u;
                controls.reflectionTolerance = 1.0e-4;
                hase::data::PhiAseResult reference;
                for(auto const workers : {1u, 2u, 5u})
                {
                    CAPTURE(position, multipleDomains, workers);
                    controls.enableDiagnostics = workers != 2u;
                    hase::core::ForwardPhiAseContext
                        context(std::vector(workers, device), executor, controls, mesh, graph);
                    hase::core::ExecutionPolicy compute(
                        1u,
                        0u,
                        workers,
                        0u,
                        "srm-population-test",
                        hase::core::ParallelMode::SINGLE,
                        false,
                        {0u},
                        0u,
                        mesh.numberOfCells,
                        137u);
                    hase::data::PhiAseResult result;
                    auto const evaluation
                        = context.evaluate(controls, compute, mesh, context.primaryBetaVolume(), result);
                    REQUIRE(evaluation.rayCount == controls.forwardRayCount);
                    REQUIRE(result.boundaryStatus == hase::data::BoundaryStatus::converged);
                    REQUIRE(result.phiAse.size() == mesh.numberOfCells);
                    if(workers == 1u)
                        reference = result;
                    for(std::uint32_t cell = 0u; cell < mesh.numberOfCells; ++cell)
                    {
                        CHECK(result.phiAse[cell] > 0.0f); // Includes every unpumped receiving cell.
                        CHECK(result.phiAse[cell] == Catch::Approx(reference.phiAse[cell]).epsilon(1.0e-6));
                        CHECK(
                            result.standardError[cell]
                            == Catch::Approx(reference.standardError[cell]).epsilon(1.0e-10));
                        CHECK(
                            result.relativeStandardError[cell]
                            == Catch::Approx(reference.relativeStandardError[cell]).epsilon(1.0e-10));
                        if(!controls.enableDiagnostics)
                            CHECK(result.totalRays[cell] == 0u);
                    }
                    CHECK(result.boundaryPasses == reference.boundaryPasses);
                }
            }
}

TEMPLATE_LIST_TEST_CASE(
    "zero-source adaptive budgets are independent of cell visits",
    "[forward][populations][zero-source]",
    TestBackends)
{
    auto selector = alpaka::onHost::makeDeviceSelector(TestType::makeDict());
    if(!selector.isAvailable())
        SKIP("Requested test backend has no available device");
    auto device = selector.makeDevice(0u);
    auto const executor = alpaka::getExecutor(TestType::makeDict());
    auto modes = std::vector{hase::core::ParallelMode::SINGLE};
#if defined(MPI_FOUND) && !defined(DISABLE_MPI)
    modes.push_back(hase::core::ParallelMode::MPI);
#endif
    for(auto const& mode : modes)
        for(auto const* reflectionMode : {"direct", "srm"})
        {
            CAPTURE(mode, reflectionMode);
            auto mesh = hase::test::populationCube(1u, 0.0f);
            std::fill(mesh.betaVolume.begin(), mesh.betaVolume.end(), 0.0);
            mesh.rebuildStaticPrefixes();
            auto graph = populationGraph(mesh);
            hase::core::AseTraceControls controls;
            controls.minRays = 1u;
            controls.maxRays = 8u;
            controls.forwardRayCount = 0u;
            controls.numIndependentRayPopulations = 8u;
            controls.relativeStandardErrorThreshold = 0.0;
            controls.enableDiagnostics = true;
            controls.useReflections = false;
            controls.reflectionMode = reflectionMode;
            hase::core::ForwardPhiAseContext context(std::vector{device}, executor, controls, mesh, graph);
            hase::core::ExecutionPolicy
                compute(1u, 3u, 1u, 0u, "zero-source-test", mode, false, {0u}, 0u, mesh.numberOfCells, 1234u);
            hase::data::PhiAseResult result;
            auto const evaluation = context.evaluate(controls, compute, mesh, context.primaryBetaVolume(), result);
            CHECK(evaluation.rayCount == 8u);
            CHECK(evaluation.adaptiveLaunches == 4u);
            REQUIRE(result.phiAse.size() == mesh.numberOfCells);
            for(std::uint32_t cell = 0u; cell < mesh.numberOfCells; ++cell)
            {
                CHECK(result.phiAse[cell] == 0.0);
                CHECK(result.totalRays[cell] == 0u);
                CHECK(std::isnan(result.relativeStandardError[cell]));
            }
        }
}

TEMPLATE_LIST_TEST_CASE(
    "an unpumped domain retains exactly zero source strength",
    "[forward][populations][source]",
    TestBackends)
{
    auto selector = alpaka::onHost::makeDeviceSelector(TestType::makeDict());
    if(!selector.isAvailable())
        SKIP("Requested test backend has no available device");
    auto device = selector.makeDevice(0u);
    auto const executor = alpaka::getExecutor(TestType::makeDict());
    auto mesh = hase::test::populationCube(3u, 0.25f);
    auto graph = slabPopulationGraph(mesh);
    REQUIRE(graph.domainSourceStrengthTotals[2u] == 0.0);
    auto resident = mesh.makeResident(device);
    hase::core::ResidentAseDomainSources sources(device, graph);
    hase::alpakaUtils::DevBundle bundle(device, executor);
    auto queue = device.makeQueue(alpaka::queueKind::nonBlocking);
    resident.toDevice(queue);
    sources.toDevice(queue);
    sources.rebuild(bundle, queue, resident.view());
    auto const totals = sources.downloadSourceStrengthTotals(queue);
    CAPTURE(totals);
    REQUIRE(totals.size() == 3u);
    CHECK(totals[0u] > 0.0);
    CHECK(totals[1u] > 0.0);
    CHECK(totals[2u] == 0.0);
}
