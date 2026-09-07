#include <alpaka/alpaka.hpp>

#include <alpakaUtils/backendNames.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <core/boundaryRayBuffer.hpp>
#include <core/calcForwardPhiAse.hpp>

#include <type_traits>
#include <vector>

namespace
{
    using TestBackends = std::decay_t<
        decltype(alpaka::onHost::allBackends(alpaka::onHost::enabledDeviceSpecs, alpaka::exec::enabledExecutors))>;
} // namespace

TEMPLATE_LIST_TEST_CASE(
    "domain source CDF follows current inversion in grouped cell order",
    "[forward][domain][backend]",
    TestBackends)
{
    auto const backend = TestType::makeDict();
    auto selector = alpaka::onHost::makeDeviceSelector(backend);
    if(!selector.isAvailable())
    {
        SUCCEED("No device available for " << alpaka::onHost::DeviceSpec{backend}.getName());
        return;
    }

    auto device = selector.makeDevice(0);
    auto const executor = alpaka::getExecutor(backend);
    auto queue = device.makeQueue(alpaka::queueKind::nonBlocking);
    hase::alpakaUtils::DevBundle devBundle(device, executor);

    hase::data::TraceData trace;
    trace.numberOfCells = 3u;
    trace.numberOfMaterials = 1u;
    trace.betaVolume = {1.0, 2.0, 3.0};
    trace.cellVolumes = {1.0, 1.0, 1.0};
    trace.cellMaterialIds = {0u, 0u, 0u};
    trace.materialActive = {1u};
    trace.materialActiveIonDensities = {2.0};
    trace.materialFluorescenceLifetimes = {1.0};
    auto residentTrace = trace.makeResident(device);
    residentTrace.toDevice(queue);

    hase::data::AseDomainGraph graph;
    graph.domainCellOffsets = {0u, 2u, 3u};
    graph.domainGlobalCells = {2u, 0u, 1u};
    graph.domainSourceStrengthPrefix.resize(3u, 0.0);
    graph.domainSourceStrengthTotals.resize(2u, 0.0);
    hase::core::ResidentAseDomainSources sources(device, graph);
    sources.toDevice(queue);
    sources.rebuild(devBundle, queue, residentTrace.view());

    auto const totals = sources.downloadSourceStrengthTotals(queue);
    sources.sourceStrengthPrefix.toHost(queue);
    auto const prefix = sources.sourceStrengthPrefix.getHostView();
    REQUIRE(totals.size() == 2u);
    CHECK(totals[0u] == Catch::Approx(8.0));
    CHECK(totals[1u] == Catch::Approx(4.0));
    REQUIRE(prefix.getExtents().product() == 3u);
    CHECK(prefix[0u] == Catch::Approx(6.0));
    CHECK(prefix[1u] == Catch::Approx(8.0));
    // The second domain starts a separate CDF; its physical total is still 4.
    CHECK(prefix[2u] == Catch::Approx(4.0));
}

TEMPLATE_LIST_TEST_CASE("domain CDFs preserve weak zero and empty sources", "[forward][domain][backend]", TestBackends)
{
    auto selector = alpaka::onHost::makeDeviceSelector(TestType::makeDict());
    if(!selector.isAvailable())
        SKIP("Requested test backend has no available device");
    auto device = selector.makeDevice(0u);
    auto queue = device.makeQueue(alpaka::queueKind::nonBlocking);
    hase::alpakaUtils::DevBundle bundle(device, alpaka::getExecutor(TestType::makeDict()));
    hase::data::TraceData trace;
    trace.numberOfCells = 3u;
    trace.numberOfMaterials = 1u;
    trace.betaVolume = {1.0, 1.0e-30, 0.0};
    trace.cellVolumes = {1.0, 1.0, 1.0};
    trace.cellMaterialIds = {0u, 0u, 0u};
    trace.materialActive = {1u};
    trace.materialActiveIonDensities = {2.0};
    trace.materialFluorescenceLifetimes = {1.0};
    auto resident = trace.makeResident(device);
    resident.toDevice(queue);
    hase::data::AseDomainGraph graph;
    graph.domainCellOffsets = {0u, 1u, 2u, 3u, 3u};
    graph.domainGlobalCells = {0u, 1u, 2u};
    graph.domainSourceStrengthPrefix.resize(3u, 0.0);
    graph.domainSourceStrengthTotals.resize(4u, 0.0);
    hase::core::ResidentAseDomainSources sources(device, graph);
    sources.toDevice(queue);
    sources.rebuild(bundle, queue, resident.view());
    auto const totals = sources.downloadSourceStrengthTotals(queue);
    REQUIRE(totals.size() == 4u);
    CHECK(totals[0u] == 2.0);
    CHECK(totals[1u] == 2.0e-30);
    CHECK(totals[2u] == 0.0);
    CHECK(totals[3u] == 0.0);
    sources.sourceStrengthPrefix.toHost(queue);
    auto const prefix = sources.sourceStrengthPrefix.getHostView();
    CHECK(prefix[0u] == 2.0);
    CHECK(prefix[1u] == 2.0e-30);
    CHECK(prefix[2u] == 0.0);
}

TEMPLATE_LIST_TEST_CASE("empty source topology clears domain totals", "[forward][domain][backend]", TestBackends)
{
    auto selector = alpaka::onHost::makeDeviceSelector(TestType::makeDict());
    if(!selector.isAvailable())
        SKIP("Requested test backend has no available device");
    auto device = selector.makeDevice(0u);
    auto queue = device.makeQueue(alpaka::queueKind::nonBlocking);
    hase::alpakaUtils::DevBundle bundle(device, alpaka::getExecutor(TestType::makeDict()));
    hase::data::AseDomainGraph graph;
    graph.domainCellOffsets = {0u, 0u, 0u};
    graph.domainSourceStrengthTotals = {1.0, 2.0};
    hase::core::ResidentAseDomainSources sources(device, graph);
    sources.toDevice(queue);
    sources.rebuild(bundle, queue, hase::data::TraceView{});
    auto const totals = sources.downloadSourceStrengthTotals(queue);
    REQUIRE(totals.size() == 2u);
    CHECK(totals[0u] == 0.0);
    CHECK(totals[1u] == 0.0);
}

TEMPLATE_LIST_TEST_CASE(
    "direct boundary evaluator owns persistent device buffers",
    "[forward][boundary][backend]",
    TestBackends)
{
    auto const backend = TestType::makeDict();
    auto selector = alpaka::onHost::makeDeviceSelector(backend);
    if(!selector.isAvailable())
    {
        SUCCEED("No device available for " << alpaka::onHost::DeviceSpec{backend}.getName());
        return;
    }

    auto device = selector.makeDevice(0);
    auto const executor = alpaka::getExecutor(backend);
    hase::data::TraceData trace;
    hase::core::AseTraceControls controls;
    controls.minRays = 1u;
    controls.maxRays = 1u;
    controls.useReflections = true;
    controls.reflectionMode = "direct";
    hase::core::ForwardPhiAseDeviceContext context(device, executor, controls, trace);

    hase::core::ForwardPhiAseRawResult raw;
    float runtime = -1.0f;
    context.begin(hase::data::TraceView{}, 0u, 17u, 0u, 0.0, controls);
    context.finish(raw, runtime);

    CHECK(raw.rayCount == 0u);
    CHECK(runtime == 0.0f);
}
