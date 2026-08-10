#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/memory.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <core/mesh.hpp>
#include <core/simulationRunControl.hpp>
#include <kernels/generalPump.hpp>
#include <kernels/timeIntegrationUpdateKernels.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <numeric>
#include <type_traits>
#include <vector>

namespace
{
    using TestBackends = std::decay_t<
        decltype(alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors))>;

    hase::core::HostMesh makeSingleTetMesh()
    {
        // Unit right tetrahedron. Local face i is opposite local vertex i.
        std::vector<double> const meshPoints = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
        std::vector<double> const samplePoints = {0.25, 0.25, 0.25};
        return hase::core::HostMesh{
            {0u, 1u, 2u, 3u},
            {10u},
            {1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2},
            {-1, -1, -1, -1},
            {-1, -1, -1, -1},
            {7, 8, 9, 10},
            {1.0f / 6.0f},
            meshPoints,
            samplePoints,
            {0.25, 0.25, 0.25},
            {0.0},
            {10u},
            {1.0f, 1.0f, 1.0f, 1.0f},
            {0.0f, 0.0f},
            std::vector<float>(11u, 0.0f),
            std::vector<float>(11u, 1.0f),
            std::vector<float>(11u, 1.0f),
            1.0e20f,
            1.0f,
            99u,
            0.0,
            4u,
            1u,
            0.0f,
            false};
    }

    hase::core::HostMesh makeTwoTetMesh()
    {
        // Two unit-height tetrahedra share the triangle (0, 1, 2).
        std::vector<double> const meshPoints
            = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0};
        std::vector<double> const centers = {0.25, 0.25, 0.25, 0.25, 0.25, -0.25};
        return hase::core::HostMesh{
            {0u, 1u, 2u, 3u, 0u, 2u, 1u, 4u},
            {10u, 10u},
            {1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2, 2, 1, 4, 0, 1, 4, 0, 2, 4, 0, 2, 1},
            {-1, -1, -1, 1, -1, -1, -1, 0},
            {-1, -1, -1, 3, -1, -1, -1, 3},
            {7, 8, 9, 0, 7, 8, 9, 0},
            {1.0f / 6.0f, 1.0f / 6.0f},
            meshPoints,
            centers,
            centers,
            {0.0, 0.0},
            {10u, 10u},
            {1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
            {0.0f, 0.0f},
            std::vector<float>(11u, 0.0f),
            std::vector<float>(11u, 1.0f),
            std::vector<float>(11u, 1.0f),
            1.0e20f,
            1.0f,
            99u,
            0.0,
            5u,
            1u,
            0.0f,
            false};
    }

    hase::core::PumpSourceParameters uniformSource(unsigned const surface)
    {
        hase::core::PumpSourceParameters source;
        source.surfaces = {static_cast<int>(surface)};
        source.totalPower = 80.0;
        source.wavelengths = {900e-9, 1000e-9};
        source.spectralWeights = {1.0, 3.0};
        source.sigmaAbsorption = {1.0e-20, 2.0e-20};
        source.sigmaEmission = {3.0e-20, 4.0e-20};
        source.polarAngles = {0.0};
        source.azimuthalAngles = {0.0};
        source.angularWeights = {1.0};
        return source;
    }

    template<typename T_Buffer>
    std::vector<double> copyDoubleBuffer(auto const& queue, T_Buffer const& deviceBuffer)
    {
        auto host = alpaka::onHost::allocHostLike(deviceBuffer);
        alpaka::onHost::memcpy(queue, host, deviceBuffer);
        alpaka::onHost::wait(queue);
        auto const size = static_cast<std::size_t>(host.getMdSpan().getExtents().x());
        return {alpaka::onHost::data(host), alpaka::onHost::data(host) + size};
    }

    template<typename T_Buffer>
    std::vector<unsigned> copyUnsignedBuffer(auto const& queue, T_Buffer const& deviceBuffer)
    {
        auto host = alpaka::onHost::allocHostLike(deviceBuffer);
        alpaka::onHost::memcpy(queue, host, deviceBuffer);
        alpaka::onHost::wait(queue);
        auto const size = static_cast<std::size_t>(host.getMdSpan().getExtents().x());
        return {alpaka::onHost::data(host), alpaka::onHost::data(host) + size};
    }
} // namespace

TEST_CASE("lumped pump vertex volumes exclude cladding cells", "[pump][projection]")
{
    auto mesh = makeTwoTetMesh();
    double const cellShare = static_cast<double>(mesh.cellVolumes[0]) / 4.0;
    auto const gainVolumes = hase::kernels::makeLumpedGainVertexVolumes(mesh);
    CHECK(gainVolumes[0] == Catch::Approx(2.0 * cellShare));
    CHECK(gainVolumes[3] == Catch::Approx(cellShare));
    CHECK(gainVolumes[4] == Catch::Approx(cellShare));

    mesh.claddingCellTypes[1] = mesh.claddingNumber;
    auto const interfaceVolumes = hase::kernels::makeLumpedGainVertexVolumes(mesh);
    CHECK(interfaceVolumes[0] == Catch::Approx(cellShare));
    CHECK(interfaceVolumes[3] == Catch::Approx(cellShare));
    CHECK(interfaceVolumes[4] == 0.0);
}

TEMPLATE_LIST_TEST_CASE(
    "vertex pump projection smooths shared gain cells and preserves interfaces",
    "[pump][backend][projection]",
    TestBackends)
{
    auto const backend = TestType::makeDict();
    auto deviceSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
    if(!deviceSelector.isAvailable())
    {
        SUCCEED("No device available for " << backend[alpaka::object::deviceSpec].getName());
        return;
    }
    auto device = deviceSelector.makeDevice(0);
    auto const executor = backend[alpaka::object::exec];
    auto queue = device.makeQueue(alpaka::queueKind::blocking);
    hase::alpakaUtils::DevBundle devBundle(device, executor);

    auto project = [&](hase::core::HostMesh& mesh)
    {
        auto deviceMesh = mesh.toDevice(device);
        std::vector<double> vertexValues(mesh.numberOfMeshPoints, 0.0);
        vertexValues[0] = 1.0;
        auto vertexIntegral = hase::alpakaUtils::toDevice(queue, vertexValues);
        auto lumpedVertexVolume = hase::alpakaUtils::toDevice(queue, hase::kernels::makeLumpedGainVertexVolumes(mesh));
        auto cellRate = hase::alpakaUtils::toDevice(queue, std::vector<double>(mesh.numberOfCells, 0.0));
        hase::kernels::enqueueProjectVertexPumpRateToCells(
            devBundle,
            queue,
            deviceMesh.toView(),
            vertexIntegral,
            lumpedVertexVolume,
            cellRate);
        return copyDoubleBuffer(queue, cellRate);
    };

    auto mesh = makeTwoTetMesh();
    auto const smoothedRate = project(mesh);
    REQUIRE(smoothedRate.size() == 2u);
    CHECK(smoothedRate[0] == Catch::Approx(smoothedRate[1]));
    double const smoothedIntegral = smoothedRate[0] * mesh.cellVolumes[0] + smoothedRate[1] * mesh.cellVolumes[1];
    CHECK(smoothedIntegral == Catch::Approx(1.0));

    mesh.claddingCellTypes[1] = mesh.claddingNumber;
    auto const interfaceRate = project(mesh);
    REQUIRE(interfaceRate.size() == 2u);
    CHECK(interfaceRate[1] == 0.0);
    CHECK(interfaceRate[0] * mesh.cellVolumes[0] == Catch::Approx(1.0));
}

TEST_CASE("ray walk SRM policies select boundary-position storage at compile time", "[ray][policy]")
{
    using namespace hase::kernels::forward::ray;
    using AseStorage = SrmPositionStorage<typename std::remove_cvref_t<ALPAKA_TYPEOF(aseSrmPolicy)>::PositionPolicy>;
    using PumpStorage = SrmPositionStorage<typename std::remove_cvref_t<ALPAKA_TYPEOF(pumpSrmPolicy)>::PositionPolicy>;

    STATIC_CHECK_FALSE(std::derived_from<AseStorage, BarycentricSrmPositionStorage>);
    STATIC_CHECK(std::derived_from<PumpStorage, BarycentricSrmPositionStorage>);
    STATIC_CHECK(concepts::BoundaryBehaviour<ALPAKA_TYPEOF(aseSrmPolicy)>);
    STATIC_CHECK(SrmBoundaryBehaviour<ALPAKA_TYPEOF(pumpSrmPolicy)>);
    STATIC_CHECK_FALSE(std::remove_cvref_t<ALPAKA_TYPEOF(aseSrmPolicy)>::PositionPolicy::storesBarycentric);
    STATIC_CHECK(std::remove_cvref_t<ALPAKA_TYPEOF(pumpSrmPolicy)>::PositionPolicy::storesBarycentric);
}

TEST_CASE("general pump samples tagged faces deterministically with conserved source power", "[pump][source]")
{
    auto mesh = makeSingleTetMesh();
    auto const faces = hase::kernels::pumpBoundaryFaces(mesh, {10});
    REQUIRE(faces.size() == 1u);
    CHECK(faces.front().cell == 0u);
    CHECK(faces.front().localFace == 3u);
    CHECK(faces.front().area == Catch::Approx(0.5));
    CHECK(faces.front().normal.z == Catch::Approx(-1.0));

    auto const source = uniformSource(10u);
    auto const first = hase::kernels::samplePumpSource(mesh, source, 80u, 1234u);
    auto const repeated = hase::kernels::samplePumpSource(mesh, source, 80u, 1234u);

    REQUIRE(first.size() == 80u);
    REQUIRE(repeated.size() == first.size());
    double sampledPower = 0.0;
    std::array<unsigned, 4u> entryQuarterVisits{};
    for(std::size_t ray = 0u; ray < first.size(); ++ray)
    {
        auto const& sampled = first[ray];
        auto const& repeatedRay = repeated[ray];
        sampledPower += sampled.power;
        CHECK(sampled.position.x == repeatedRay.position.x);
        CHECK(sampled.position.y == repeatedRay.position.y);
        CHECK(sampled.wavelength == repeatedRay.wavelength);
        CHECK(sampled.position.z == Catch::Approx(0.0));
        CHECK(sampled.position.x >= 0.0);
        CHECK(sampled.position.y >= 0.0);
        CHECK(sampled.position.x + sampled.position.y <= 1.0);
        CHECK(sampled.direction.x == Catch::Approx(0.0).margin(1.0e-14));
        CHECK(sampled.direction.y == Catch::Approx(0.0).margin(1.0e-14));
        CHECK(sampled.direction.z == Catch::Approx(1.0));
        CHECK(sampled.cell == 0u);
        CHECK(sampled.forbiddenFace == 3);
        if(sampled.position.x + sampled.position.y <= 0.5)
            ++entryQuarterVisits[0];
        else if(sampled.position.x >= 0.5)
            ++entryQuarterVisits[1];
        else if(sampled.position.y >= 0.5)
            ++entryQuarterVisits[2];
        else
            ++entryQuarterVisits[3];
        if(sampled.wavelength == source.wavelengths[0])
        {
            CHECK(sampled.sigmaAbsorption == source.sigmaAbsorption[0]);
            CHECK(sampled.sigmaEmission == source.sigmaEmission[0]);
        }
        else
        {
            CHECK(sampled.wavelength == source.wavelengths[1]);
            CHECK(sampled.sigmaAbsorption == source.sigmaAbsorption[1]);
            CHECK(sampled.sigmaEmission == source.sigmaEmission[1]);
        }
    }
    CHECK(sampledPower == Catch::Approx(source.totalPower));
    for(auto const visits : entryQuarterVisits)
        CHECK(visits == 20u);
}

TEST_CASE("general pump stratifies faces by the spatial entry distribution", "[pump][source]")
{
    auto const mesh = makeSingleTetMesh();
    auto source = uniformSource(10u);
    source.surfaces = {7, 8, 9, 10};
    source.profile.kind = 1u;
    source.profile.radiusU = 0.35;
    source.profile.radiusV = 0.55;
    source.profile.exponent = 2.0;
    source.profile.center[0] = 0.15;
    source.profile.center[1] = 0.1;

    auto const faces = hase::kernels::pumpBoundaryFaces(mesh, source.surfaces);
    REQUIRE(faces.size() == 4u);
    std::array<double, 4u> entryWeights{};
    double totalEntryWeight = 0.0;
    auto const regions = hase::kernels::pumpEntryRegions(faces);
    for(auto const& region : regions)
    {
        auto const face = std::ranges::find_if(
            faces,
            [&](auto const& candidate)
            { return candidate.cell == region.cell && candidate.localFace == region.localFace; });
        REQUIRE(face != faces.cend());
        double const weight = hase::kernels::pumpEntryWeight(region, source.profile);
        entryWeights[static_cast<std::size_t>(std::distance(faces.cbegin(), face))] += weight;
        totalEntryWeight += weight;
    }

    constexpr unsigned rayCount = 1000u;
    auto const rays = hase::kernels::samplePumpSource(mesh, source, rayCount, 91u);
    std::array<unsigned, 4u> visits{};
    for(auto const& ray : rays)
    {
        auto const face = std::ranges::find_if(
            faces,
            [&](auto const& candidate)
            { return candidate.cell == ray.cell && static_cast<int>(candidate.localFace) == ray.forbiddenFace; });
        REQUIRE(face != faces.cend());
        ++visits[static_cast<std::size_t>(std::distance(faces.cbegin(), face))];
    }

    for(std::size_t face = 0u; face < faces.size(); ++face)
    {
        double const expected = rayCount * entryWeights[face] / totalEntryWeight;
        CHECK(std::abs(static_cast<double>(visits[face]) - expected) < 1.0);
    }
}

TEST_CASE("pump entry stratification is stable across ray partitions", "[pump][source]")
{
    auto const mesh = makeSingleTetMesh();
    auto source = uniformSource(10u);
    source.surfaces = {7, 8, 9, 10};

    auto const complete = hase::kernels::samplePumpSource(mesh, source, 80u, 17u);
    auto const partition = hase::kernels::samplePumpSource(mesh, source, 80u, 17u, 31u, 23u);
    REQUIRE(partition.size() == 23u);
    for(std::size_t ray = 0u; ray < partition.size(); ++ray)
    {
        auto const& expected = complete[31u + ray];
        auto const& sampled = partition[ray];
        CHECK(sampled.position.x == expected.position.x);
        CHECK(sampled.position.y == expected.position.y);
        CHECK(sampled.position.z == expected.position.z);
        CHECK(sampled.cell == expected.cell);
        CHECK(sampled.forbiddenFace == expected.forbiddenFace);
        CHECK(sampled.wavelength == expected.wavelength);
    }
}

TEST_CASE("pump entry CDF lookup clamps a rounded upper boundary", "[pump][source]")
{
    std::vector<double> const cdf{0.25, 0.5, 1.0};
    CHECK(hase::kernels::pumpEntryRegionForTarget(cdf, 0.0) == 0u);
    CHECK(hase::kernels::pumpEntryRegionForTarget(cdf, 0.25) == 1u);
    CHECK(hase::kernels::pumpEntryRegionForTarget(cdf, 1.0) == 2u);
}

TEST_CASE("pump entry sampling accepts an empty global ray range", "[pump][source]")
{
    auto const rays = hase::kernels::samplePumpSource(makeSingleTetMesh(), uniformSource(10u), 0u, 17u);
    CHECK(rays.empty());
}

TEST_CASE("pump relay preparation encodes physical exit and entry surfaces", "[pump][relay]")
{
    auto const mesh = makeSingleTetMesh();
    hase::core::PumpRelayParameters relay;
    relay.exitSurfaces = {7};
    relay.entrySurfaces = {10};
    relay.transmission = 0.4;
    relay.rotation = 0.25;
    relay.offset[0] = 0.1;
    relay.tilt[1] = -0.2;

    auto const geometry = hase::kernels::preparePumpRelayGeometry(mesh, {relay});
    REQUIRE(geometry.descriptors.size() == 1u);
    REQUIRE(geometry.exitMask.size() == mesh.numberOfCells * mesh.numberOfFacesPerCell);
    REQUIRE(geometry.entryFaceIds.size() == 1u);
    CHECK(geometry.exitMask[0u] == 1u);
    CHECK(geometry.exitMask[3u] == 0u);
    CHECK(geometry.entryFaceIds[0u] == 3u);
    auto const& descriptor = geometry.descriptors.front();
    CHECK(descriptor.transmission == Catch::Approx(0.4));
    CHECK(descriptor.cosine == Catch::Approx(std::cos(0.25)));
    CHECK(descriptor.offsetU == Catch::Approx(0.1));
    CHECK(descriptor.tiltV == Catch::Approx(-0.2));
    CHECK(descriptor.entryFaceEnd - descriptor.entryFaceBegin == 1u);
}

TEMPLATE_LIST_TEST_CASE(
    "general pump reuses device relay target cache without changing deposition",
    "[pump][backend][integration]",
    TestBackends)
{
    auto const backend = TestType::makeDict();
    auto deviceSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
    if(!deviceSelector.isAvailable())
    {
        SUCCEED("No device available for " << backend[alpaka::object::deviceSpec].getName());
        return;
    }
    auto device = deviceSelector.makeDevice(0);
    auto const executor = backend[alpaka::object::exec];
    auto queue = device.makeQueue(alpaka::queueKind::blocking);
    hase::alpakaUtils::DevBundle devBundle(device, executor);

    auto mesh = makeSingleTetMesh();
    auto deviceMesh = mesh.toDevice(device);
    auto betaVolume = hase::alpakaUtils::toDevice(queue, std::vector<double>{0.0});
    auto vertexIntegral = hase::alpakaUtils::toDevice(queue, std::vector<double>(mesh.numberOfMeshPoints, 0.0));

    auto source = uniformSource(10u);
    source.wavelengths = {940e-9};
    source.spectralWeights = {1.0};
    source.sigmaAbsorption = {1.0e-20};
    source.sigmaEmission = {0.0};
    hase::core::PumpRelayParameters relay;
    relay.exitSurfaces = {7};
    relay.entrySurfaces = {7};
    relay.transmission = 0.5;
    source.relays = {relay};
    hase::core::PumpParameters pump;
    pump.rayCount = 512u;
    pump.rngSeed = 17u;
    pump.sources = {source};
    auto prepared
        = hase::kernels::prepareGeneralPumpDeviceSources<decltype(device)>(queue, mesh, pump, 0u, pump.rayCount);
    REQUIRE(prepared.size() == 1u);

    hase::kernels::enqueueGeneralPumpIntegrals(
        devBundle,
        queue,
        deviceMesh.toView(),
        prepared,
        betaVolume,
        vertexIntegral,
        hase::kernels::pumpRelayPolicy);
    auto const firstVertexValues = copyDoubleBuffer(queue, vertexIntegral);
    auto const cacheValues = copyUnsignedBuffer(queue, prepared.front().cacheState());
    REQUIRE(firstVertexValues.size() == mesh.numberOfMeshPoints);
    CHECK(std::accumulate(firstVertexValues.cbegin(), firstVertexValues.cend(), 0.0) > 0.0);
    CHECK(std::count(cacheValues.begin(), cacheValues.end(), 1u) == pump.rayCount);

    hase::kernels::enqueueGeneralPumpIntegrals(
        devBundle,
        queue,
        deviceMesh.toView(),
        prepared,
        betaVolume,
        vertexIntegral,
        hase::kernels::pumpRelayPolicy);
    auto const repeatedVertexValues = copyDoubleBuffer(queue, vertexIntegral);
    REQUIRE(repeatedVertexValues.size() == firstVertexValues.size());
    for(std::size_t point = 0u; point < firstVertexValues.size(); ++point)
        CHECK(repeatedVertexValues[point] == Catch::Approx(firstVertexValues[point]).epsilon(2.0e-6));
}

TEST_CASE("general pump super-Gaussian profile and angular sampling use physical coordinates", "[pump][source]")
{
    hase::core::PumpProfileParameters profile;
    profile.kind = 1u;
    profile.radiusU = 2.0;
    profile.radiusV = 1.0;
    profile.exponent = 2.0;
    profile.center[0] = 0.25;
    CHECK(hase::kernels::pumpProfileWeight(profile, {0.25, 0.0, 0.0}) == Catch::Approx(1.0));
    CHECK(hase::kernels::pumpProfileWeight(profile, {2.25, 0.0, 0.0}) == Catch::Approx(std::exp(-1.0)));
    CHECK(hase::kernels::pumpProfileWeight(profile, {0.25, 1.0, 0.0}) == Catch::Approx(std::exp(-1.0)));

    auto source = uniformSource(10u);
    constexpr double polar = 0.4;
    source.polarAngles = {polar};
    source.azimuthalAngles = {0.7};
    auto const rays = hase::kernels::samplePumpSource(makeSingleTetMesh(), source, 16u, 5u);
    for(auto const& ray : rays)
    {
        double const norm = ray.direction.euclidLength();
        CHECK(norm == Catch::Approx(1.0));
        CHECK(ray.direction.z == Catch::Approx(std::cos(polar)));
    }
}

TEMPLATE_LIST_TEST_CASE(
    "general pump orchestration conserves cell-centered deposition",
    "[pump][backend][integration]",
    TestBackends)
{
    auto const backend = TestType::makeDict();
    auto deviceSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
    if(!deviceSelector.isAvailable())
    {
        SUCCEED("No device available for " << backend[alpaka::object::deviceSpec].getName());
        return;
    }
    auto device = deviceSelector.makeDevice(0);
    auto const executor = backend[alpaka::object::exec];
    auto queue = device.makeQueue(alpaka::queueKind::blocking);
    hase::alpakaUtils::DevBundle devBundle(device, executor);

    auto mesh = makeSingleTetMesh();
    auto deviceMesh = mesh.toDevice(device);
    auto betaVolume = hase::alpakaUtils::toDevice(queue, std::vector<double>{0.0});
    auto vertexIntegral = hase::alpakaUtils::toDevice(queue, std::vector<double>(mesh.numberOfMeshPoints, 0.0));
    auto const lumpedVolumes = hase::kernels::makeLumpedGainVertexVolumes(mesh);
    auto lumpedVertexVolume = hase::alpakaUtils::toDevice(queue, lumpedVolumes);
    auto cellRate = hase::alpakaUtils::toDevice(queue, std::vector<double>{0.0});

    auto source = uniformSource(10u);
    source.wavelengths = {940e-9};
    source.spectralWeights = {1.0};
    source.sigmaAbsorption = {1.0e-20};
    source.sigmaEmission = {0.0};
    hase::core::PumpParameters pump;
    pump.rayCount = 1024u;
    pump.rngSeed = 42u;
    pump.sources = {source};
    auto prepared
        = hase::kernels::prepareGeneralPumpDeviceSources<decltype(device)>(queue, mesh, pump, 0u, pump.rayCount);

    hase::kernels::enqueueGeneralPump(
        devBundle,
        queue,
        deviceMesh.toView(),
        prepared,
        betaVolume,
        vertexIntegral,
        lumpedVertexVolume,
        cellRate);

    auto const vertices = copyDoubleBuffer(queue, vertexIntegral);
    auto const rates = copyDoubleBuffer(queue, cellRate);
    REQUIRE(vertices.size() == mesh.numberOfMeshPoints);
    REQUIRE(rates.size() == 1u);
    double const depositedIntegral = std::accumulate(vertices.cbegin(), vertices.cend(), 0.0);
    CHECK(depositedIntegral > 0.0);
    CHECK(std::isfinite(rates[0]));
    CHECK(rates[0] >= 0.0);
    CHECK(rates[0] * static_cast<double>(mesh.cellVolumes[0]) == Catch::Approx(depositedIntegral).epsilon(2.0e-6));
}

TEMPLATE_LIST_TEST_CASE(
    "device relay and barycentric SRM produce beta fields within five percent",
    "[pump][backend][integration][policy]",
    TestBackends)
{
    auto const backend = TestType::makeDict();
    auto deviceSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
    if(!deviceSelector.isAvailable())
    {
        SUCCEED("No device available for " << backend[alpaka::object::deviceSpec].getName());
        return;
    }
    auto device = deviceSelector.makeDevice(0);
    auto const executor = backend[alpaka::object::exec];
    auto queue = device.makeQueue(alpaka::queueKind::blocking);
    hase::alpakaUtils::DevBundle devBundle(device, executor);
    auto mesh = makeSingleTetMesh();
    auto deviceMesh = mesh.toDevice(device);

    auto source = uniformSource(10u);
    source.wavelengths = {940e-9};
    source.spectralWeights = {1.0};
    source.sigmaAbsorption = {1.0e-20};
    source.sigmaEmission = {0.0};
    hase::core::PumpRelayParameters returnPass;
    returnPass.exitSurfaces = {7};
    returnPass.entrySurfaces = {7};
    returnPass.transmission = 0.8;
    source.relays = {returnPass};
    hase::core::PumpParameters pump;
    pump.rayCount = 4096u;
    pump.rngSeed = 71u;
    pump.sources = {source};

    auto runOneStep = [&](auto boundaryPolicy)
    {
        auto prepared
            = hase::kernels::prepareGeneralPumpDeviceSources<decltype(device)>(queue, mesh, pump, 0u, pump.rayCount);
        auto betaVolume = hase::alpakaUtils::toDevice(queue, std::vector<double>{0.1});
        auto vertexIntegral = hase::alpakaUtils::toDevice(queue, std::vector<double>(mesh.numberOfMeshPoints, 0.0));
        auto lumpedVertexVolume = hase::alpakaUtils::toDevice(queue, hase::kernels::makeLumpedGainVertexVolumes(mesh));
        auto cellRate = hase::alpakaUtils::toDevice(queue, std::vector<double>{0.0});
        hase::kernels::enqueueGeneralPump(
            devBundle,
            queue,
            deviceMesh.toView(),
            prepared,
            betaVolume,
            vertexIntegral,
            lumpedVertexVolume,
            cellRate,
            boundaryPolicy);

        auto beta = hase::alpakaUtils::toDevice(queue, std::vector<double>{0.1});
        auto betaNext = hase::alpakaUtils::toDevice(queue, std::vector<double>{0.0});
        auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfCells});
        queue.enqueue(
            frameSpec,
            alpaka::KernelBundle{hase::kernels::AddScaled{1.0e-3}, deviceMesh.toView(), beta, cellRate, betaNext});
        queue.enqueue(frameSpec, alpaka::KernelBundle{hase::kernels::ClipBeta{}, deviceMesh.toView(), betaNext});
        return copyDoubleBuffer(queue, betaNext);
    };

    auto const relayBeta = runOneStep(hase::kernels::pumpRelayPolicy);
    auto const srmBeta = runOneStep(hase::kernels::pumpSrmBarycentricPolicy);
    REQUIRE(relayBeta.size() == srmBeta.size());
    REQUIRE_FALSE(relayBeta.empty());
    double absoluteDifference = 0.0;
    double relayMagnitude = 0.0;
    for(std::size_t sample = 0u; sample < relayBeta.size(); ++sample)
    {
        CHECK(std::isfinite(relayBeta[sample]));
        CHECK(std::isfinite(srmBeta[sample]));
        CHECK(relayBeta[sample] >= 0.1);
        CHECK(srmBeta[sample] >= 0.1);
        absoluteDifference += std::abs(relayBeta[sample] - srmBeta[sample]);
        relayMagnitude += std::abs(relayBeta[sample] - 0.1);
    }
    REQUIRE(relayMagnitude > 0.0);
    CHECK(absoluteDifference / relayMagnitude < 0.05);
}
