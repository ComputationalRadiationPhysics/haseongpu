#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/memory.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <core/mesh.hpp>
#include <core/simulationRunControl.hpp>
#include <kernels/generalPump.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
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
        return hase::core::HostMesh{
            {0u, 1u, 2u, 3u},
            {10u},
            {1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2},
            {-1, -1, -1, -1},
            {-1, -1, -1, -1},
            {7, 8, 9, 10},
            {1.0f / 6.0f},
            // Structure-of-arrays point coordinates.
            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},
            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},
            {0.25, 0.25, 0.25},
            {0.0},
            {0.0, 0.0, 0.0, 0.0},
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
            true};
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
} // namespace

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
    CHECK(first.originX == repeated.originX);
    CHECK(first.originY == repeated.originY);
    CHECK(first.wavelength == repeated.wavelength);
    CHECK(std::accumulate(first.power.begin(), first.power.end(), 0.0) == Catch::Approx(source.totalPower));
    for(std::size_t ray = 0u; ray < first.size(); ++ray)
    {
        CHECK(first.originZ[ray] == Catch::Approx(0.0));
        CHECK(first.originX[ray] >= 0.0);
        CHECK(first.originY[ray] >= 0.0);
        CHECK(first.originX[ray] + first.originY[ray] <= 1.0);
        CHECK(first.directionX[ray] == Catch::Approx(0.0).margin(1.0e-14));
        CHECK(first.directionY[ray] == Catch::Approx(0.0).margin(1.0e-14));
        CHECK(first.directionZ[ray] == Catch::Approx(1.0));
        CHECK(first.cell[ray] == 0u);
        CHECK(first.forbiddenFace[ray] == 3);
        if(first.wavelength[ray] == source.wavelengths[0])
        {
            CHECK(first.sigmaAbsorption[ray] == source.sigmaAbsorption[0]);
            CHECK(first.sigmaEmission[ray] == source.sigmaEmission[0]);
        }
        else
        {
            CHECK(first.wavelength[ray] == source.wavelengths[1]);
            CHECK(first.sigmaAbsorption[ray] == source.sigmaAbsorption[1]);
            CHECK(first.sigmaEmission[ray] == source.sigmaEmission[1]);
        }
    }
}

TEST_CASE("planar pump relay retroreflects, transmits, and vignettes", "[pump][relay]")
{
    auto const mesh = makeSingleTetMesh();
    hase::kernels::PumpRayBatch exiting;
    exiting.originX = {0.2};
    exiting.originY = {0.3};
    exiting.originZ = {0.0};
    exiting.directionX = {0.0};
    exiting.directionY = {0.0};
    exiting.directionZ = {-1.0};
    exiting.power = {5.0};
    exiting.wavelength = {940e-9};
    exiting.sigmaAbsorption = {1.0e-20};
    exiting.sigmaEmission = {0.0};
    exiting.cell = {0u};
    exiting.forbiddenFace = {3};
    exiting.exitFace = {3};

    auto const frame = hase::kernels::makeRelayFrame(mesh, {10});
    REQUIRE(frame.faces.size() == 1u);
    CHECK(hase::kernels::pointInTriangle({0.2, 0.3, 0.0}, frame.faces.front(), frame.u, frame.v));

    hase::core::PumpRelayParameters relay;
    relay.exitSurfaces = {10};
    relay.entrySurfaces = {10};
    relay.transmission = 0.4;
    auto returned = hase::kernels::applyPumpRelay(mesh, exiting, relay);
    REQUIRE(returned.size() == 1u);
    CHECK(returned.originX[0] == Catch::Approx(0.2));
    CHECK(returned.originY[0] == Catch::Approx(0.3));
    CHECK(returned.originZ[0] == Catch::Approx(0.0).margin(1.0e-14));
    CHECK(returned.directionZ[0] == Catch::Approx(1.0));
    CHECK(returned.power[0] == Catch::Approx(2.0));
    CHECK(returned.cell[0] == 0u);
    CHECK(returned.forbiddenFace[0] == 3);

    relay.offset[0] = 2.0;
    CHECK(hase::kernels::applyPumpRelay(mesh, exiting, relay).size() == 0u);
}

TEMPLATE_LIST_TEST_CASE(
    "general pump device trace matches Beer-Lambert power and conservative deposition",
    "[pump][backend]",
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
    auto cellIntegral = hase::alpakaUtils::toDevice(queue, std::vector<double>{0.0});
    auto sampleIntegral = hase::alpakaUtils::toDevice(queue, std::vector<double>(4u, 0.0));

    hase::kernels::PumpRayBatch ray;
    ray.originX = {0.2};
    ray.originY = {0.2};
    ray.originZ = {0.0};
    ray.directionX = {0.0};
    ray.directionY = {0.0};
    ray.directionZ = {1.0};
    ray.power = {10.0};
    ray.wavelength = {940e-9};
    ray.sigmaAbsorption = {1.0e-20};
    ray.sigmaEmission = {0.0};
    ray.cell = {0u};
    ray.forbiddenFace = {3};
    ray.exitFace = {-1};

    auto traced = hase::kernels::tracePumpBatch(
        devBundle,
        queue,
        deviceMesh.toView(),
        betaVolume,
        cellIntegral,
        sampleIntegral,
        std::move(ray));

    constexpr double pathLength = 0.6;
    double const expectedPower = 10.0 * std::exp(-pathLength);
    double const expectedIntegral = (10.0 - expectedPower) * 940e-9 / (6.62607015e-34 * 299792458.0 * 1.0e20);
    REQUIRE(traced.size() == 1u);
    CHECK(traced.power[0] == Catch::Approx(expectedPower).epsilon(2.0e-6));
    CHECK(traced.exitFace[0] == 0);

    auto const cellValues = copyDoubleBuffer(queue, cellIntegral);
    auto const sampleValues = copyDoubleBuffer(queue, sampleIntegral);
    REQUIRE(cellValues.size() == 1u);
    CHECK(cellValues[0] == Catch::Approx(expectedIntegral).epsilon(2.0e-6));
    CHECK(
        std::accumulate(sampleValues.begin(), sampleValues.end(), 0.0)
        == Catch::Approx(cellValues[0]).epsilon(2.0e-6));
    CHECK(std::ranges::all_of(sampleValues, [](double const value) { return value >= 0.0; }));
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
    for(std::size_t ray = 0u; ray < rays.size(); ++ray)
    {
        double const norm = std::sqrt(
            rays.directionX[ray] * rays.directionX[ray] + rays.directionY[ray] * rays.directionY[ray]
            + rays.directionZ[ray] * rays.directionZ[ray]);
        CHECK(norm == Catch::Approx(1.0));
        CHECK(rays.directionZ[ray] == Catch::Approx(std::cos(polar)));
    }
}

TEMPLATE_LIST_TEST_CASE(
    "general pump orchestration conserves cell-to-sample deposition",
    "[pump][backend]",
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
    auto cellIntegral = hase::alpakaUtils::toDevice(queue, std::vector<double>{0.0});
    double const lumpedShare = static_cast<double>(mesh.cellVolumes[0]) / 4.0;
    auto lumpedVolume = hase::alpakaUtils::toDevice(queue, std::vector<double>(4u, lumpedShare));
    auto sampleRate = hase::alpakaUtils::toDevice(queue, std::vector<double>(4u, 0.0));

    auto source = uniformSource(10u);
    source.wavelengths = {940e-9};
    source.spectralWeights = {1.0};
    source.sigmaAbsorption = {1.0e-20};
    source.sigmaEmission = {0.0};
    hase::core::PumpParameters pump;
    pump.rayCount = 1024u;
    pump.rngSeed = 42u;
    pump.sources = {source};

    hase::kernels::enqueueGeneralPump(
        devBundle,
        queue,
        mesh,
        deviceMesh.toView(),
        pump,
        betaVolume,
        cellIntegral,
        lumpedVolume,
        sampleRate);

    auto const cells = copyDoubleBuffer(queue, cellIntegral);
    auto const samples = copyDoubleBuffer(queue, sampleRate);
    REQUIRE(cells.size() == 1u);
    REQUIRE(samples.size() == 4u);
    CHECK(cells[0] > 0.0);
    CHECK(std::ranges::all_of(samples, [](double const value) { return std::isfinite(value) && value >= 0.0; }));
    double const integratedSamples = std::accumulate(samples.begin(), samples.end(), 0.0) * lumpedShare;
    CHECK(integratedSamples == Catch::Approx(cells[0]).epsilon(2.0e-6));
}
