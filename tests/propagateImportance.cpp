/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <core/geometry.hpp>
#include <core/mesh.hpp>
#include <kernels/importanceSampling.hpp>
#include <kernels/propagateRay.hpp>

#include <cmath>
#include <iostream>
#include <type_traits>
#include <vector>

using TestApis
    = std::decay_t<decltype(alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors))>;

constexpr unsigned propagationBatchSize = 32u;
constexpr unsigned propagationResultWidth = 3u;

hase::core::HostMesh constructDummyMesh(double betaVolume, float topReflectivity = 1.0f)
{
    hase::core::HostMesh mesh;
    mesh.numberOfTriangles = 1;
    mesh.numberOfLevels = 2;
    mesh.numberOfPoints = 3;
    mesh.thickness = 1.0f;
    mesh.points = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    mesh.triangleCenterX = {1.0 / 3.0};
    mesh.triangleCenterY = {1.0 / 3.0};
    mesh.trianglePointIndices = {0u, 1u, 2u};
    mesh.triangleNormalPoint = {0u, 1u, 2u};
    mesh.triangleNormalsX = {0.0, 1.0, -1.0};
    mesh.triangleNormalsY = {-1.0, 1.0, 0.0};
    mesh.forbiddenEdge = {-1, -1, -1};
    mesh.triangleNeighbors = {-1, -1, -1};
    mesh.triangleSurfaces = {0.5f};
    mesh.betaVolume = {betaVolume};
    mesh.betaCells = {betaVolume};
    mesh.claddingCellTypes = {0u};
    mesh.refractiveIndices = {1.5f, 1.0f, 1.5f, 1.0f};
    mesh.reflectivities = {1.0f, topReflectivity};
    mesh.nTot = 1.0f;
    mesh.crystalTFluo = 1.0f;
    mesh.claddingNumber = 99u;
    mesh.claddingAbsorption = 0.0;
    mesh.calcTotalReflectionAngles();
    return mesh;
}

struct PropagationKernel
{
    // this needs adjustment since only z is the only DOF for the rays - not comparable to real simulation
    ALPAKA_FN_ACC void operator()(auto const& acc, hase::core::DeviceMeshView const mesh, auto result) const
    {
        for(auto [id] : alpaka::onAcc::makeIdxMap(
                acc,
                alpaka::onAcc::worker::threadsInGrid,
                alpaka::IdxRange{hase::alpakaUtils::Vec1D{propagationBatchSize}}))
        {
            unsigned level = 0u;
            unsigned triangle = 0u;
            double const zStart = 0.05 + 0.01 * static_cast<double>(id);
            hase::core::Ray ray{{1.0 / 3.0, 1.0 / 3.0, zStart}, {0.0, 0.0, 0.4}, 0.4f};

            auto const offset = id * propagationResultWidth;
            result[offset] = hase::kernels::propagateRay(ray, &level, &triangle, mesh, 0.0, 0.0);
            result[offset + 1u] = static_cast<double>(level);
            result[offset + 2u] = static_cast<double>(triangle);
        }
    }
};

struct ReflectionKernel
{
    // this needs adjustment since only the z is the only DOF for the rays - not comparable to real simulation
    ALPAKA_FN_ACC void operator()(auto const& acc, hase::core::DeviceMeshView const mesh, auto result) const
    {
        for(auto [id] : alpaka::onAcc::makeIdxMap(
                acc,
                alpaka::onAcc::worker::threadsInGrid,
                alpaka::IdxRange{hase::alpakaUtils::Vec1D{propagationBatchSize}}))
        {
            double const z = 0.1 + 0.01 * static_cast<double>(id);
            hase::core::Point startPoint{1.0 / 3.0, 1.0 / 3.0, z};
            hase::core::Point endPoint{1.0 / 3.0, 1.0 / 3.0, z};

            result[id] = hase::kernels::propagateRayWithReflection(
                startPoint,
                endPoint,
                1u,
                hase::core::TOP_REFLECTION,
                0u,
                0u,
                mesh,
                0.0,
                0.0);
        }
    }
};

template<typename T_Device, typename T_Executor>
std::vector<double> runPropagationKernel(hase::core::HostMesh& hostMesh, T_Device& device, T_Executor const& executor)
{
    auto queue = device.makeQueue();
    auto deviceMesh = hostMesh.toDevice(device);
    auto result = alpaka::onHost::alloc<double>(device, propagationBatchSize * propagationResultWidth);
    queue.enqueue(
        executor,
        alpaka::onHost::FrameSpec{hase::alpakaUtils::Vec1D{4u}, hase::alpakaUtils::Vec1D{8u}},
        alpaka::KernelBundle{PropagationKernel{}, deviceMesh.toView(), result});
    auto hostResult = alpaka::onHost::allocHostLike(result);
    alpaka::onHost::memcpy(queue, hostResult, result);
    alpaka::onHost::wait(queue);

    auto* data = alpaka::onHost::data(hostResult);
    return {data, data + propagationBatchSize * propagationResultWidth};
}

template<typename T_Device, typename T_Executor>
std::vector<double> runReflectionKernel(hase::core::HostMesh& hostMesh, T_Device& device, T_Executor const& executor)
{
    auto queue = device.makeQueue();
    auto deviceMesh = hostMesh.toDevice(device);
    auto result = alpaka::onHost::alloc<double>(device, propagationBatchSize);
    queue.enqueue(
        executor,
        alpaka::onHost::FrameSpec{hase::alpakaUtils::Vec1D{4u}, hase::alpakaUtils::Vec1D{8u}},
        alpaka::KernelBundle{ReflectionKernel{}, deviceMesh.toView(), result});
    auto hostResult = alpaka::onHost::allocHostLike(result);
    alpaka::onHost::memcpy(queue, hostResult, result);
    alpaka::onHost::wait(queue);
    auto* data = alpaka::onHost::data(hostResult);
    return {data, data + propagationBatchSize};
}

TEMPLATE_LIST_TEST_CASE("propagateRay preserves neutral gain inside a prism", "", TestApis)
{
    auto cfg = TestType::makeDict();
    auto deviceSpec = cfg[alpaka::object::deviceSpec];
    auto exec = cfg[alpaka::object::exec];

    auto devSelector = alpaka::onHost::makeDeviceSelector(deviceSpec);
    if(!devSelector.isAvailable())
    {
        std::cout << "No device available for " << deviceSpec.getName() << std::endl;
        return;
    }

    auto device = devSelector.makeDevice(0);
    auto mesh = constructDummyMesh(0.0);
    auto result = runPropagationKernel(mesh, device, exec);

    REQUIRE(result.size() == propagationBatchSize * propagationResultWidth);
    for(unsigned id = 0u; id < propagationBatchSize; ++id)
    {
        auto const offset = id * propagationResultWidth;
        CAPTURE(id);
        REQUIRE(std::isfinite(result[offset]));
        CHECK(result[offset] == Catch::Approx(1.0));
        CHECK(result[offset + 1u] == Catch::Approx(0.0));
        CHECK(result[offset + 2u] == Catch::Approx(0.0));
    }
}

TEMPLATE_LIST_TEST_CASE("hase::kernels::propagateRayWithReflection responds to surface reflectivity", "", TestApis)
{
    auto cfg = TestType::makeDict();
    auto deviceSpec = cfg[alpaka::object::deviceSpec];
    auto exec = cfg[alpaka::object::exec];

    auto devSelector = alpaka::onHost::makeDeviceSelector(deviceSpec);
    if(!devSelector.isAvailable())
    {
        std::cout << "No device available for " << deviceSpec.getName() << std::endl;
        return;
    }

    auto device = devSelector.makeDevice(0);
    auto opaqueMesh = constructDummyMesh(0.0, 0.0f);
    auto reflectiveMesh = constructDummyMesh(0.0, 1.0f);

    auto const opaqueGain = runReflectionKernel(opaqueMesh, device, exec);
    auto const reflectiveGain = runReflectionKernel(reflectiveMesh, device, exec);

    REQUIRE(opaqueGain.size() == propagationBatchSize);
    REQUIRE(reflectiveGain.size() == propagationBatchSize);
    for(unsigned id = 0u; id < propagationBatchSize; ++id)
    {
        CAPTURE(id);
        REQUIRE(std::isfinite(opaqueGain[id]));
        REQUIRE(std::isfinite(reflectiveGain[id]));
        CHECK(opaqueGain[id] == Catch::Approx(0.0));
        CHECK(reflectiveGain[id] > opaqueGain[id]);
    }
}

TEMPLATE_LIST_TEST_CASE(
    "importance sampling propagation and distribution produce a usable ray distribution",
    "",
    TestApis)
{
    auto cfg = TestType::makeDict();
    auto deviceSpec = cfg[alpaka::object::deviceSpec];
    auto exec = cfg[alpaka::object::exec];

    auto devSelector = alpaka::onHost::makeDeviceSelector(deviceSpec);
    if(!devSelector.isAvailable())
    {
        std::cout << "No device available for " << deviceSpec.getName() << std::endl;
        return;
    }

    auto device = devSelector.makeDevice(0);
    auto hostMesh = constructDummyMesh(2.0);
    auto deviceMesh = hostMesh.toDevice(device);
    hase::alpakaUtils::DevBundle devBundle{device, exec};

    auto queue = device.makeQueue();
    constexpr unsigned reflectionSlices = 1u;
    constexpr unsigned raysPerSample = 16u;
    auto preImportance = alpaka::onHost::alloc<double>(device, reflectionSlices * deviceMesh.numberOfPrisms);
    auto importance = alpaka::onHost::alloc<double>(device, reflectionSlices * deviceMesh.numberOfPrisms);
    auto raysPerPrism = alpaka::onHost::alloc<unsigned>(device, reflectionSlices * deviceMesh.numberOfPrisms);
    unsigned droppedRays = 0u;
    hase::core::InfiniteRaySnapshot infiniteRaySnapshot{};
    auto droppedRaysView = alpaka::makeView(alpaka::api::host, &droppedRays, alpaka::Vec{1u});
    auto infiniteRaySnapshotView = alpaka::makeView(alpaka::api::host, &infiniteRaySnapshot, alpaka::Vec{1u});
    auto deviceDroppedRays = hase::alpakaUtils::toDevice(queue, droppedRaysView);
    auto deviceInfiniteRaySnapshots = hase::alpakaUtils::toDevice(queue, infiniteRaySnapshotView);

    hase::kernels::importanceSamplingPropagation(
        devBundle,
        0u,
        reflectionSlices,
        deviceMesh.toView(),
        0.0,
        0.0,
        preImportance,
        deviceDroppedRays,
        deviceInfiniteRaySnapshots);

    auto hostPreImportance = alpaka::onHost::allocHostLike(preImportance);
    alpaka::onHost::memcpy(queue, hostPreImportance, preImportance);
    alpaka::onHost::wait(queue);
    double const sumPhi = alpaka::onHost::data(hostPreImportance)[0];

    unsigned rngStride = 0u;
    unsigned const distributedRays = hase::kernels::importanceSamplingDistribution(
        devBundle,
        reflectionSlices,
        deviceMesh.toView(),
        raysPerSample,
        preImportance,
        importance,
        raysPerPrism,
        sumPhi,
        rngStride);

    auto hostImportance = alpaka::onHost::allocHostLike(importance);
    auto hostRaysPerPrism = alpaka::onHost::allocHostLike(raysPerPrism);
    alpaka::onHost::memcpy(queue, hostImportance, importance);
    alpaka::onHost::memcpy(queue, hostRaysPerPrism, raysPerPrism);
    alpaka::onHost::wait(queue);

    REQUIRE(std::isfinite(sumPhi));
    REQUIRE(sumPhi > 0.0);
    CHECK(distributedRays == raysPerSample);
    CHECK(alpaka::onHost::data(hostRaysPerPrism)[0] == raysPerSample);
    CHECK(std::isfinite(alpaka::onHost::data(hostImportance)[0]));
    CHECK(alpaka::onHost::data(hostImportance)[0] > 0.0);
}

TEMPLATE_LIST_TEST_CASE("importance sampling distribution skips zero total pre-importance", "", TestApis)
{
    auto cfg = TestType::makeDict();
    auto deviceSpec = cfg[alpaka::object::deviceSpec];
    auto exec = cfg[alpaka::object::exec];

    auto devSelector = alpaka::onHost::makeDeviceSelector(deviceSpec);
    if(!devSelector.isAvailable())
    {
        std::cout << "No device available for " << deviceSpec.getName() << std::endl;
        return;
    }

    auto device = devSelector.makeDevice(0);
    auto hostMesh = constructDummyMesh(0.0);
    auto deviceMesh = hostMesh.toDevice(device);
    hase::alpakaUtils::DevBundle devBundle{device, exec};

    auto queue = device.makeQueue();
    constexpr unsigned reflectionSlices = 1u;
    constexpr unsigned raysPerSample = 16u;
    auto preImportance = alpaka::onHost::alloc<double>(device, reflectionSlices * deviceMesh.numberOfPrisms);
    auto importance = alpaka::onHost::alloc<double>(device, reflectionSlices * deviceMesh.numberOfPrisms);
    auto raysPerPrism = alpaka::onHost::alloc<unsigned>(device, reflectionSlices * deviceMesh.numberOfPrisms);

    alpaka::onHost::fill(queue, preImportance, 0.0);
    alpaka::onHost::fill(queue, importance, -1.0);
    alpaka::onHost::fill(queue, raysPerPrism, 99u);
    alpaka::onHost::wait(queue);

    unsigned rngStride = 0u;
    unsigned const distributedRays = hase::kernels::importanceSamplingDistribution(
        devBundle,
        reflectionSlices,
        deviceMesh.toView(),
        raysPerSample,
        preImportance,
        importance,
        raysPerPrism,
        0.0,
        rngStride);

    auto hostImportance = alpaka::onHost::allocHostLike(importance);
    auto hostRaysPerPrism = alpaka::onHost::allocHostLike(raysPerPrism);
    alpaka::onHost::memcpy(queue, hostImportance, importance);
    alpaka::onHost::memcpy(queue, hostRaysPerPrism, raysPerPrism);
    alpaka::onHost::wait(queue);

    CHECK(distributedRays == 0u);
    CHECK(rngStride == 0u);
    CHECK(alpaka::onHost::data(hostRaysPerPrism)[0] == 0u);
    CHECK(alpaka::onHost::data(hostImportance)[0] == Catch::Approx(0.0));
}
