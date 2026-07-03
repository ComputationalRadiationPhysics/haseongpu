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
#include <kernels/reflection.hpp>

#include <cmath>
#include <iostream>
#include <type_traits>
#include <vector>

using TestApis
    = std::decay_t<decltype(alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors))>;

constexpr unsigned propagationBatchSize = 32u;
constexpr unsigned propagationResultWidth = 2u;

hase::core::HostMesh constructDummyMesh(double betaVolume)
{
    hase::core::HostMesh mesh;
    mesh.numberOfCells = 2u;
    mesh.numberOfPrisms = mesh.numberOfCells;
    mesh.numberOfMeshPoints = 5u;
    mesh.numberOfPoints = 1u;
    mesh.numberOfSamples = 1u;
    mesh.numberOfTriangles = mesh.numberOfCells;
    mesh.numberOfLevels = 1u;
    mesh.numberOfFacesPerCell = hase::core::tet4FaceCount;
    mesh.numberOfCellVertices = hase::core::tet4VertexCount;
    mesh.thickness = 1.0f;
    mesh.points = {
        1.0 / 3.0,
        0.0,
        1.0,
        0.0,
        1.0 / 3.0,
        1.0 / 3.0,
        0.0,
        0.0,
        1.0,
        1.0 / 3.0,
        0.0,
        0.5,
        0.5,
        0.5,
        1.0,
    };
    mesh.samplePoints = {1.0 / 3.0, 1.0 / 3.0, 0.25};
    mesh.cellCenters = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.375, 0.625};
    mesh.cellPointIndices = {0u, 1u, 2u, 3u, 4u, 1u, 2u, 3u};
    mesh.cellTypes = {hase::core::vtkTetraCellType, hase::core::vtkTetraCellType};
    mesh.cellFaces = {
        0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3, 4, 1, 2, 4, 2, 3, 4, 3, 1, 1, 3, 2,
    };
    mesh.cellNeighborCells = {-1, -1, -1, 1, -1, -1, -1, 0};
    mesh.cellNeighborLocalFaces = {-1, -1, -1, 3, -1, -1, -1, 3};
    mesh.cellFaceBoundaries = {-1, -1, -1, 0, -1, -1, -1, 0};
    mesh.cellVolumes = {1.0f / 6.0f, 1.0f / 6.0f};
    mesh.betaVolume = {betaVolume, betaVolume};
    mesh.betaCells = {betaVolume};
    mesh.claddingCellTypes = {0u, 0u};
    mesh.refractiveIndices = {1.5f, 1.0f, 1.5f, 1.0f};
    mesh.reflectivities = {1.0f, 1.0f};
    mesh.nTot = 1.0f;
    mesh.crystalTFluo = 1.0f;
    mesh.claddingNumber = 99u;
    mesh.claddingAbsorption = 0.0;
    mesh.calcCellVolumePrefix();
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
                alpaka::IdxRange{propagationBatchSize}))
        {
            unsigned level = 0u;
            double const zStart = 0.05 + 0.01 * static_cast<double>(id);
            hase::core::Ray ray{{1.0 / 3.0, 1.0 / 3.0, zStart}, {0.0, 0.0, 0.4}, 0.4f};

            auto const offset = id * propagationResultWidth;
            result[offset] = hase::kernels::propagateRay(ray, &level, mesh, 0.0, 0.0);
            result[offset + 1u] = static_cast<double>(level);
        }
    }
};

struct ReflectionPlaneKernel
{
    ALPAKA_FN_ACC void operator()(auto const& acc, hase::core::DeviceMeshView const mesh, auto result) const
    {
        for(auto [id] : alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{1u}))
        {
            hase::core::Point reflectionPoint{0.0, 0.0, 0.0};
            double reflectionAngle = 0.0;
            result[id * 2u] = static_cast<double>(hase::kernels::calcNextReflection(
                {1.0 / 3.0, 1.0 / 3.0, 0.25},
                {1.0 / 3.0, 1.0 / 3.0, 0.25},
                1u,
                hase::core::TOP_REFLECTION,
                &reflectionPoint,
                &reflectionAngle,
                mesh));
            result[id * 2u + 1u] = reflectionPoint.z;
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
                alpaka::IdxRange{propagationBatchSize}))
        {
            double const z = 0.1 + 0.01 * static_cast<double>(id);
            hase::core::Point startPoint{1.0 / 3.0, 1.0 / 3.0, z};
            hase::core::Point endPoint{1.0 / 3.0, 1.0 / 3.0, z + 0.2};

            result[id] = hase::kernels::propagateRayWithReflection(
                startPoint,
                endPoint,
                0u,
                hase::core::TOP_REFLECTION,
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
    auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(device, executor, alpaka::Vec{propagationBatchSize});
    queue.enqueue(frameSpec, alpaka::KernelBundle{PropagationKernel{}, deviceMesh.toView(), result});
    auto hostResult = alpaka::onHost::allocHostLike(result);
    alpaka::onHost::memcpy(queue, hostResult, result);
    alpaka::onHost::wait(queue);

    auto* data = alpaka::onHost::data(hostResult);
    return {data, data + propagationBatchSize * propagationResultWidth};
}

template<typename T_Device, typename T_Executor>
std::vector<double> runReflectionPlaneKernel(
    hase::core::HostMesh& hostMesh,
    T_Device& device,
    T_Executor const& executor)
{
    auto queue = device.makeQueue();
    auto deviceMesh = hostMesh.toDevice(device);
    auto result = alpaka::onHost::alloc<double>(device, 2u);
    auto frameSpec = alpaka::onHost::getFrameSpec(device, executor, 1u);
    queue.enqueue(frameSpec, alpaka::KernelBundle{ReflectionPlaneKernel{}, deviceMesh.toView(), result});
    auto hostResult = alpaka::onHost::allocHostLike(result);
    alpaka::onHost::memcpy(queue, hostResult, result);
    alpaka::onHost::wait(queue);
    auto* data = alpaka::onHost::data(hostResult);
    return {data, data + 2u};
}

template<typename T_Device, typename T_Executor>
std::vector<double> runReflectionKernel(hase::core::HostMesh& hostMesh, T_Device& device, T_Executor const& executor)
{
    auto queue = device.makeQueue();
    auto deviceMesh = hostMesh.toDevice(device);
    auto result = alpaka::onHost::alloc<double>(device, propagationBatchSize);
    auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(device, executor, alpaka::Vec{propagationBatchSize});
    queue.enqueue(frameSpec, alpaka::KernelBundle{ReflectionKernel{}, deviceMesh.toView(), result});
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
        double const zStart = 0.05 + 0.01 * static_cast<double>(id);
        double const expectedCell = zStart + 0.4 >= 0.5 ? 1.0 : 0.0;
        CHECK(result[offset + 1u] == Catch::Approx(expectedCell));
    }
}

TEMPLATE_LIST_TEST_CASE("hase::kernels::propagateRayWithReflection preserves no-reflection gain", "", TestApis)
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
    auto const gain = runReflectionKernel(mesh, device, exec);

    REQUIRE(gain.size() == propagationBatchSize);
    for(unsigned id = 0u; id < propagationBatchSize; ++id)
    {
        CAPTURE(id);
        REQUIRE(std::isfinite(gain[id]));
        CHECK(gain[id] == Catch::Approx(25.0));
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
        queue,
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
    auto const* preImportanceData = alpaka::onHost::data(hostPreImportance);
    double sumPhi = 0.0;
    for(unsigned i = 0u; i < reflectionSlices * deviceMesh.numberOfPrisms; ++i)
    {
        sumPhi += preImportanceData[i];
    }

    unsigned rngStride = 0u;
    unsigned const distributedRays = hase::kernels::importanceSamplingDistribution(
        devBundle,
        queue,
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
    auto const* raysPerPrismData = alpaka::onHost::data(hostRaysPerPrism);
    auto const* importanceData = alpaka::onHost::data(hostImportance);
    unsigned distributedRayCount = 0u;
    unsigned cellsWithRays = 0u;
    for(unsigned i = 0u; i < reflectionSlices * deviceMesh.numberOfPrisms; ++i)
    {
        distributedRayCount += raysPerPrismData[i];
        cellsWithRays += raysPerPrismData[i] > 0u ? 1u : 0u;
        CHECK(std::isfinite(importanceData[i]));
        CHECK(importanceData[i] > 0.0);
    }
    CHECK(distributedRayCount == raysPerSample);
    CHECK(cellsWithRays > 0u);
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
        queue,
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
