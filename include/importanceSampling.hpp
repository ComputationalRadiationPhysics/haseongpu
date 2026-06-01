/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once

#include <detail/importanceSampling.hpp>
#include <mesh.hpp>

// struct ValidateEmpty
// {
//     ALPAKA_FN_ACC constexpr void operator() (auto const & acc,alpaka::concepts::IMdSpan auto importance) const
//     {
//         for(auto
//         id:alpaka::onAcc::makeIdxMap(acc,alpaka::onAcc::worker::threadsInGrid,alpaka::IdxRange{importance}))
//         {
//
//         }
//     }
// };
/**
 * @brief Calculates preImportance which needs only to be done once.
 *
 * @param sample_i         Index of sample point for which importance sampling should be done.
 * @param reflectionSlices 1 + (2 * maxReflections) - Coded information about how many
 *                         reflections a ray should do and on which surface to start.
 * @param dMesh            All information about triangles, points, contants.
 *                         Is located on device memory. See mesh.h for details.
 * @param sigmaA           Absorption value of the ray.
 * @param sigmaE           Emission value of the ray.
 * @param preImportance    Importance based on gain from test rays (will be returned as pointer).
 * @param blockDim         Number of threads per block.
 * @param gridDim          Number of blocks per grid.
 *
 */
void importanceSamplingPropagation(
    auto& devBundle, // makes this function a template -> hence the function definition was moved into this header file
    unsigned sample_i,
    unsigned const reflectionSlices,
    DeviceMeshView const& deviceMesh,
    double const sigmaA,
    double const sigmaE,
    alpaka::concepts::IMdSpan auto preImportance)
{
    auto queue = devBundle.device.makeQueue();
    auto validateMeshFrameSpec
        = alpaka::onHost::getFrameSpec<Vec1D::type>(queue.getDevice(), Vec1D{deviceMesh.numberOfSamples});
    queue.enqueue(
        devBundle.executor,
        validateMeshFrameSpec,
        alpaka::KernelBundle{ValidateMeshKernel{}, deviceMesh, sample_i});
    alpaka::onHost::wait(queue);
    auto propagateFrameSpec = alpaka::onHost::FrameSpec{
        Vec3D{reflectionSlices, 1, validateMeshFrameSpec.getNumFrames().x()},
        Vec3D{1, 1, validateMeshFrameSpec.getFrameExtents().x()}};
    queue.enqueue(
        devBundle.executor,
        propagateFrameSpec,
        alpaka::KernelBundle{
            PropagateFromTriangleCenter{},
            deviceMesh,
            preImportance,
            reflectionSlices,
            deviceMesh.numberOfPrisms,
            sample_i,
            sigmaA,
            sigmaE});
    alpaka::onHost::wait(queue);
}

/**
 * @brief Calculates importance and ray distribution on prisms.
 *        Based on preImportance and triangle surfaces.
 *
 * @param reflectionSlices   1 + (2 * maxReflections) - Coded information about how many
 *                           reflections a ray should do and on which surface to start.
 * @param dMesh              All information about triangles, points, contants.
 *                           Is located on device memory. See mesh.h for details.
 * @param raysPerSample      Number of rays that should be distributed on gain medium.
 * @param preImportance      Values calculated from importanceSamplingPropagation.
 * @param importance         Final importance values for this number of rays per sample.
 * @param raysPerPrism       Information about how many rays will start in each prism.
 * @param hSumPhi            Global memory where all threads can sum their phiAse on.
 * @param distributeRandomly In case that not all rays can be distributed by importance
 *                           sampling, the rest will be distributed randomly.
 * @param blockDim           Number of threads per block.
 * @param gridDim            Number of blocks per grid.
 *
 */
unsigned importanceSamplingDistribution(
    auto& devBundle,
    unsigned reflectionSlices,
    DeviceMeshView deviceMesh,
    unsigned raysPerSample,
    alpaka::concepts::IMdSpan auto preImportance,
    alpaka::concepts::IMdSpan auto importance,
    alpaka::concepts::IMdSpan auto raysPerPrism,
    double hSumPhi,
    bool const distributeRandomly,
    unsigned& threadLocalStridingRNG)
{
    auto queue = devBundle.device.makeQueue();
    unsigned hRaysDump = 0;
    auto hRaysDumpView = alpaka::makeView(alpaka::api::host, &hRaysDump, alpaka::Vec{1u});
    auto dRaysDump = hase::alpakaUtils::toDevice(queue, hRaysDumpView);
    using Vec3D = alpaka::Vec<uint32_t, 3>;
    auto meshFrameSpec
        = alpaka::onHost::getFrameSpec<uint32_t>(devBundle.device, Vec3D{reflectionSlices, 1, raysPerSample});
    queue.enqueue(
        devBundle.executor,
        meshFrameSpec,
        alpaka::KernelBundle{
            DistributeRaysByImportance{},
            deviceMesh,
            raysPerPrism,
            preImportance,
            hSumPhi, // we can copy trivial types directly as part of the kernel arguments
            reflectionSlices,
            raysPerSample});
    using ValueType = alpaka::trait::GetValueType_t<ALPAKA_TYPEOF(raysPerPrism)>;
    alpaka::onHost::reduce(queue, devBundle.executor, ValueType{0}, dRaysDump, std::plus{}, raysPerPrism);
    // copy number of distributed rays after distribution
    alpaka::onHost::memcpy(queue, hRaysDumpView, dRaysDump);

    alpaka::onHost::wait(queue);
    // assert that the overall distibuted number of rays doesnt exceed our maximum
    assert(hRaysDump <= raysPerSample);
    // Distribute remaining rays randomly if wanted
    if(distributeRandomly)
    {
        unsigned raysLeft = raysPerSample - hRaysDump;
        auto distributeFrameSpec = alpaka::onHost::getFrameSpec<uint32_t>(devBundle.device, Vec1D{raysLeft});
        // since frameSpec is always > threadSpec we can assume that this always traverses over the whole grid/workers
        auto const threadLocalStridingIndex = threadLocalStridingRNG;
        queue.enqueue(
            devBundle.executor,
            distributeFrameSpec,
            alpaka::KernelBundle{
                DistributeRemaingRaysRandomly{},
                deviceMesh,
                raysPerPrism,
                raysLeft,
                dRaysDump,
                threadLocalStridingIndex});
        threadLocalStridingRNG
            += (distributeFrameSpec.getNumFrames().product() * distributeFrameSpec.getFrameExtents().product());
        alpaka::onHost::memcpy(queue, hRaysDumpView, dRaysDump);
#ifndef NDEBUG
        alpaka::onHost::wait(queue);
#endif
        assert(hRaysDump <= raysPerSample);
    }
    queue.enqueue(
        devBundle.executor,
        meshFrameSpec,
        alpaka::KernelBundle{
            RecalculateImportance{},
            deviceMesh,
            raysPerPrism,
            dRaysDump,
            importance,
            reflectionSlices});
    alpaka::onHost::wait(queue);
    return hRaysDump;
}
