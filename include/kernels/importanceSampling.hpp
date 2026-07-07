/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 * Copyright 2026 Tim Hanel
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

#include <alpakaUtils/DevBundle.hpp>
#include <concepts/concepts.hpp>
#include <core/mesh.hpp>
#include <kernels/detail/importanceSampling.hpp>

#include <cmath>

namespace hase::kernels
{
    using hase::alpakaUtils::Vec1D;
    using hase::alpakaUtils::Vec3D;

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
        auto& devBundle, // makes this function a template -> hence the function definition was moved into this header
        concepts::Queue auto const&
            queue, // we need to hand in the queue to enable synchronization being orchestrated from the outside
        unsigned sample_i,
        unsigned const reflectionSlices,
        core::DeviceMeshView const& deviceMesh,
        double const sigmaA,
        double const sigmaE,
        alpaka::concepts::IMdSpan auto preImportance,
        alpaka::concepts::IMdSpan auto droppedRays,
        alpaka::concepts::IMdSpan auto infiniteRaySnapshots)
    {
        // wrapper for NDEBUG only be compiled in debug mode
        ALPAKA_ASSERT((
            [&]
            {
                auto validateMeshFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                    queue.getDevice(),
                    devBundle.executor,
                    alpaka::Vec{deviceMesh.numberOfSamples});
                queue.enqueue(
                    validateMeshFrameSpec,
                    alpaka::KernelBundle{hase::kernels::ValidateMeshKernel{}, deviceMesh, sample_i});

                alpaka::onHost::wait(queue);

                return true;
            }()));
        alpaka::concepts::IMdSpan auto preImportanceSpan = preImportance.getMdSpan();
        alpaka::concepts::IMdSpan auto droppedRaysSpan = droppedRays.getMdSpan();
        alpaka::concepts::IMdSpan auto infiniteRaySnapshotsSpan = infiniteRaySnapshots.getMdSpan();
        auto propagateFrameSpec = alpakaUtils::getFrameSpec<uint32_t>(
            queue.getDevice(),
            devBundle.executor,
            alpaka::Vec{static_cast<uint32_t>(reflectionSlices * deviceMesh.numberOfPrisms)});
        queue.enqueue(
            propagateFrameSpec,
            alpaka::KernelBundle{
                PropagateFromTriangleCenter{},
                deviceMesh,
                sample_i,
                sigmaA,
                sigmaE,
                preImportanceSpan,
                droppedRaysSpan,
                infiniteRaySnapshotsSpan});
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
        concepts::Queue auto const&
            queue, // we need to hand in the queue to enable synchronization being orchestrated from the outside
        unsigned reflectionSlices,
        core::DeviceMeshView deviceMesh,
        unsigned raysPerSample,
        alpaka::concepts::IBuffer auto preImportance,
        alpaka::concepts::IBuffer auto importance,
        alpaka::concepts::IBuffer auto raysPerPrism,
        double hSumPhi,
        unsigned& threadLocalStridingRNG)
    {
        alpaka::onHost::memset(queue, raysPerPrism, 0);
        unsigned hRaysDump = 0;
        auto hRaysDumpView = alpaka::makeView(alpaka::api::host, &hRaysDump, alpaka::Vec{1u});
        auto dRaysDump = alpakaUtils::toDevice(queue, hRaysDumpView);
        if(!std::isfinite(hSumPhi) || hSumPhi <= 0.0)
        {
            alpaka::onHost::memset(queue, importance, 0);
            alpaka::onHost::wait(queue);
            return 0u;
        }
        alpaka::onHost::inclusiveScan(queue, devBundle.executor, importance, preImportance);

        auto const drawFrameSpec
            = alpakaUtils::getFrameSpec<uint32_t>(devBundle.device, devBundle.executor, alpaka::Vec{raysPerSample});
        auto const drawThreadLocalStridingIndex = threadLocalStridingRNG;
        queue.enqueue(
            drawFrameSpec,
            alpaka::KernelBundle{
                DrawRaysByImportance{},
                importance,
                raysPerPrism,
                hSumPhi,
                deviceMesh.numberOfPrisms * reflectionSlices,
                raysPerSample,
                drawThreadLocalStridingIndex});
        threadLocalStridingRNG += (drawFrameSpec.getNumFrames().product() * drawFrameSpec.getFrameExtents().product());

        alpaka::onHost::concurrent<double>(
            queue,
            devBundle.executor,
            preImportance.getExtents(),
            RecalculateUnbiasedImportance{deviceMesh, hSumPhi},
            preImportance,
            importance);
        return raysPerSample;
    }

} // namespace hase::kernels
