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
                         // file
        unsigned sample_i,
        unsigned const reflectionSlices,
        hase::core::DeviceMeshView const& deviceMesh,
        double const sigmaA,
        double const sigmaE,
        alpaka::concepts::IMdSpan auto preImportance)
    {
        auto queue = devBundle.device.makeQueue();
        auto validateMeshFrameSpec
            = alpaka::onHost::getFrameSpec<Vec1D::type>(queue.getDevice(), Vec1D{deviceMesh.numberOfSamples});
        // wrapper for NDEBUG only be compiled in debug mode
        ALPAKA_ASSERT((
            [&]
            {
                queue.enqueue(
                    devBundle.executor,
                    validateMeshFrameSpec,
                    alpaka::KernelBundle{hase::kernels::ValidateMeshKernel{}, deviceMesh, sample_i});

                alpaka::onHost::wait(queue);

                return true;
            }()));
        alpaka::onHost::wait(queue);
        auto propagateFrameSpec = alpaka::onHost::FrameSpec{
            Vec3D{reflectionSlices, 1, validateMeshFrameSpec.getNumFrames().x()},
            Vec3D{1, 1, validateMeshFrameSpec.getFrameExtents().x()}};
        alpaka::onHost::concurrent<double>(
            queue,
            devBundle.executor,
            preImportance.getExtents(),
            hase::kernels::PropagateFromTriangleCenter{deviceMesh, sample_i, sigmaA, sigmaE},
            preImportance);
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
        hase::core::DeviceMeshView deviceMesh,
        unsigned raysPerSample,
        alpaka::concepts::IBuffer auto preImportance,
        alpaka::concepts::IBuffer auto importance,
        alpaka::concepts::IBuffer auto raysPerPrism,
        double hSumPhi,
        unsigned& threadLocalStridingRNG)
    {
        auto queue = devBundle.device.makeQueue();
        alpaka::onHost::fill(queue, raysPerPrism, 0u);
        unsigned hRaysDump = 0;
        auto hRaysDumpView = alpaka::makeView(alpaka::api::host, &hRaysDump, alpaka::Vec{1u});
        auto dRaysDump = hase::alpakaUtils::toDevice(queue, hRaysDumpView);
        using Vec3D = alpaka::Vec<uint32_t, 3>;
        auto meshFrameSpec
            = alpaka::onHost::getFrameSpec<uint32_t>(devBundle.device, Vec3D{reflectionSlices, 1, raysPerSample});
        if(!std::isfinite(hSumPhi) || hSumPhi <= 0.0)
        {
            alpaka::onHost::fill(queue, importance, 0.0);
            alpaka::onHost::wait(queue);
            return 0u;
        }
        alpaka::onHost::inclusiveScan(queue, devBundle.executor, importance, preImportance);

        auto const drawFrameSpec = alpaka::onHost::getFrameSpec<uint32_t>(devBundle.device, Vec1D{raysPerSample});
        auto const drawThreadLocalStridingIndex = threadLocalStridingRNG;
        queue.enqueue(
            devBundle.executor,
            drawFrameSpec,
            alpaka::KernelBundle{
                hase::kernels::DrawRaysByImportance{},
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
            hase::kernels::RecalculateUnbiasedImportance{deviceMesh, hSumPhi},
            preImportance,
            importance);

        alpaka::onHost::wait(queue);
        return raysPerSample;
    }

} // namespace hase::kernels
