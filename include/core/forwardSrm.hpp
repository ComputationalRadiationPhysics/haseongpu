/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpakaUtils/DevBundle.hpp>
#include <concepts/concepts.hpp>
#include <core/mesh.hpp>
#include <core/srm.hpp>
#include <core/types.hpp>
#include <kernels/forward/accumulation.hpp>

#include <limits>
#include <vector>

namespace hase::core
{
    struct ForwardPhiAseRawResult
    {
        std::vector<double> scoreSum;
        std::vector<double> scoreSquareSum;
        std::vector<unsigned> totalRays;
        std::vector<unsigned> droppedRays;
        unsigned rayCount = 0u;
        SrmStatus srmStatus = SrmStatus::DISABLED;
        unsigned srmPasses = 0u;
        double srmRemainingFraction = 0.0;
        unsigned srmMaxIterations = 0u;
        unsigned srmDivergenceStreak = 0u;
    };

    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    void runForwardSrm(
        alpakaUtils::DevBundle<T_Device, T_Exec>& devBundle,
        concepts::Queue auto const& queue,
        DeviceMeshView const mesh,
        ExperimentParameters const& experiment,
        ForwardPhiAseRawResult& result,
        unsigned const rayCount,
        double const betaVolumeTotal,
        alpaka::concepts::IBuffer auto& phi,
        alpaka::concepts::IBuffer auto& phiSquare,
        alpaka::concepts::IBuffer auto& volumeRayVisits,
        alpaka::concepts::IBuffer auto& droppedRays,
        alpaka::concepts::IBuffer auto const& sigmaA,
        alpaka::concepts::IBuffer auto const& sigmaE,
        unsigned const lambdaResolution,
        unsigned const threadLocalStridingRNG,
        SrmControls const srmControls)
    {
        auto accumulationSpans = hase::kernels::forward::ForwardAccumulationSpans{
            phi,
            phiSquare,
            volumeRayVisits,
            droppedRays};
        auto spectrumSpans = hase::kernels::forward::ForwardSpectrumSpans{sigmaA, sigmaE, lambdaResolution};
        unsigned const faceCount = mesh.numberOfCells * mesh.numberOfFacesPerCell;
        unsigned const reservoirSlots = faceCount * experiment.surfaceReservoirSize;
        alpaka::concepts::IBuffer auto countsA = alpaka::onHost::alloc<unsigned>(devBundle.device, faceCount);
        alpaka::concepts::IBuffer auto countsB = alpaka::onHost::alloc<unsigned>(devBundle.device, faceCount);
        alpaka::concepts::IBuffer auto dirXA = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto dirXB = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto dirYA = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto dirYB = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto dirZA = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto dirZB = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto weightsA = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto weightsB = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto sigmaIndicesA = alpaka::onHost::alloc<unsigned>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto sigmaIndicesB = alpaka::onHost::alloc<unsigned>(devBundle.device, reservoirSlots);
        alpaka::concepts::IBuffer auto faceWeightsA = alpaka::onHost::alloc<double>(devBundle.device, faceCount);
        alpaka::concepts::IBuffer auto faceWeightsB = alpaka::onHost::alloc<double>(devBundle.device, faceCount);
        alpaka::concepts::IBuffer auto samplingCdf = alpaka::onHost::alloc<double>(devBundle.device, faceCount);
        alpaka::concepts::IBuffer auto samplingTotalWeight = alpaka::onHost::alloc<double>(devBundle.device, 1u);
        alpaka::concepts::IBuffer auto systematicOffset = alpaka::onHost::alloc<double>(devBundle.device, 1u);
        alpaka::concepts::IBuffer auto stratifiedRayCounts = alpaka::onHost::alloc<unsigned>(devBundle.device, faceCount);
        alpaka::concepts::IBuffer auto stratifiedRayOffsets = alpaka::onHost::alloc<unsigned>(devBundle.device, faceCount);
        alpaka::concepts::IBuffer auto stratifiedRayFaces = alpaka::onHost::alloc<unsigned>(devBundle.device, rayCount);
        alpaka::concepts::IBuffer auto samplingCdfScanBuffer = alpaka::onHost::alloc<char>(
            devBundle.device,
            alpaka::onHost::getScanBufferSize<double>(alpaka::Vec{faceCount}));
        alpaka::concepts::IBuffer auto stratifiedCountScanBuffer = alpaka::onHost::alloc<char>(
            devBundle.device,
            alpaka::onHost::getScanBufferSize<unsigned>(alpaka::Vec{faceCount}));
        auto reservoirSpansA = hase::kernels::forward::SurfaceReservoirSpans{
            countsA,
            dirXA,
            dirYA,
            dirZA,
            weightsA,
            sigmaIndicesA,
            faceWeightsA,
            experiment.surfaceReservoirSize};
        auto samplingCdfSpans = hase::kernels::forward::SurfaceReservoirSamplingCdfSpans{
            samplingCdf,
            samplingTotalWeight,
            stratifiedRayFaces,
            faceCount <= rayCount};
        auto reservoirSpansB = hase::kernels::forward::SurfaceReservoirSpans{
            countsB,
            dirXB,
            dirYB,
            dirZB,
            weightsB,
            sigmaIndicesB,
            faceWeightsB,
            experiment.surfaceReservoirSize};

        auto clearReservoir = [&](auto& counts, auto& faceWeights)
        {
            alpaka::onHost::fill(queue, counts, 0u, alpaka::Vec{faceCount});
            alpaka::onHost::fill(queue, faceWeights, 0.0, alpaka::Vec{faceCount});
        };

        clearReservoir(countsA, faceWeightsA);
        clearReservoir(countsB, faceWeightsB);
        alpaka::onHost::wait(queue);

        auto const directFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{static_cast<unsigned int>(rayCount)});
        queue.enqueue(
            directFrameSpec,
            alpaka::KernelBundle{
                hase::kernels::forward::AccumulateForwardPhiAseReservoir{},
                mesh,
                rayCount,
                betaVolumeTotal,
                accumulationSpans,
                reservoirSpansA,
                spectrumSpans,
                threadLocalStridingRNG});
        alpaka::onHost::wait(queue);

        auto const reflectedFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{static_cast<unsigned int>(rayCount)});
        bool inputA = true;
        auto const faceFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{faceCount});
        auto const scalarFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{1u});
        auto samplingTotalWeightHost = alpaka::onHost::allocHostLike(samplingTotalWeight);
        auto updateSamplingCdf = [&](auto const& reservoir, unsigned const seed)
        {
            alpaka::onHost::inclusiveScan(
                queue,
                devBundle.executor,
                samplingCdfScanBuffer,
                samplingCdf,
                reservoir.faceWeights);
            queue.enqueue(
                scalarFrameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::CaptureSurfaceReservoirSamplingTotalWeight{},
                    faceCount,
                    samplingCdfSpans});
            queue.enqueue(
                faceFrameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::NormalizeSurfaceReservoirSamplingCdf{},
                    faceCount,
                    samplingCdfSpans});
            if(samplingCdfSpans.useFaceStratification)
            {
                queue.enqueue(
                    scalarFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::GenerateSurfaceReservoirSystematicOffset{},
                        systematicOffset,
                        seed});
                queue.enqueue(
                    faceFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::AssignSurfaceReservoirStratifiedRayCounts{},
                        faceCount,
                        rayCount,
                        samplingCdfSpans,
                        systematicOffset,
                        stratifiedRayCounts});
                alpaka::onHost::exclusiveScan(
                    queue,
                    devBundle.executor,
                    stratifiedCountScanBuffer,
                    stratifiedRayOffsets,
                    stratifiedRayCounts);
                queue.enqueue(
                    faceFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::ScatterSurfaceReservoirStratifiedRayFaces{},
                        faceCount,
                        stratifiedRayCounts,
                        stratifiedRayOffsets,
                        stratifiedRayFaces});
            }
            alpaka::onHost::wait(queue);
            alpaka::onHost::memcpy(queue, samplingTotalWeightHost, samplingTotalWeight);
            alpaka::onHost::wait(queue);
            return alpaka::onHost::data(samplingTotalWeightHost)[0u];
        };
        double const initialWeight = updateSamplingCdf(reservoirSpansA, threadLocalStridingRNG);
        if(initialWeight == 0.0)
        {
            result.srmStatus = SrmStatus::CONVERGED;
            return;
        }

        double previousWeight = initialWeight;
        unsigned growCount = 0u;
        result.srmStatus = SrmStatus::MAX_ITERATIONS;
        result.srmRemainingFraction = 1.0;
        for(unsigned pass = 1u; pass <= srmControls.maxIterations; ++pass)
        {
            double const sourceWeight = previousWeight / static_cast<double>(rayCount);
            // Buffer switch for reservoir weights.
            if(inputA)
            {
                clearReservoir(countsB, faceWeightsB);
                alpaka::onHost::wait(queue);
                queue.enqueue(
                    reflectedFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::AccumulateReflectedForwardPhiAse{},
                        mesh,
                        rayCount,
                        sourceWeight,
                        accumulationSpans,
                        reservoirSpansA,
                        samplingCdfSpans,
                        reservoirSpansB,
                        spectrumSpans,
                        threadLocalStridingRNG + pass * rayCount});
            }
            else
            {
                clearReservoir(countsA, faceWeightsA);
                alpaka::onHost::wait(queue);
                queue.enqueue(
                    reflectedFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::AccumulateReflectedForwardPhiAse{},
                        mesh,
                        rayCount,
                        sourceWeight,
                        accumulationSpans,
                        reservoirSpansB,
                        samplingCdfSpans,
                        reservoirSpansA,
                        spectrumSpans,
                        threadLocalStridingRNG + pass * rayCount});
            }
            alpaka::onHost::wait(queue);
            inputA = !inputA;

            double const currentWeight = updateSamplingCdf(
                inputA ? reservoirSpansA : reservoirSpansB,
                threadLocalStridingRNG + pass * rayCount);
            result.srmPasses = pass;
            result.srmRemainingFraction = currentWeight / initialWeight;
            if(currentWeight > previousWeight)
            {
                ++growCount;
                if(growCount >= srmControls.divergenceStreak)
                {
                    result.srmStatus = SrmStatus::DIVERGED;
                    break;
                }
            }
            else
            {
                growCount = 0u;
                if(alpaka::math::abs(currentWeight - previousWeight) / alpaka::math::max(currentWeight, 1.0e-30)
                   < experiment.reflectionTolerance)
                {
                    result.srmStatus = SrmStatus::STABLE;
                    break;
                }
            }
            if(result.srmRemainingFraction < experiment.reflectionTolerance)
            {
                result.srmStatus = SrmStatus::CONVERGED;
                break;
            }
            previousWeight = currentWeight;
        }
    }
} // namespace hase::core
