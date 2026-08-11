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

#include <array>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hase::core
{
    struct ForwardPhiAseRawResult
    {
        std::vector<double> vertexBatchScoreSum;
        std::vector<unsigned> rseBatchRayCounts;
        std::vector<unsigned> totalRays;
        std::vector<unsigned> droppedRays;
        unsigned rayCount = 0u;
        SrmStatus srmStatus = SrmStatus::DISABLED;
        unsigned srmPasses = 0u;
        double srmRemainingFraction = 0.0;
        unsigned srmMaxIterations = 0u;
        unsigned srmDivergenceStreak = 0u;
    };

    template<alpaka::onHost::concepts::Device T_Device>
    class ForwardSrmWorkspace
    {
        using T_DoubleBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<double>(std::declval<T_Device&>(), unsigned{1}));
        using T_UnsignedBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<unsigned>(std::declval<T_Device&>(), unsigned{1}));
        using T_CharBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<char>(std::declval<T_Device&>(), unsigned{1}));

    public:
        ForwardSrmWorkspace(T_Device& device, unsigned faceCount, unsigned reservoirSize, unsigned maxRayCount)
            : countsA(alpaka::onHost::alloc<unsigned>(device, faceCount))
            , countsB(alpaka::onHost::alloc<unsigned>(device, faceCount))
            , dirXA(alpaka::onHost::alloc<double>(device, faceCount * reservoirSize))
            , dirXB(alpaka::onHost::alloc<double>(device, faceCount * reservoirSize))
            , dirYA(alpaka::onHost::alloc<double>(device, faceCount * reservoirSize))
            , dirYB(alpaka::onHost::alloc<double>(device, faceCount * reservoirSize))
            , dirZA(alpaka::onHost::alloc<double>(device, faceCount * reservoirSize))
            , dirZB(alpaka::onHost::alloc<double>(device, faceCount * reservoirSize))
            , weightsA(alpaka::onHost::alloc<double>(device, faceCount * reservoirSize))
            , weightsB(alpaka::onHost::alloc<double>(device, faceCount * reservoirSize))
            , sigmaIndicesA(alpaka::onHost::alloc<unsigned>(device, faceCount * reservoirSize))
            , sigmaIndicesB(alpaka::onHost::alloc<unsigned>(device, faceCount * reservoirSize))
            , faceWeightsA(alpaka::onHost::alloc<double>(device, faceCount))
            , faceWeightsB(alpaka::onHost::alloc<double>(device, faceCount))
            , samplingCdf(alpaka::onHost::alloc<double>(device, faceCount))
            , samplingTotalWeight(alpaka::onHost::alloc<double>(device, 1u))
            , systematicOffset(alpaka::onHost::alloc<double>(device, 1u))
            , stratifiedRayCounts(alpaka::onHost::alloc<unsigned>(device, faceCount))
            , stratifiedRayOffsets(alpaka::onHost::alloc<unsigned>(device, faceCount))
            , stratifiedRayFaces(alpaka::onHost::alloc<unsigned>(device, maxRayCount))
            , samplingCdfScanBuffer(
                  alpaka::onHost::alloc<char>(
                      device,
                      alpaka::onHost::getScanBufferSize<double>(alpaka::Vec{faceCount})))
            , stratifiedCountScanBuffer(
                  alpaka::onHost::alloc<char>(
                      device,
                      alpaka::onHost::getScanBufferSize<unsigned>(alpaka::Vec{faceCount})))
            , faceCount(faceCount)
            , reservoirSize(reservoirSize)
            , maxRayCount(maxRayCount)
        {
        }

        T_UnsignedBuffer countsA;
        T_UnsignedBuffer countsB;
        T_DoubleBuffer dirXA;
        T_DoubleBuffer dirXB;
        T_DoubleBuffer dirYA;
        T_DoubleBuffer dirYB;
        T_DoubleBuffer dirZA;
        T_DoubleBuffer dirZB;
        T_DoubleBuffer weightsA;
        T_DoubleBuffer weightsB;
        T_UnsignedBuffer sigmaIndicesA;
        T_UnsignedBuffer sigmaIndicesB;
        T_DoubleBuffer faceWeightsA;
        T_DoubleBuffer faceWeightsB;
        T_DoubleBuffer samplingCdf;
        T_DoubleBuffer samplingTotalWeight;
        T_DoubleBuffer systematicOffset;
        T_UnsignedBuffer stratifiedRayCounts;
        T_UnsignedBuffer stratifiedRayOffsets;
        T_UnsignedBuffer stratifiedRayFaces;
        T_CharBuffer samplingCdfScanBuffer;
        T_CharBuffer stratifiedCountScanBuffer;
        unsigned faceCount;
        unsigned reservoirSize;
        unsigned maxRayCount;
    };

    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    void runForwardSrm(
        alpakaUtils::DevBundle<T_Device, T_Exec>& devBundle,
        concepts::Queue auto const& queue,
        DeviceMeshView const mesh,
        ExperimentParameters const& experiment,
        ForwardPhiAseRawResult& result,
        unsigned const rayCount,
        unsigned const rseBatch,
        double const betaVolumeTotal,
        alpaka::concepts::IBuffer auto& vertexBatchScoreSum,
        alpaka::concepts::IBuffer auto& volumeRayVisits,
        alpaka::concepts::IBuffer auto& droppedRays,
        alpaka::concepts::IBuffer auto const& sigmaA,
        alpaka::concepts::IBuffer auto const& sigmaE,
        unsigned const lambdaResolution,
        unsigned const threadLocalStridingRNG,
        SrmControls const srmControls,
        ForwardSrmWorkspace<T_Device>& workspace)
    {
        auto accumulationSpans = hase::kernels::forward::ForwardAccumulationSpans{
            vertexBatchScoreSum.getMdSpan(),
            volumeRayVisits.getMdSpan(),
            droppedRays.getMdSpan()};
        auto spectrumSpans = hase::kernels::forward::ForwardSpectrumSpans{sigmaA, sigmaE, lambdaResolution};
        unsigned const faceCount = mesh.numberOfCells * mesh.numberOfFacesPerCell;
        if(workspace.faceCount != faceCount || workspace.reservoirSize != experiment.surfaceReservoirSize
           || rayCount > workspace.maxRayCount)
            throw std::runtime_error("persistent forward SRM workspace does not match this launch");
        auto& countsA = workspace.countsA;
        auto& countsB = workspace.countsB;
        auto& dirXA = workspace.dirXA;
        auto& dirXB = workspace.dirXB;
        auto& dirYA = workspace.dirYA;
        auto& dirYB = workspace.dirYB;
        auto& dirZA = workspace.dirZA;
        auto& dirZB = workspace.dirZB;
        auto& weightsA = workspace.weightsA;
        auto& weightsB = workspace.weightsB;
        auto& sigmaIndicesA = workspace.sigmaIndicesA;
        auto& sigmaIndicesB = workspace.sigmaIndicesB;
        auto& faceWeightsA = workspace.faceWeightsA;
        auto& faceWeightsB = workspace.faceWeightsB;
        auto& samplingCdf = workspace.samplingCdf;
        auto& samplingTotalWeight = workspace.samplingTotalWeight;
        auto& systematicOffset = workspace.systematicOffset;
        auto& stratifiedRayCounts = workspace.stratifiedRayCounts;
        auto& stratifiedRayOffsets = workspace.stratifiedRayOffsets;
        auto& stratifiedRayFaces = workspace.stratifiedRayFaces;
        auto& samplingCdfScanBuffer = workspace.samplingCdfScanBuffer;
        auto& stratifiedCountScanBuffer = workspace.stratifiedCountScanBuffer;
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

        constexpr unsigned rayFrameExtent = 128u;
        auto const rayFrameSpec = alpaka::onHost::FrameSpec{
            alpaka::Vec{(rayCount + rayFrameExtent - 1u) / rayFrameExtent},
            alpaka::Vec{rayFrameExtent},
            devBundle.executor};
        queue.enqueue(
            rayFrameSpec,
            alpaka::KernelBundle{
                hase::kernels::forward::AccumulateForwardPhiAseReservoir{},
                mesh,
                rayCount,
                rseBatch,
                betaVolumeTotal,
                accumulationSpans,
                reservoirSpansA,
                spectrumSpans,
                threadLocalStridingRNG});
        alpaka::onHost::wait(queue);

        bool inputA = true;
        auto const faceFrameSpec
            = hase::alpakaUtils::getFrameSpec<uint32_t>(devBundle.device, devBundle.executor, alpaka::Vec{faceCount});
        auto const scalarFrameSpec
            = hase::alpakaUtils::getFrameSpec<uint32_t>(devBundle.device, devBundle.executor, alpaka::Vec{1u});
        auto samplingTotalWeightHost = alpaka::onHost::allocHostLike(samplingTotalWeight);
        auto updateSamplingCdf = [&](auto const& reservoir, unsigned const pass)
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
                        threadLocalStridingRNG,
                        pass});
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
        double const initialWeight = updateSamplingCdf(reservoirSpansA, 0u);
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
                    rayFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::AccumulateReflectedForwardPhiAse{},
                        mesh,
                        rayCount,
                        rseBatch,
                        sourceWeight,
                        accumulationSpans,
                        reservoirSpansA,
                        samplingCdfSpans,
                        reservoirSpansB,
                        spectrumSpans,
                        threadLocalStridingRNG,
                        pass});
            }
            else
            {
                clearReservoir(countsA, faceWeightsA);
                alpaka::onHost::wait(queue);
                queue.enqueue(
                    rayFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::AccumulateReflectedForwardPhiAse{},
                        mesh,
                        rayCount,
                        rseBatch,
                        sourceWeight,
                        accumulationSpans,
                        reservoirSpansB,
                        samplingCdfSpans,
                        reservoirSpansA,
                        spectrumSpans,
                        threadLocalStridingRNG,
                        pass});
            }
            alpaka::onHost::wait(queue);
            inputA = !inputA;

            double const currentWeight = updateSamplingCdf(inputA ? reservoirSpansA : reservoirSpansB, pass);
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
