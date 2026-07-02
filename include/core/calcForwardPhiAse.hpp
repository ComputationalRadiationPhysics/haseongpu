/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpakaUtils/DevBundle.hpp>
#include <core/logging.hpp>
#include <core/mesh.hpp>
#include <core/types.hpp>
#include <kernels/forward/accumulation.hpp>

#include <algorithm>
#include <cmath>
#include <ctime>
#include <limits>
#include <stdexcept>
#include <vector>

namespace hase::core
{
    struct BetaVolumeContribution
    {
        constexpr auto operator()(alpaka::concepts::Simd auto const& beta, alpaka::concepts::Simd auto const& volume)
            const
        {
            return beta * alpaka::pCast<double>(volume);
        }
    };

    struct ForwardPhiAseRawResult
    {
        std::vector<double> scoreSum;
        std::vector<double> scoreSquareSum;
        std::vector<unsigned> totalRays;
        std::vector<unsigned> droppedRays;
        unsigned rayCount = 0u;
    };

    [[nodiscard]] inline ForwardPhiAseRawResult makeForwardRawResult(unsigned const volumeCount)
    {
        return ForwardPhiAseRawResult{
            std::vector<double>(volumeCount, 0.0),
            std::vector<double>(volumeCount, 0.0),
            std::vector<unsigned>(volumeCount, 0u),
            std::vector<unsigned>(volumeCount, 0u),
            0u};
    }

    [[nodiscard]] inline double calcForwardBetaVolumeTotal(HostMesh const& hostMesh)
    {
        double total = 0.0;
        unsigned const count = std::min(
            static_cast<unsigned>(hostMesh.betaVolume.size()),
            static_cast<unsigned>(hostMesh.cellVolumes.size()));
        for(unsigned volume = 0u; volume < count; ++volume)
        {
            total += hostMesh.betaVolume.at(volume) * static_cast<double>(hostMesh.cellVolumes.at(volume));
        }
        return total;
    }

    inline void mergeForwardRawResult(ForwardPhiAseRawResult& target, ForwardPhiAseRawResult const& source)
    {
        if(target.scoreSum.empty())
        {
            target = source;
            return;
        }

        target.rayCount += source.rayCount;
        for(unsigned volume = 0u; volume < target.scoreSum.size(); ++volume)
        {
            target.scoreSum.at(volume) += source.scoreSum.at(volume);
            target.scoreSquareSum.at(volume) += source.scoreSquareSum.at(volume);
            target.totalRays.at(volume) += source.totalRays.at(volume);
            target.droppedRays.at(volume) += source.droppedRays.at(volume);
        }
    }

    [[nodiscard]] inline double calcForwardStandardError(
        double const scoreSum,
        double const scoreSquareSum,
        unsigned const rayCount,
        double const normalizationVolume,
        double const volumeSize)
    {
        if(rayCount < 2u || volumeSize <= 0.0 || normalizationVolume < 0.0 || !std::isfinite(scoreSum)
           || !std::isfinite(scoreSquareSum))
        {
            return std::numeric_limits<double>::max();
        }

        double const n = static_cast<double>(rayCount);
        double const varianceOfMean = (scoreSquareSum - scoreSum * scoreSum / n) / (n * (n - 1.0));
        double const volumeScale = normalizationVolume / volumeSize;
        return std::sqrt(std::max(0.0, varianceOfMean)) * volumeScale;
    }

    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    float calcForwardPhiAseRaw(
        alpakaUtils::DevBundle<T_Device, T_Exec>& devBundle,
        ExperimentParameters const& experiment,
        HostMesh const& hostMesh,
        DeviceMeshContainer<T_Device> const& meshContainer,
        ForwardPhiAseRawResult& result,
        float& runtime,
        unsigned rayCount,
        unsigned threadLocalStridingRNG)
    {
        if(experiment.forwardRayLength <= 0.0)
        {
            throw std::runtime_error("Forward volume propagation requires forwardRayLength > 0.");
        }
        if(experiment.useReflections && experiment.surfaceReservoirSize == 0u)
        {
            throw std::runtime_error("Forward reflections require surfaceReservoirSize > 0.");
        }

        time_t starttime = time(0);
        auto queue = devBundle.device.makeQueue();
        DeviceMeshView mesh = meshContainer.toView();
        unsigned const volumeCount = mesh.numberOfCells;
        double const betaVolumeTotal = calcForwardBetaVolumeTotal(hostMesh);

        result = makeForwardRawResult(volumeCount);
        result.rayCount = rayCount;
        if(rayCount == 0u)
        {
            runtime = difftime(time(0), starttime);
            return runtime;
        }

        auto dPhiAccumulator = alpaka::onHost::alloc<double>(devBundle.device, volumeCount);
        auto dPhiSquareAccumulator = alpaka::onHost::alloc<double>(devBundle.device, volumeCount);
        auto dVolumeRayVisits = alpaka::onHost::alloc<unsigned>(devBundle.device, volumeCount);
        auto dDroppedRays = alpaka::onHost::alloc<unsigned>(devBundle.device, volumeCount);
        auto dSigmaA = hase::alpakaUtils::toDevice(queue, experiment.sigmaA);
        auto dSigmaE = hase::alpakaUtils::toDevice(queue, experiment.sigmaE);
        auto accumulationSpans = hase::kernels::forward::ForwardAccumulationSpans{
            dPhiAccumulator,
            dPhiSquareAccumulator,
            dVolumeRayVisits,
            dDroppedRays};
        auto spectrumSpans = hase::kernels::forward::ForwardSpectrumSpans{
            dSigmaA,
            dSigmaE,
            static_cast<unsigned>(experiment.sigmaA.size())};

        alpaka::onHost::fill(queue, dPhiAccumulator, double{0}, alpaka::Vec{volumeCount});
        alpaka::onHost::fill(queue, dPhiSquareAccumulator, double{0}, alpaka::Vec{volumeCount});
        alpaka::onHost::fill(queue, dVolumeRayVisits, 0u, alpaka::Vec{volumeCount});
        alpaka::onHost::fill(queue, dDroppedRays, 0u, alpaka::Vec{volumeCount});
        alpaka::onHost::wait(queue);

        if(!experiment.useReflections)
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{static_cast<unsigned int>(rayCount)});
            queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::AccumulateForwardPhiAse{},
                    mesh,
                    rayCount,
                    experiment.forwardRayLength,
                    betaVolumeTotal,
                    accumulationSpans,
                    spectrumSpans,
                    threadLocalStridingRNG});
            alpaka::onHost::wait(queue);
        }
        else
        {
            unsigned const faceCount = mesh.numberOfCells * mesh.numberOfFacesPerCell;
            unsigned const reservoirSlots = faceCount * experiment.surfaceReservoirSize;
            auto countsA = alpaka::onHost::alloc<unsigned>(devBundle.device, faceCount);
            auto countsB = alpaka::onHost::alloc<unsigned>(devBundle.device, faceCount);
            auto dirXA = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
            auto dirXB = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
            auto dirYA = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
            auto dirYB = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
            auto dirZA = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
            auto dirZB = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
            auto weightsA = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
            auto weightsB = alpaka::onHost::alloc<double>(devBundle.device, reservoirSlots);
            auto sigmaIndicesA = alpaka::onHost::alloc<unsigned>(devBundle.device, reservoirSlots);
            auto sigmaIndicesB = alpaka::onHost::alloc<unsigned>(devBundle.device, reservoirSlots);
            auto totalWeightA = alpaka::onHost::alloc<double>(devBundle.device, 1u);
            auto totalWeightB = alpaka::onHost::alloc<double>(devBundle.device, 1u);
            auto reservoirSpansA = hase::kernels::forward::SurfaceReservoirSpans{
                countsA,
                dirXA,
                dirYA,
                dirZA,
                weightsA,
                sigmaIndicesA,
                totalWeightA,
                experiment.surfaceReservoirSize};
            auto reservoirSpansB = hase::kernels::forward::SurfaceReservoirSpans{
                countsB,
                dirXB,
                dirYB,
                dirZB,
                weightsB,
                sigmaIndicesB,
                totalWeightB,
                experiment.surfaceReservoirSize};

            auto clearReservoir = [&](auto& counts, auto& totalWeight)
            {
                alpaka::onHost::fill(queue, counts, 0u, alpaka::Vec{faceCount});
                alpaka::onHost::fill(queue, totalWeight, 0.0, alpaka::Vec{1u});
            };

            clearReservoir(countsA, totalWeightA);
            clearReservoir(countsB, totalWeightB);
            alpaka::onHost::wait(queue);

            auto directFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{static_cast<unsigned int>(rayCount)});
            queue.enqueue(
                directFrameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::AccumulateForwardPhiAseReservoir{},
                    mesh,
                    rayCount,
                    experiment.forwardRayLength,
                    betaVolumeTotal,
                    accumulationSpans,
                    reservoirSpansA,
                    spectrumSpans,
                    threadLocalStridingRNG});
            alpaka::onHost::wait(queue);

            std::vector<double> hostTotalWeight(1u, 0.0);
            auto reflectedFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{static_cast<unsigned int>(reservoirSlots)});
            bool inputA = true;
            for(unsigned iteration = 0u; iteration < experiment.reflectionMaxIterations; ++iteration)
            {
                alpaka::onHost::memcpy(queue, hostTotalWeight, inputA ? totalWeightA : totalWeightB);
                alpaka::onHost::wait(queue);
                double const sourceFraction = rayCount > 0u ? hostTotalWeight[0] / static_cast<double>(rayCount) : 0.0;
                if(sourceFraction < experiment.reflectionTolerance)
                {
                    break;
                }

                if(inputA)
                {
                    clearReservoir(countsB, totalWeightB);
                    alpaka::onHost::wait(queue);
                    queue.enqueue(
                        reflectedFrameSpec,
                        alpaka::KernelBundle{
                            hase::kernels::forward::AccumulateReflectedForwardPhiAse{},
                            mesh,
                            reservoirSlots,
                            experiment.forwardRayLength,
                            accumulationSpans,
                            reservoirSpansA,
                            reservoirSpansB,
                            spectrumSpans});
                }
                else
                {
                    clearReservoir(countsA, totalWeightA);
                    alpaka::onHost::wait(queue);
                    queue.enqueue(
                        reflectedFrameSpec,
                        alpaka::KernelBundle{
                            hase::kernels::forward::AccumulateReflectedForwardPhiAse{},
                            mesh,
                            reservoirSlots,
                            experiment.forwardRayLength,
                            accumulationSpans,
                            reservoirSpansB,
                            reservoirSpansA,
                            spectrumSpans});
                }
                alpaka::onHost::wait(queue);
                inputA = !inputA;
            }
        }

        alpaka::onHost::memcpy(queue, result.scoreSum, dPhiAccumulator);
        alpaka::onHost::memcpy(queue, result.scoreSquareSum, dPhiSquareAccumulator);
        alpaka::onHost::memcpy(queue, result.totalRays, dVolumeRayVisits);
        alpaka::onHost::memcpy(queue, result.droppedRays, dDroppedRays);
        alpaka::onHost::wait(queue);

        runtime = difftime(time(0), starttime);
        return runtime;
    }

    inline void finalizeForwardPhiAse(
        ExperimentParameters const& experiment,
        HostMesh const& hostMesh,
        ForwardPhiAseRawResult const& rawResult,
        Result& result)
    {
        unsigned const volumeCount = static_cast<unsigned>(rawResult.scoreSum.size());
        double const betaVolumeTotal = calcForwardBetaVolumeTotal(hostMesh);

        result = Result(
            std::vector(volumeCount, 0.0f),
            std::vector(volumeCount, 0.0),
            rawResult.totalRays,
            std::vector(volumeCount, 0.0),
            rawResult.droppedRays);
        for(unsigned volume = 0u; volume < volumeCount; ++volume)
        {
            double const volumeSize = hostMesh.cellVolumes.at(volume);
            double const scoreSum = rawResult.scoreSum.at(volume);
            if(volumeSize > 0.0 && rawResult.rayCount > 0u)
            {
                double const estimate
                    = scoreSum * betaVolumeTotal / (static_cast<double>(rawResult.rayCount) * volumeSize);
                result.phiAse.at(volume) = static_cast<float>(estimate);
                result.mse.at(volume) = result.droppedRays[volume] == 0u ? calcForwardStandardError(
                                                                               scoreSum,
                                                                               rawResult.scoreSquareSum.at(volume),
                                                                               rawResult.rayCount,
                                                                               betaVolumeTotal,
                                                                               volumeSize)
                                                                         : std::numeric_limits<double>::max();
            }
            else
            {
                result.phiAse.at(volume) = 0.0f;
                result.mse.at(volume) = std::numeric_limits<double>::max();
            }
        }
    }

    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    float calcForwardPhiAse(
        alpakaUtils::DevBundle<T_Device, T_Exec>& devBundle,
        ExperimentParameters const& experiment,
        HostMesh const& hostMesh,
        DeviceMeshContainer<T_Device> const& meshContainer,
        Result& result,
        float& runtime,
        unsigned threadLocalStridingRNG)
    {
        ForwardPhiAseRawResult rawResult;
        calcForwardPhiAseRaw(
            devBundle,
            experiment,
            hostMesh,
            meshContainer,
            rawResult,
            runtime,
            experiment.resolvedForwardRayCount(),
            threadLocalStridingRNG);
        finalizeForwardPhiAse(experiment, hostMesh, rawResult, result);
        return runtime;
    }
} // namespace hase::core
