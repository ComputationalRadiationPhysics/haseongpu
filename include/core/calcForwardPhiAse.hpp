/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/forwardSrm.hpp>
#include <core/mesh.hpp>

#include <ctime>
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

    [[nodiscard]] ForwardPhiAseRawResult makeForwardRawResult(unsigned volumeCount);

    [[nodiscard]] double calcForwardBetaVolumeTotal(HostMesh const& hostMesh);

    void mergeForwardRawResult(ForwardPhiAseRawResult& target, ForwardPhiAseRawResult const& source);

    [[nodiscard]] double calcForwardStandardError(
        double scoreSum,
        double scoreSquareSum,
        unsigned rayCount,
        double normalizationVolume,
        double volumeSize);

    void finalizeForwardPhiAse(
        HostMesh const& hostMesh,
        ForwardPhiAseRawResult const& rawResult,
        Result& result);

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
        SrmControls const srmControls = resolveSrmControls(experiment);
        result.srmMaxIterations = experiment.useReflections ? srmControls.maxIterations : 0u;
        result.srmDivergenceStreak = experiment.useReflections ? srmControls.divergenceStreak : 0u;
        if(rayCount == 0u)
        {
            result.srmStatus = experiment.useReflections ? SrmStatus::CONVERGED : SrmStatus::DISABLED;
            runtime = difftime(time(0), starttime);
            return runtime;
        }

        alpaka::concepts::IBuffer auto dPhiAccumulator = alpaka::onHost::alloc<double>(devBundle.device, volumeCount);
        alpaka::concepts::IBuffer auto dPhiSquareAccumulator
            = alpaka::onHost::alloc<double>(devBundle.device, volumeCount);
        alpaka::concepts::IBuffer auto dVolumeRayVisits = alpaka::onHost::alloc<unsigned>(devBundle.device, volumeCount);
        alpaka::concepts::IBuffer auto dDroppedRays = alpaka::onHost::alloc<unsigned>(devBundle.device, volumeCount);
        alpaka::concepts::IBuffer auto dSigmaA = hase::alpakaUtils::toDevice(queue, experiment.sigmaA);
        alpaka::concepts::IBuffer auto dSigmaE = hase::alpakaUtils::toDevice(queue, experiment.sigmaE);
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
                    betaVolumeTotal,
                    accumulationSpans,
                    spectrumSpans,
                    threadLocalStridingRNG});
            alpaka::onHost::wait(queue);
        }
        else
        {
            runForwardSrm(
                devBundle,
                queue,
                mesh,
                experiment,
                result,
                rayCount,
                betaVolumeTotal,
                dPhiAccumulator,
                dPhiSquareAccumulator,
                dVolumeRayVisits,
                dDroppedRays,
                dSigmaA,
                dSigmaE,
                static_cast<unsigned>(experiment.sigmaA.size()),
                threadLocalStridingRNG,
                srmControls);
        }

        alpaka::onHost::memcpy(queue, result.scoreSum, dPhiAccumulator);
        alpaka::onHost::memcpy(queue, result.scoreSquareSum, dPhiSquareAccumulator);
        alpaka::onHost::memcpy(queue, result.totalRays, dVolumeRayVisits);
        alpaka::onHost::memcpy(queue, result.droppedRays, dDroppedRays);
        alpaka::onHost::wait(queue);

        runtime = difftime(time(0), starttime);
        return runtime;
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
        finalizeForwardPhiAse(hostMesh, rawResult, result);
        return runtime;
    }
} // namespace hase::core
