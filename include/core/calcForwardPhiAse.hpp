/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <benchmark.hpp>
#include <core/forwardSrm.hpp>
#include <core/mesh.hpp>
#include <kernels/forwardPhiAseMapping.hpp>
#include <kernels/vertexAccumulation.hpp>
#include <random/random.hpp>

#include <array>
#include <chrono>
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

    [[nodiscard]] ForwardPhiAseRawResult makeForwardRawResult(unsigned volumeCount, unsigned vertexCount);

    [[nodiscard]] double calcForwardBetaVolumeTotal(HostMesh const& hostMesh);

    void mergeForwardRawResult(ForwardPhiAseRawResult& target, ForwardPhiAseRawResult const& source);

    [[nodiscard]] double calcForwardRelativeStandardError(double scoreSum, double scoreSquareSum, unsigned rayCount);

    [[nodiscard]] double calcForwardStandardError(
        double scoreSum,
        double scoreSquareSum,
        unsigned rayCount,
        double normalizationVolume,
        double volumeSize);

    void finalizeForwardPhiAse(HostMesh const& hostMesh, ForwardPhiAseRawResult const& rawResult, Result& result);

    void finalizeForwardPhiAse(
        HostMesh const& hostMesh,
        ForwardPhiAseRawResult const& rawResult,
        double betaVolumeTotal,
        Result& result);

    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    class ForwardPhiAseDeviceContext
    {
        using T_Queue = ALPAKA_TYPEOF(std::declval<T_Device>().makeQueue(alpaka::queueKind::nonBlocking));
        using T_DoubleBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<double>(std::declval<T_Device&>(), std::size_t{1}));
        using T_FloatBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<float>(std::declval<T_Device&>(), std::size_t{1}));
        using T_UnsignedBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<unsigned>(std::declval<T_Device&>(), std::size_t{1}));

    public:
        ForwardPhiAseDeviceContext(
            T_Device const& device,
            T_Exec const& executor,
            ExperimentParameters const& experiment,
            HostMesh const& hostMesh)
            : m_devBundle(device, executor)
            , m_queue(m_devBundle.device.makeQueue(alpaka::queueKind::nonBlocking))
            , m_vertexBatchScoreSum(
                  alpaka::onHost::alloc<double>(
                      m_devBundle.device,
                      hase::kernels::forward::forwardRseBatchCount * 2u
                          * static_cast<std::size_t>(hostMesh.numberOfMeshPoints)))
            , m_volumeRayVisits(
                  alpaka::onHost::alloc<unsigned>(
                      m_devBundle.device,
                      static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_droppedRays(
                  alpaka::onHost::alloc<unsigned>(
                      m_devBundle.device,
                      static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_sigmaA(hase::alpakaUtils::toDevice(m_queue, experiment.sigmaA))
            , m_sigmaE(hase::alpakaUtils::toDevice(m_queue, experiment.sigmaE))
            , m_lumpedMaterialVertexVolume(
                  hase::alpakaUtils::toDevice(m_queue, hase::kernels::makeLumpedMaterialVertexVolumes(hostMesh)))
            , m_volumePhiAse(
                  alpaka::onHost::alloc<float>(m_devBundle.device, static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_standardError(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_relativeStandardError(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_volumeDndtAse(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_betaVolumeTotal(alpaka::onHost::alloc<double>(m_devBundle.device, std::size_t{1}))
            , m_volumeCount(hostMesh.numberOfCells)
            , m_vertexCount(hostMesh.numberOfMeshPoints)
            , m_spectralCount(static_cast<unsigned>(experiment.sigmaA.size()))
        {
            if(experiment.useReflections)
            {
                m_srmWorkspace = std::make_unique<ForwardSrmWorkspace<T_Device>>(
                    m_devBundle.device,
                    m_volumeCount * tet4FaceCount,
                    experiment.surfaceReservoirSize,
                    std::max(experiment.maxRays, experiment.resolvedForwardRayCount()));
            }
        }

        void begin(
            DeviceMeshView const mesh,
            unsigned rayCount,
            unsigned rngSeed,
            unsigned globalRayOffset,
            unsigned globalRayCount,
            double sourceStratificationOffset,
            unsigned spectrumStratificationPhase,
            double betaVolumeTotal,
            ExperimentParameters const& experiment,
            bool resetAccumulators = true)
        {
            m_started = std::chrono::steady_clock::now();
            m_rayCount = rayCount;
            if(resetAccumulators)
            {
                m_accumulatedRayCount = 0u;
                m_rseBatchRayCounts.fill(0u);
            }
            m_accumulatedRayCount += rayCount;
            for(unsigned batch = 0u; batch < hase::kernels::forward::forwardRseBatchCount; ++batch)
            {
                m_rseBatchRayCounts.at(batch)
                    += hase::kernels::forward::rseBatchRayCount(globalRayOffset, rayCount, batch);
            }
            if(rayCount == 0u)
                return;

            if(resetAccumulators)
            {
                alpaka::onHost::fill(
                    m_queue,
                    m_vertexBatchScoreSum,
                    0.0,
                    alpaka::Vec{
                        hase::kernels::forward::forwardRseBatchCount * 2u * static_cast<std::size_t>(m_vertexCount)});
                alpaka::onHost::fill(
                    m_queue,
                    m_volumeRayVisits,
                    0u,
                    alpaka::Vec{static_cast<std::size_t>(m_volumeCount)});
                alpaka::onHost::fill(m_queue, m_droppedRays, 0u, alpaka::Vec{static_cast<std::size_t>(m_volumeCount)});
            }

            auto accumulation = hase::kernels::forward::ForwardAccumulationSpans{
                m_vertexBatchScoreSum.getMdSpan(),
                m_volumeRayVisits.getMdSpan(),
                m_droppedRays.getMdSpan()};
            auto spectrum = hase::kernels::forward::ForwardSpectrumSpans{m_sigmaA, m_sigmaE, m_spectralCount};
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{rayCount});
            m_srmResult = makeForwardRawResult(m_volumeCount, m_vertexCount);
            m_srmResult.rayCount = rayCount;
            if(experiment.useReflections)
            {
                if(!m_srmWorkspace)
                    throw std::runtime_error("persistent forward SRM workspace was not initialized");
                auto const controls = resolveSrmControls(experiment);
                m_srmResult.srmMaxIterations = controls.maxIterations;
                m_srmResult.srmDivergenceStreak = controls.divergenceStreak;
                runForwardSrm(
                    m_devBundle,
                    m_queue,
                    mesh,
                    experiment,
                    m_srmResult,
                    rayCount,
                    globalRayOffset,
                    globalRayCount,
                    sourceStratificationOffset,
                    spectrumStratificationPhase,
                    betaVolumeTotal,
                    m_vertexBatchScoreSum,
                    m_volumeRayVisits,
                    m_droppedRays,
                    m_sigmaA,
                    m_sigmaE,
                    m_spectralCount,
                    rngSeed,
                    controls,
                    *m_srmWorkspace);
            }
            else
            {
                BENCH_SYNC(m_queue, AccumulateForwardPhiAse);
                m_queue.enqueue(
                    frameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::forward::AccumulateForwardPhiAse{},
                        mesh,
                        rayCount,
                        globalRayOffset,
                        globalRayCount,
                        sourceStratificationOffset,
                        spectrumStratificationPhase,
                        betaVolumeTotal,
                        accumulation,
                        spectrum,
                        rngSeed});
            }
        }

        void finish(ForwardPhiAseRawResult& result, float& runtime, bool downloadAccumulators = true)
        {
            result = makeForwardRawResult(m_volumeCount, m_vertexCount);
            result.rayCount = m_accumulatedRayCount;
            result.rseBatchRayCounts = m_rseBatchRayCounts;
            result.srmStatus = m_srmResult.srmStatus;
            result.srmPasses = m_srmResult.srmPasses;
            result.srmRemainingFraction = m_srmResult.srmRemainingFraction;
            result.srmMaxIterations = m_srmResult.srmMaxIterations;
            result.srmDivergenceStreak = m_srmResult.srmDivergenceStreak;
            if(m_rayCount == 0u)
            {
                runtime = 0.0f;
                return;
            }
            if(downloadAccumulators)
            {
                alpaka::onHost::memcpy(m_queue, result.vertexBatchScoreSum, m_vertexBatchScoreSum);
                alpaka::onHost::memcpy(m_queue, result.totalRays, m_volumeRayVisits);
                alpaka::onHost::memcpy(m_queue, result.droppedRays, m_droppedRays);
            }
            alpaka::onHost::wait(m_queue);
            runtime = static_cast<float>(
                std::chrono::duration<double>(std::chrono::steady_clock::now() - m_started).count());
        }

        void evaluate(
            DeviceMeshView const mesh,
            ForwardPhiAseRawResult& result,
            float& runtime,
            unsigned rayCount,
            unsigned rngSeed,
            unsigned globalRayOffset,
            unsigned globalRayCount,
            double sourceStratificationOffset,
            unsigned spectrumStratificationPhase,
            double betaVolumeTotal,
            ExperimentParameters const& experiment)
        {
            begin(
                mesh,
                rayCount,
                rngSeed,
                globalRayOffset,
                globalRayCount,
                sourceStratificationOffset,
                spectrumStratificationPhase,
                betaVolumeTotal,
                experiment);
            finish(result, runtime);
        }

        void finalizeCellPhiAse(
            DeviceMeshView const mesh,
            unsigned rayCount,
            double betaVolumeTotal,
            double fluorescenceRate,
            double sigmaA,
            double sigmaE)
        {
            alpaka::Vec<unsigned, hase::kernels::forward::forwardRseBatchCount> rseBatchRayCounts{};
            for(unsigned batch = 0u; batch < hase::kernels::forward::forwardRseBatchCount; ++batch)
                rseBatchRayCounts[batch] = m_rseBatchRayCounts.at(batch);
            hase::kernels::enqueueFinalizeForwardCellPhiAse(
                m_devBundle,
                m_queue,
                mesh,
                m_vertexBatchScoreSum,
                m_lumpedMaterialVertexVolume,
                m_droppedRays,
                m_volumePhiAse,
                m_standardError,
                m_relativeStandardError,
                m_volumeDndtAse,
                rayCount,
                rseBatchRayCounts,
                betaVolumeTotal,
                fluorescenceRate,
                sigmaA,
                sigmaE);
            alpaka::onHost::wait(m_queue);
        }

        double rebuildBetaVolumePrefix(DeviceMeshContainer<T_Device>& meshContainer, auto const& betaVolume)
        {
            auto mesh = meshContainer.toView();
            mesh.betaVolume = std::span<double const>(betaVolume.data(), betaVolume.getMdSpan().getExtents().x());
            hase::kernels::enqueueBuildBetaVolumePrefix(
                m_devBundle,
                m_queue,
                mesh,
                betaVolume,
                meshContainer.betaVolumePrefix,
                m_betaVolumeTotal);
            std::vector<double> hostTotal(1u, 0.0);
            alpaka::onHost::memcpy(m_queue, hostTotal, m_betaVolumeTotal);
            alpaka::onHost::wait(m_queue);
            return hostTotal.front();
        }

        std::vector<double> downloadVolumeDndtAse()
        {
            std::vector<double> result(m_volumeCount, 0.0);
            alpaka::onHost::memcpy(m_queue, result, m_volumeDndtAse);
            alpaka::onHost::wait(m_queue);
            return result;
        }

        Result downloadFinalizedResult(
            bool includePhiAse,
            bool includeStandardError,
            bool includeRelativeStandardError,
            bool includeTotalRays)
        {
            Result result;
            if(includePhiAse)
            {
                result.phiAse.resize(m_volumeCount);
                alpaka::onHost::memcpy(m_queue, result.phiAse, m_volumePhiAse);
            }
            if(includeStandardError)
            {
                result.standardError.resize(m_volumeCount);
                alpaka::onHost::memcpy(m_queue, result.standardError, m_standardError);
            }
            if(includeRelativeStandardError)
            {
                result.relativeStandardError.resize(m_volumeCount);
                alpaka::onHost::memcpy(m_queue, result.relativeStandardError, m_relativeStandardError);
            }
            if(includeTotalRays)
            {
                result.totalRays.resize(m_volumeCount);
                result.droppedRays.resize(m_volumeCount);
                alpaka::onHost::memcpy(m_queue, result.totalRays, m_volumeRayVisits);
                alpaka::onHost::memcpy(m_queue, result.droppedRays, m_droppedRays);
            }
            alpaka::onHost::wait(m_queue);
            return result;
        }

        [[nodiscard]] auto& volumePhiAse()
        {
            return m_volumePhiAse;
        }

        [[nodiscard]] auto& volumeDndtAse()
        {
            return m_volumeDndtAse;
        }

    private:
        hase::alpakaUtils::DevBundle<T_Device, T_Exec> m_devBundle;
        T_Queue m_queue;
        T_DoubleBuffer m_vertexBatchScoreSum;
        T_UnsignedBuffer m_volumeRayVisits;
        T_UnsignedBuffer m_droppedRays;
        T_DoubleBuffer m_sigmaA;
        T_DoubleBuffer m_sigmaE;
        T_DoubleBuffer m_lumpedMaterialVertexVolume;
        T_FloatBuffer m_volumePhiAse;
        T_DoubleBuffer m_standardError;
        T_DoubleBuffer m_relativeStandardError;
        T_DoubleBuffer m_volumeDndtAse;
        T_DoubleBuffer m_betaVolumeTotal;
        std::unique_ptr<ForwardSrmWorkspace<T_Device>> m_srmWorkspace;
        ForwardPhiAseRawResult m_srmResult;
        unsigned m_volumeCount;
        unsigned m_vertexCount;
        unsigned m_spectralCount;
        unsigned m_rayCount = 0u;
        unsigned m_accumulatedRayCount = 0u;
        std::array<unsigned, hase::kernels::forward::forwardRseBatchCount> m_rseBatchRayCounts{};
        std::chrono::steady_clock::time_point m_started;
    };

    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    float calcForwardPhiAseRaw(
        alpakaUtils::DevBundle<T_Device, T_Exec>& devBundle,
        ExperimentParameters const& experiment,
        HostMesh const& hostMesh,
        DeviceMeshContainer<T_Device> const& meshContainer,
        ForwardPhiAseRawResult& result,
        float& runtime,
        unsigned rayCount,
        unsigned threadLocalStridingRNG,
        unsigned const globalRayOffset = 0u,
        unsigned const globalRayCount = 0u,
        double const sourceStratificationOffset = 0.0,
        unsigned const spectrumStratificationPhase = 0u)
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
        unsigned const resolvedGlobalRayCount = globalRayCount == 0u ? rayCount : globalRayCount;
        double const resolvedSourceStratificationOffset
            = globalRayCount == 0u ? random::stratifiedUnitOffset(threadLocalStridingRNG) : sourceStratificationOffset;
        unsigned const resolvedSpectrumStratificationPhase = globalRayCount == 0u
                                                                 ? random::stratifiedSpectrumPhase(
                                                                       threadLocalStridingRNG,
                                                                       static_cast<unsigned>(experiment.sigmaA.size()))
                                                                 : spectrumStratificationPhase;

        result = makeForwardRawResult(volumeCount, mesh.numberOfMeshPoints);
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

        alpaka::concepts::IBuffer auto dVertexBatchScoreSum = alpaka::onHost::alloc<double>(
            devBundle.device,
            hase::kernels::forward::forwardRseBatchCount * 2u * mesh.numberOfMeshPoints);
        alpaka::concepts::IBuffer auto dVolumeRayVisits
            = alpaka::onHost::alloc<unsigned>(devBundle.device, volumeCount);
        alpaka::concepts::IBuffer auto dDroppedRays = alpaka::onHost::alloc<unsigned>(devBundle.device, volumeCount);
        alpaka::concepts::IBuffer auto dSigmaA = hase::alpakaUtils::toDevice(queue, experiment.sigmaA);
        alpaka::concepts::IBuffer auto dSigmaE = hase::alpakaUtils::toDevice(queue, experiment.sigmaE);
        auto accumulationSpans = hase::kernels::forward::ForwardAccumulationSpans{
            dVertexBatchScoreSum.getMdSpan(),
            dVolumeRayVisits.getMdSpan(),
            dDroppedRays.getMdSpan()};
        auto spectrumSpans = hase::kernels::forward::ForwardSpectrumSpans{
            dSigmaA,
            dSigmaE,
            static_cast<unsigned>(experiment.sigmaA.size())};

        alpaka::onHost::fill(
            queue,
            dVertexBatchScoreSum,
            double{0},
            alpaka::Vec{hase::kernels::forward::forwardRseBatchCount * 2u * mesh.numberOfMeshPoints});
        alpaka::onHost::fill(queue, dVolumeRayVisits, 0u, alpaka::Vec{volumeCount});
        alpaka::onHost::fill(queue, dDroppedRays, 0u, alpaka::Vec{volumeCount});
        alpaka::onHost::wait(queue);

        if(!experiment.useReflections)
        {
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{static_cast<unsigned int>(rayCount)});
            BENCH_SYNC(queue, AccumulateForwardPhiAse);
            queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::AccumulateForwardPhiAse{},
                    mesh,
                    rayCount,
                    globalRayOffset,
                    resolvedGlobalRayCount,
                    resolvedSourceStratificationOffset,
                    resolvedSpectrumStratificationPhase,
                    betaVolumeTotal,
                    accumulationSpans,
                    spectrumSpans,
                    threadLocalStridingRNG});
            alpaka::onHost::wait(queue);
        }
        else
        {
            ForwardSrmWorkspace<T_Device> workspace{
                devBundle.device,
                volumeCount * tet4FaceCount,
                experiment.surfaceReservoirSize,
                rayCount};
            runForwardSrm(
                devBundle,
                queue,
                mesh,
                experiment,
                result,
                rayCount,
                globalRayOffset,
                resolvedGlobalRayCount,
                resolvedSourceStratificationOffset,
                resolvedSpectrumStratificationPhase,
                betaVolumeTotal,
                dVertexBatchScoreSum,
                dVolumeRayVisits,
                dDroppedRays,
                dSigmaA,
                dSigmaE,
                static_cast<unsigned>(experiment.sigmaA.size()),
                threadLocalStridingRNG,
                srmControls,
                workspace);
        }

        alpaka::onHost::memcpy(queue, result.vertexBatchScoreSum, dVertexBatchScoreSum);
        alpaka::onHost::memcpy(queue, result.totalRays, dVolumeRayVisits);
        alpaka::onHost::memcpy(queue, result.droppedRays, dDroppedRays);
        alpaka::onHost::wait(queue);

        for(unsigned batch = 0u; batch < hase::kernels::forward::forwardRseBatchCount; ++batch)
            result.rseBatchRayCounts.at(batch)
                = hase::kernels::forward::rseBatchRayCount(globalRayOffset, rayCount, batch);

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
