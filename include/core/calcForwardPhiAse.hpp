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

#include <algorithm>
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

    [[nodiscard]] ForwardPhiAseRawResult makeForwardRawResult(
        unsigned volumeCount,
        unsigned vertexCount,
        unsigned batchCount = hase::kernels::forward::defaultForwardRseBatchCount);

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
        using T_CharBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<char>(std::declval<T_Device&>(), std::size_t{1}));

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
                      hase::kernels::forward::defaultForwardRseBatchCount * 2u
                          * static_cast<std::size_t>(hostMesh.numberOfMeshPoints)))
            , m_rseBatchRayCountsDevice(
                  alpaka::onHost::alloc<unsigned>(
                      m_devBundle.device,
                      static_cast<std::size_t>(hase::kernels::forward::defaultForwardRseBatchCount)))
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
            , m_betaVolumeWeights(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_betaVolumePrefixScanBuffer(
                  alpaka::onHost::alloc<char>(
                      m_devBundle.device,
                      alpaka::onHost::getScanBufferSize<double>(
                          alpaka::Vec{static_cast<std::size_t>(hostMesh.numberOfCells)})))
            , m_volumeCount(hostMesh.numberOfCells)
            , m_vertexCount(hostMesh.numberOfMeshPoints)
            , m_spectralCount(static_cast<unsigned>(experiment.sigmaA.size()))
            , m_batchCount(hase::kernels::forward::defaultForwardRseBatchCount)
            , m_rseBatchRayCounts(m_batchCount, 0u)
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

        /** @brief Resize persistent batch accumulation storage for a worker group. */
        void configureBatchCount(unsigned const batchCount)
        {
            if(batchCount == 0u)
                throw std::invalid_argument("forward ASE batch count must be positive");
            if(batchCount == m_batchCount)
                return;
            alpaka::onHost::wait(m_queue);
            m_vertexBatchScoreSum = alpaka::onHost::alloc<double>(
                m_devBundle.device,
                batchCount * 2u * static_cast<std::size_t>(m_vertexCount));
            m_rseBatchRayCountsDevice
                = alpaka::onHost::alloc<unsigned>(m_devBundle.device, static_cast<std::size_t>(batchCount));
            m_batchCount = batchCount;
            m_rseBatchRayCounts.assign(batchCount, 0u);
        }

        void begin(
            DeviceMeshView const mesh,
            unsigned rayCount,
            unsigned rngSeed,
            unsigned rseBatch,
            double betaVolumeTotal,
            ExperimentParameters const& experiment,
            bool resetAccumulators = true)
        {
            m_started = std::chrono::steady_clock::now();
            m_rayCount = rayCount;
            if(resetAccumulators)
            {
                m_accumulatedRayCount = 0u;
                std::ranges::fill(m_rseBatchRayCounts, 0u);
            }
            m_accumulatedRayCount += rayCount;
            if(rseBatch >= m_batchCount)
                throw std::out_of_range("forward ASE RSE batch index is out of range");
            m_rseBatchRayCounts.at(rseBatch) += rayCount;
            if(rayCount == 0u)
                return;

            if(resetAccumulators)
            {
                alpaka::onHost::fill(
                    m_queue,
                    m_vertexBatchScoreSum,
                    0.0,
                    alpaka::Vec{m_batchCount * 2u * static_cast<std::size_t>(m_vertexCount)});
                alpaka::onHost::fill(
                    m_queue,
                    m_volumeRayVisits,
                    0u,
                    alpaka::Vec{static_cast<std::size_t>(m_volumeCount)});
                alpaka::onHost::fill(m_queue, m_droppedRays, 0u, alpaka::Vec{static_cast<std::size_t>(m_volumeCount)});
            }

            auto accumulation = kernels::forward::ForwardAccumulationSpans{
                m_vertexBatchScoreSum.getMdSpan(),
                m_volumeRayVisits.getMdSpan(),
                m_droppedRays.getMdSpan()};
            auto spectrum = hase::kernels::forward::ForwardSpectrumSpans{m_sigmaA, m_sigmaE, m_spectralCount};
            auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                m_devBundle.device,
                m_devBundle.executor,
                alpaka::Vec{rayCount});
            m_srmResult = makeForwardRawResult(m_volumeCount, m_vertexCount, m_batchCount);
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
                    rseBatch,
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
                        rseBatch,
                        betaVolumeTotal,
                        accumulation,
                        spectrum,
                        rngSeed});
            }
        }

        void finish(ForwardPhiAseRawResult& result, float& runtime, bool downloadAccumulators = true)
        {
            result = makeForwardRawResult(m_volumeCount, m_vertexCount, m_batchCount);
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
            unsigned rseBatch,
            double betaVolumeTotal,
            ExperimentParameters const& experiment)
        {
            begin(mesh, rayCount, rngSeed, rseBatch, betaVolumeTotal, experiment);
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
            alpaka::onHost::memcpy(m_queue, m_rseBatchRayCountsDevice, m_rseBatchRayCounts);
            hase::kernels::enqueueFinalizeForwardCellPhiAse(
                m_devBundle,
                m_queue,
                mesh,
                m_vertexBatchScoreSum,
                m_rseBatchRayCountsDevice,
                m_lumpedMaterialVertexVolume,
                m_droppedRays,
                m_volumePhiAse,
                m_standardError,
                m_relativeStandardError,
                m_volumeDndtAse,
                rayCount,
                m_batchCount,
                betaVolumeTotal,
                fluorescenceRate,
                sigmaA,
                sigmaE);
            alpaka::onHost::wait(m_queue);
        }

        /**
         * @brief Upload gathered raw batches and finalize them on this device.
         *
         * This is the explicit host/device boundary after a thread or MPI
         * gather. Normalization, RSE evaluation, and ASE derivative generation
         * remain device-side.
         */
        void uploadAndFinalize(
            DeviceMeshView const mesh,
            ForwardPhiAseRawResult const& rawResult,
            double const betaVolumeTotal,
            double const fluorescenceRate,
            double const sigmaA,
            double const sigmaE)
        {
            if(rawResult.rseBatchRayCounts.size() != m_batchCount
               || rawResult.vertexBatchScoreSum.size() != m_vertexBatchScoreSum.getExtents().product()
               || rawResult.totalRays.size() != m_volumeCount || rawResult.droppedRays.size() != m_volumeCount)
                throw std::runtime_error("gathered forward ASE result does not match the device context");
            alpaka::onHost::memcpy(m_queue, m_vertexBatchScoreSum, rawResult.vertexBatchScoreSum);
            alpaka::onHost::memcpy(m_queue, m_volumeRayVisits, rawResult.totalRays);
            alpaka::onHost::memcpy(m_queue, m_droppedRays, rawResult.droppedRays);
            m_rseBatchRayCounts = rawResult.rseBatchRayCounts;
            m_accumulatedRayCount = rawResult.rayCount;
            finalizeCellPhiAse(mesh, rawResult.rayCount, betaVolumeTotal, fluorescenceRate, sigmaA, sigmaE);
        }

        double rebuildBetaVolumePrefix(DeviceMeshContainer<T_Device>& meshContainer, auto const& betaVolume)
        {
            auto mesh = meshContainer.toView();
            mesh.betaVolume = std::span<double const>(betaVolume.data(), betaVolume.getMdSpan().getExtents().x());
            if(mesh.numberOfCells == 0u)
            {
                alpaka::onHost::fill(m_queue, m_betaVolumeTotal, 0.0, alpaka::Vec{std::size_t{1}});
            }
            else
            {
                auto const cellFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                    m_devBundle.device,
                    m_devBundle.executor,
                    alpaka::Vec{mesh.numberOfCells});
                m_queue.enqueue(
                    cellFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::BuildBetaVolumeWeights{},
                        mesh,
                        betaVolume,
                        m_betaVolumeWeights});
                alpaka::onHost::inclusiveScan(
                    m_queue,
                    m_devBundle.executor,
                    m_betaVolumePrefixScanBuffer,
                    meshContainer.betaVolumePrefix,
                    m_betaVolumeWeights);
                auto const scalarFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
                    m_devBundle.device,
                    m_devBundle.executor,
                    alpaka::Vec{1u});
                m_queue.enqueue(
                    scalarFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::CaptureBetaVolumeTotal{},
                        mesh.numberOfCells,
                        meshContainer.betaVolumePrefix,
                        m_betaVolumeTotal});
            }
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
        T_UnsignedBuffer m_rseBatchRayCountsDevice;
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
        T_DoubleBuffer m_betaVolumeWeights;
        T_CharBuffer m_betaVolumePrefixScanBuffer;
        std::unique_ptr<ForwardSrmWorkspace<T_Device>> m_srmWorkspace;
        ForwardPhiAseRawResult m_srmResult;
        unsigned m_volumeCount;
        unsigned m_vertexCount;
        unsigned m_spectralCount;
        unsigned m_batchCount;
        unsigned m_rayCount = 0u;
        unsigned m_accumulatedRayCount = 0u;
        std::vector<unsigned> m_rseBatchRayCounts;
        std::chrono::steady_clock::time_point m_started;
    };

} // namespace hase::core
