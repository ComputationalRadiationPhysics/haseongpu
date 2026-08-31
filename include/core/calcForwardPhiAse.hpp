/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpakaUtils/HybridBuffer.hpp>
#include <alpakaUtils/memory.hpp>
#include <benchmark.hpp>
#include <core/forwardSrm.hpp>
#include <data/TraceData.hpp>
#include <kernels/forwardPhiAseMapping.hpp>
#include <kernels/vertexAccumulation.hpp>
#include <random/random.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <stdexcept>
#include <vector>

namespace hase::core
{
    /** @brief Multiply excitation fraction by cell volume for a device reduction. */
    struct BetaVolumeContribution
    {
        constexpr auto operator()(alpaka::concepts::Simd auto const& beta, alpaka::concepts::Simd auto const& volume)
            const
        {
            return beta * alpaka::pCast<double>(volume);
        }
    };

    /**
     * @brief Allocate zero-initialized host accumulators for a forward trace.
     * @param volumeCount Number of cells represented by volume accumulators.
     * @param vertexCount Number of material-vertex accumulation entries per batch.
     * @param batchCount Number of independent RSE batches.
     * @return Raw result with consistently sized, zero-filled arrays.
     */
    [[nodiscard]] ForwardPhiAseRawResult makeForwardRawResult(
        unsigned volumeCount,
        unsigned vertexCount,
        unsigned batchCount = kernels::forward::defaultForwardRseBatchCount);

    /**
     * @param trace Prepared host trace containing the cumulative source-strength prefix.
     * @return Total spontaneous source strength, or zero for an empty trace.
     */
    [[nodiscard]] double calcForwardSourceStrengthTotal(hase::data::TraceData const& trace);

    /**
     * @brief Add one raw batch result into a compatible aggregate.
     * @param target Aggregate updated in place.
     * @param source Raw accumulators and counters to add.
     * @throws std::invalid_argument If the accumulator layouts differ.
     */
    void mergeForwardRawResult(ForwardPhiAseRawResult& target, ForwardPhiAseRawResult const& source);

    /**
     * @param scoreSum Sum of per-ray scores for one cell.
     * @param scoreSquareSum Sum of squared per-ray scores for the same cell.
     * @param rayCount Number of histories represented by the moments.
     * @return Relative standard error, or infinity when it is undefined.
     */
    [[nodiscard]] double calcForwardRelativeStandardError(double scoreSum, double scoreSquareSum, unsigned rayCount);

    /**
     * @param scoreSum Sum of per-ray scores for one cell.
     * @param scoreSquareSum Sum of squared per-ray scores for the same cell.
     * @param rayCount Number of histories represented by the moments.
     * @param normalizationVolume Total excited-volume normalization.
     * @param volumeSize Physical volume of the receiving cell.
     * @return Absolute standard error after physical normalization.
     */
    [[nodiscard]] double calcForwardStandardError(
        double scoreSum,
        double scoreSquareSum,
        unsigned rayCount,
        double normalizationVolume,
        double volumeSize);

    /**
     * @brief Normalize raw accumulators using the source total stored in `hostMesh`.
     * @param hostMesh Prepared geometry, material, and excitation data.
     * @param rawResult Unnormalized batch accumulators.
     * @param result Final cell-ordered values replaced in place.
     */
    void finalizeForwardPhiAse(
        hase::data::TraceData const& hostMesh,
        ForwardPhiAseRawResult const& rawResult,
        data::PhiAseResult& result);

    /**
     * @brief Normalize raw accumulators with an explicitly supplied source total.
     * @param hostMesh Prepared geometry and material data.
     * @param rawResult Unnormalized batch accumulators.
     * @param betaVolumeTotal Total excitation-volume normalization.
     * @param result Final cell-ordered values replaced in place.
     */
    void finalizeForwardPhiAse(
        hase::data::TraceData const& hostMesh,
        ForwardPhiAseRawResult const& rawResult,
        double betaVolumeTotal,
        data::PhiAseResult& result);

    /** @brief Persistent device allocations and queue for one forward ASE worker. */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    class ForwardPhiAseDeviceContext
    {
        using T_Queue = ALPAKA_TYPEOF(std::declval<T_Device>().makeQueue(alpaka::queueKind::nonBlocking));
        using T_DoubleBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<double>(std::declval<T_Device&>(), std::size_t{1}));
        using T_FloatBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<float>(std::declval<T_Device&>(), std::size_t{1}));
        using T_UnsignedBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<std::uint32_t>(std::declval<T_Device&>(), std::size_t{1}));
        using T_ByteBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<char>(std::declval<T_Device&>(), std::size_t{1}));
        using T_RseBatchRayCountsBuffer = alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<std::uint32_t>>;
        using T_BetaVolumeTotalBuffer = alpakaUtils::GetHybridBuffer_t<T_Device, std::array<double, 1u>>;

    public:
        /**
         * @param device Device that owns all persistent buffers and the internal queue.
         * @param executor Executor used for trace, reduction, and finalization kernels.
         * @param experiment Controls used to size optional reflection storage.
         * @param hostMesh Prepared trace whose dimensions determine buffer extents.
         */
        ForwardPhiAseDeviceContext(
            T_Device const& device,
            T_Exec const& executor,
            AseTraceControls const& experiment,
            data::TraceData const& hostMesh)
            : m_devBundle(device, executor)
            , m_queue(m_devBundle.device.makeQueue(alpaka::queueKind::nonBlocking))
            , m_rseBatchRayCounts(kernels::forward::defaultForwardRseBatchCount, 0u)
            , m_rseBatchRayCountsBuffer(alpakaUtils::getHybridBuffer(m_devBundle.device, m_rseBatchRayCounts))
            , m_sourceStrengthTotalHost{}
            , m_sourceStrengthTotal(alpakaUtils::getHybridBuffer(m_devBundle.device, m_sourceStrengthTotalHost))
            , m_vertexBatchScoreSum(
                  alpaka::onHost::alloc<double>(
                      m_devBundle.device,
                      hase::kernels::forward::defaultForwardRseBatchCount * hostMesh.numberOfMaterials
                          * static_cast<std::size_t>(hostMesh.numberOfMeshPoints)))
            , m_volumeRayVisits(
                  alpaka::onHost::alloc<std::uint32_t>(
                      m_devBundle.device,
                      static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_droppedRays(
                  alpaka::onHost::alloc<std::uint32_t>(
                      m_devBundle.device,
                      static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_volumePhiAse(
                  alpaka::onHost::alloc<float>(m_devBundle.device, static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_standardError(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_relativeStandardError(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_volumeDndtAse(
                  alpaka::onHost::alloc<double>(m_devBundle.device, static_cast<std::size_t>(hostMesh.numberOfCells)))
            , m_sourceStrengthPrefixScanBuffer(
                  alpaka::onHost::alloc<char>(
                      m_devBundle.device,
                      alpaka::onHost::getScanBufferSize<double>(
                          alpaka::Vec{static_cast<std::size_t>(hostMesh.numberOfCells)})))
            , m_volumeCount(hostMesh.numberOfCells)
            , m_materialVertexCount(hostMesh.numberOfMaterials * hostMesh.numberOfMeshPoints)
            , m_batchCount(hase::kernels::forward::defaultForwardRseBatchCount)
        {
            if(experiment.useReflections)
            {
                std::uint32_t const maxRayCount = std::max(experiment.maxRays, experiment.resolvedForwardRayCount());
                if(experiment.reflectionMode == "direct")
                {
                    m_directReflectionScratch
                        = std::make_unique<ReflectionResamplingScratch<T_Device>>(m_devBundle.device, maxRayCount);
                }
                else
                {
                    std::uint32_t const faceCount = hostMesh.numberOfCells * hostMesh.numberOfFacesPerCell;
                    auto const slotsPerFace = alpaka::Vec{experiment.surfaceReservoirSize};
                    if(experiment.srmPositionMode == "centroid")
                    {
                        m_centroidSurfaceReservoirScratch = std::make_unique<CentroidSurfaceReservoirScratch>(
                            m_devBundle.device,
                            faceCount,
                            slotsPerFace,
                            maxRayCount);
                    }
                    else
                    {
                        m_exactSurfaceReservoirScratch = std::make_unique<ExactSurfaceReservoirScratch>(
                            m_devBundle.device,
                            faceCount,
                            slotsPerFace,
                            maxRayCount);
                    }
                }
            }
        }

        /**
         * @brief Resize persistent batch accumulation storage for a worker group.
         * @param batchCount Positive number of independent statistical batches.
         */
        void configureBatchCount(unsigned const batchCount)
        {
            if(batchCount == 0u)
                throw std::invalid_argument("forward ASE batch count must be positive");
            if(batchCount == m_batchCount)
                return;
            alpaka::onHost::wait(m_queue);
            m_vertexBatchScoreSum = alpaka::onHost::alloc<double>(
                m_devBundle.device,
                batchCount * static_cast<std::size_t>(m_materialVertexCount));
            m_batchCount = batchCount;
            m_rseBatchRayCounts.assign(batchCount, 0u);
            m_rseBatchRayCountsBuffer = hase::alpakaUtils::getHybridBuffer(m_devBundle.device, m_rseBatchRayCounts);
        }

        /**
         * @brief Reset or extend accumulators and enqueue one forward-ray batch.
         * @param mesh Device-resident trace view.
         * @param rayCount Number of histories in this launch.
         * @param rngSeed Seed for this adaptive launch.
         * @param rseBatch Statistical batch receiving the scores.
         * @param betaVolumeTotal Current total excitation-volume normalization.
         * @param experiment Trace and reflection controls.
         * @param resetAccumulators Whether to clear results before this launch.
         */
        void begin(
            hase::data::TraceView const mesh,
            unsigned rayCount,
            unsigned rngSeed,
            unsigned rseBatch,
            double betaVolumeTotal,
            AseTraceControls const& experiment,
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
                    alpaka::Vec{m_batchCount * static_cast<std::size_t>(m_materialVertexCount)});
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
            auto frameSpec
                = alpakaUtils::getFrameSpec<uint32_t>(m_devBundle.device, m_devBundle.executor, alpaka::Vec{rayCount});
            m_srmResult = makeForwardRawResult(m_volumeCount, m_materialVertexCount, m_batchCount);
            m_srmResult.rayCount = rayCount;
            if(experiment.useReflections)
            {
                auto const controls = resolveSrmControls(experiment);
                m_srmResult.srmMaxIterations = controls.maxIterations;
                m_srmResult.srmDivergenceStreak = controls.divergenceStreak;
                if(experiment.reflectionMode == "direct")
                {
                    if(!m_directReflectionScratch)
                        throw std::runtime_error("direct reflection scratch was not initialized");
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
                        rngSeed,
                        controls,
                        *m_directReflectionScratch);
                }
                else
                {
                    auto runSurfaceSrm = [&](auto& scratch)
                    {
                        runForwardSurfaceSrm(
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
                            rngSeed,
                            controls,
                            scratch);
                    };
                    if(experiment.srmPositionMode == "centroid")
                    {
                        if(!m_centroidSurfaceReservoirScratch)
                            throw std::runtime_error("centroid surface reservoir scratch was not initialized");
                        runSurfaceSrm(*m_centroidSurfaceReservoirScratch);
                    }
                    else
                    {
                        if(!m_exactSurfaceReservoirScratch)
                            throw std::runtime_error("exact surface reservoir scratch was not initialized");
                        runSurfaceSrm(*m_exactSurfaceReservoirScratch);
                    }
                }
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
                        rngSeed});
            }
        }

        /**
         * @brief Wait for the current launch and optionally download raw accumulators.
         * @param result Host result replaced with current counters and SRM metadata.
         * @param runtime Wall-clock seconds since `begin`, replaced on return.
         * @param downloadAccumulators Whether to copy score, visit, and dropped-ray arrays.
         */
        void finish(ForwardPhiAseRawResult& result, float& runtime, bool downloadAccumulators = true)
        {
            result = makeForwardRawResult(m_volumeCount, m_materialVertexCount, m_batchCount);
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

        /**
         * @brief Enqueue and finish one independently reset statistical batch.
         * @param mesh Device-resident trace view.
         * @param result Downloaded raw result.
         * @param runtime Measured wall-clock seconds.
         * @param rayCount Number of histories to trace.
         * @param rngSeed Adaptive-launch seed.
         * @param rseBatch Statistical batch index.
         * @param betaVolumeTotal Current total excitation-volume normalization.
         * @param experiment Trace and reflection controls.
         */
        void evaluate(
            data::TraceView const mesh,
            ForwardPhiAseRawResult& result,
            float& runtime,
            unsigned rayCount,
            unsigned rngSeed,
            unsigned rseBatch,
            double betaVolumeTotal,
            AseTraceControls const& experiment)
        {
            begin(mesh, rayCount, rngSeed, rseBatch, betaVolumeTotal, experiment);
            finish(result, runtime);
        }

        /**
         * @brief Normalize resident raw accumulators into cell result buffers.
         * @param mesh Device-resident trace view.
         * @param rayCount Total histories represented by the raw accumulators.
         * @param sourceStrengthTotal Total source-strength normalization.
         */
        void finalizeCellPhiAse(data::TraceView const mesh, unsigned rayCount, double sourceStrengthTotal)
        {
            m_rseBatchRayCountsBuffer.toDevice(m_queue);
            kernels::enqueueFinalizeForwardCellPhiAse(
                m_devBundle,
                m_queue,
                mesh,
                m_vertexBatchScoreSum,
                m_rseBatchRayCountsBuffer.toDeviceView(),
                m_droppedRays,
                m_volumePhiAse,
                m_standardError,
                m_relativeStandardError,
                m_volumeDndtAse,
                rayCount,
                m_batchCount,
                sourceStrengthTotal);
            alpaka::onHost::wait(m_queue);
        }

        /**
         * @brief Upload gathered raw batches and finalize them on this device.
         *
         * This is the explicit host/device boundary after a thread or MPI
         * gather. Normalization, RSE evaluation, and ASE derivative generation
         * remain device-side.
         * @param mesh Device-resident trace view used for finalization.
         * @param rawResult Gathered host accumulators to upload.
         * @param sourceStrengthTotal Total source-strength normalization.
         */
        void uploadAndFinalize(
            data::TraceView const mesh,
            ForwardPhiAseRawResult const& rawResult,
            double const sourceStrengthTotal)
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
            finalizeCellPhiAse(mesh, rawResult.rayCount, sourceStrengthTotal);
        }

        /**
         * @brief Rebuild the resident cumulative source-strength array from excitation.
         * @param meshContainer Owning resident trace whose prefix buffer is updated.
         * @param betaVolume Current cell excitation view on the same device.
         * @return Total source strength copied back to the host.
         */
        double rebuildSourceStrengthPrefix(
            data::ResidentTrace<T_Device>& meshContainer,
            alpaka::concepts::IView<double> auto const& betaVolume)
        {
            auto mesh = meshContainer.view();
            mesh.betaVolume = std::span<double const>(betaVolume.data(), betaVolume.getExtents().x());
            if(mesh.numberOfCells == 0u)
            {
                alpaka::onHost::fill(m_queue, m_sourceStrengthTotal.toDeviceView(), 0.0, alpaka::Vec{std::size_t{1}});
            }
            else
            {
                auto sourceStrengthPrefix = meshContainer.sourceStrengthPrefix.toDeviceView();
                auto const cellFrameSpec = alpakaUtils::getFrameSpec<uint32_t>(
                    m_devBundle.device,
                    m_devBundle.executor,
                    alpaka::Vec{mesh.numberOfCells});
                m_queue.enqueue(
                    cellFrameSpec,
                    alpaka::KernelBundle{kernels::BuildSourceStrengthWeights{}, mesh, sourceStrengthPrefix});
                alpaka::onHost::inclusiveScanInPlace(
                    m_queue,
                    m_devBundle.executor,
                    m_sourceStrengthPrefixScanBuffer,
                    sourceStrengthPrefix);
                auto const scalarFrameSpec
                    = alpaka::onHost::getFrameSpec(m_devBundle.device, m_devBundle.executor, alpaka::Vec{1u});
                m_queue.enqueue(
                    scalarFrameSpec,
                    alpaka::KernelBundle{
                        hase::kernels::CaptureSourceStrengthTotal{},
                        mesh.numberOfCells,
                        sourceStrengthPrefix,
                        m_sourceStrengthTotal.toDeviceView()});
            }
            m_sourceStrengthTotal.toHost(m_queue);
            return m_sourceStrengthTotal.getHostView()[0u];
        }

        /** @return Cell-ordered ASE population derivative after waiting for its download. */
        std::vector<double> downloadVolumeDndtAse()
        {
            std::vector<double> result(m_volumeCount, 0.0);
            alpaka::onHost::memcpy(m_queue, result, m_volumeDndtAse);
            alpaka::onHost::wait(m_queue);
            return result;
        }

        /**
         * @param includePhiAse Whether to download ASE flux.
         * @param includeStandardError Whether to download absolute standard error.
         * @param includeRelativeStandardError Whether to download relative standard error.
         * @param includeTotalRays Whether to download visit and dropped-ray counters.
         * @return Finalized result containing only requested arrays.
         */
        data::PhiAseResult downloadFinalizedResult(
            bool includePhiAse,
            bool includeStandardError,
            bool includeRelativeStandardError,
            bool includeTotalRays)
        {
            data::PhiAseResult result;
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

        /** @return Owning device buffer containing finalized cell ASE flux. */
        [[nodiscard]] auto& volumePhiAse()
        {
            return m_volumePhiAse;
        }

        /** @return Owning device buffer containing the cell ASE population derivative. */
        [[nodiscard]] auto& volumeDndtAse()
        {
            return m_volumeDndtAse;
        }

    private:
        using ExactSurfaceReservoirScratch
            = SurfaceReservoirScratch<T_Device, hase::kernels::forward::surfaceReservoirPosition::Exact>;
        using CentroidSurfaceReservoirScratch
            = SurfaceReservoirScratch<T_Device, kernels::forward::surfaceReservoirPosition::Centroid>;

        alpakaUtils::DevBundle<T_Device, T_Exec> m_devBundle;
        T_Queue m_queue;
        std::vector<std::uint32_t> m_rseBatchRayCounts;
        T_RseBatchRayCountsBuffer m_rseBatchRayCountsBuffer;
        std::array<double, 1u> m_sourceStrengthTotalHost;
        T_BetaVolumeTotalBuffer m_sourceStrengthTotal;
        T_DoubleBuffer m_vertexBatchScoreSum;
        T_UnsignedBuffer m_volumeRayVisits;
        T_UnsignedBuffer m_droppedRays;
        T_FloatBuffer m_volumePhiAse;
        T_DoubleBuffer m_standardError;
        T_DoubleBuffer m_relativeStandardError;
        T_DoubleBuffer m_volumeDndtAse;
        T_ByteBuffer m_sourceStrengthPrefixScanBuffer;
        std::unique_ptr<ReflectionResamplingScratch<T_Device>> m_directReflectionScratch;
        std::unique_ptr<ExactSurfaceReservoirScratch> m_exactSurfaceReservoirScratch;
        std::unique_ptr<CentroidSurfaceReservoirScratch> m_centroidSurfaceReservoirScratch;
        ForwardPhiAseRawResult m_srmResult;
        unsigned m_volumeCount;
        unsigned m_materialVertexCount;
        unsigned m_batchCount;
        unsigned m_rayCount = 0u;
        unsigned m_accumulatedRayCount = 0u;
        std::chrono::steady_clock::time_point m_started;
    };

} // namespace hase::core
