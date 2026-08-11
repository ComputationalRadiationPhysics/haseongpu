/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#if !defined(DISABLE_MPI) && defined(MPI_FOUND)

#    include <core/calcPhiAseThreaded.hpp>
#    include <core/types.hpp>
#    include <mpi.h>

#    include <algorithm>
#    include <array>
#    include <limits>
#    include <memory>
#    include <stdexcept>
#    include <type_traits>
#    include <utility>
#    include <vector>

namespace hase::core
{
    namespace detail
    {
        /** @brief Process-lifetime MPI initialization owned by HASE when required. */
        class MpiLifetime
        {
        public:
            MpiLifetime()
            {
                int initialized = 0;
                MPI_Initialized(&initialized);
                m_owned = initialized == 0;
                if(m_owned)
                {
                    int provided = MPI_THREAD_SINGLE;
                    if(MPI_Init_thread(nullptr, nullptr, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS
                       || provided < MPI_THREAD_FUNNELED)
                        throw std::runtime_error("failed to initialize MPI with funneled thread support");
                }
            }

            ~MpiLifetime()
            {
                int finalized = 0;
                MPI_Finalized(&finalized);
                if(m_owned && finalized == 0)
                    MPI_Finalize();
            }

            MpiLifetime(MpiLifetime const&) = delete;
            MpiLifetime& operator=(MpiLifetime const&) = delete;

        private:
            bool m_owned = false;
        };

        inline void ensureMpiInitialized()
        {
            static MpiLifetime lifetime;
        }
    } // namespace detail

    /**
     * @brief Worker policy representing one MPI rank that owns one local device.
     *
     * Each rank executes zero or more complete batches selected by
     * `hase::mapIdx`. Device state, asynchronous queues, accumulators, and
     * reflection reservoirs remain rank-private. Only batch results and generic
     * collective values cross rank boundaries.
     *
     * @tparam T_Device Alpaka device type.
     * @tparam T_Exec Alpaka executor type.
     */
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    class MPIRank
    {
    public:
        /** @brief Bind the current rank to exactly one persistent local device context. */
        MPIRank(
            MPI_Comm const communicator,
            DeviceMeshView const mesh,
            ForwardPhiAseDeviceContext<T_Device, T_Exec>& deviceContext,
            ExperimentParameters const& experiment,
            double const betaVolumeTotal,
            unsigned const volumeCount,
            unsigned const vertexCount,
            unsigned const batchCount)
            : m_communicator(communicator)
            , m_mesh(mesh)
            , m_deviceContext(deviceContext)
            , m_experiment(experiment)
            , m_betaVolumeTotal(betaVolumeTotal)
            , m_volumeCount(volumeCount)
            , m_vertexCount(vertexCount)
            , m_batchCount(batchCount)
        {
            int rank = 0;
            int size = 0;
            MPI_Comm_rank(m_communicator, &rank);
            MPI_Comm_size(m_communicator, &size);
            m_workerIndex = static_cast<unsigned>(rank);
            m_workerCount = static_cast<unsigned>(size);
        }

    private:
        MPI_Comm m_communicator;
        unsigned m_workerIndex = 0u;
        unsigned m_workerCount = 1u;
        DeviceMeshView m_mesh;
        ForwardPhiAseDeviceContext<T_Device, T_Exec>& m_deviceContext;
        ExperimentParameters const& m_experiment;
        double m_betaVolumeTotal;
        unsigned m_volumeCount;
        unsigned m_vertexCount;
        unsigned m_batchCount;

        friend struct HaseWorkerDispatch<MPIRank<T_Device, T_Exec>>;
        friend struct HaseWorkItemDispatch<MPIRank<T_Device, T_Exec>, ForwardRayBatch>;
        friend struct HaseWorkItemDispatch<MPIRank<T_Device, T_Exec>, FinalizeForwardAse>;
    };

    /** @brief Identity and collective dispatch for one-rank/one-device workers. */
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    struct HaseWorkerDispatch<MPIRank<T_Device, T_Exec>>
    {
        using T_Policy = MPIRank<T_Device, T_Exec>;

        [[nodiscard]] static unsigned workerIndex(T_Policy const& policy)
        {
            return policy.m_workerIndex;
        }

        [[nodiscard]] static unsigned workerCount(T_Policy const& policy)
        {
            return policy.m_workerCount;
        }

        [[nodiscard]] static bool isRoot(T_Policy const& policy)
        {
            return policy.m_workerIndex == 0u;
        }

        template<typename T_Value>
        [[nodiscard]] static T_Value scatter(T_Policy& policy, T_Value value)
        {
            static_assert(std::is_trivially_copyable_v<T_Value>, "MPI scatter values must be trivially copyable");
            MPI_Bcast(std::addressof(value), static_cast<int>(sizeof(T_Value)), MPI_BYTE, 0, policy.m_communicator);
            return value;
        }

        [[nodiscard]] static std::shared_ptr<std::vector<ForwardRayBatchResults> const> gather(
            T_Policy& policy,
            ForwardRayBatchResults localResults)
        {
            auto gathered = std::make_shared<std::vector<ForwardRayBatchResults>>(1u);
            auto& globalResults = gathered->front();
            globalResults.reserve(policy.m_batchCount);
            for(unsigned batch = 0u; batch < policy.m_batchCount; ++batch)
            {
                auto const local = std::ranges::find(localResults, batch, &ForwardRayBatchResult::index);
                unsigned const localPresence = local == localResults.end() ? 0u : 1u;
                unsigned globalPresence = 0u;
                MPI_Allreduce(&localPresence, &globalPresence, 1, MPI_UNSIGNED, MPI_SUM, policy.m_communicator);
                if(globalPresence != 1u)
                    throw std::runtime_error("each forward statistical batch must be produced by exactly one rank");

                ForwardRayBatchResult global;
                global.index = batch;
                global.raw = makeForwardRawResult(policy.m_volumeCount, policy.m_vertexCount, policy.m_batchCount);
                ForwardPhiAseRawResult localRaw
                    = makeForwardRawResult(policy.m_volumeCount, policy.m_vertexCount, policy.m_batchCount);
                float localRuntime = 0.0f;
                if(localPresence != 0u)
                {
                    localRaw = local->raw;
                    localRuntime = local->runtime;
                }
                MPI_Allreduce(
                    localRaw.vertexBatchScoreSum.data(),
                    global.raw.vertexBatchScoreSum.data(),
                    static_cast<int>(global.raw.vertexBatchScoreSum.size()),
                    MPI_DOUBLE,
                    MPI_SUM,
                    policy.m_communicator);
                MPI_Allreduce(
                    localRaw.rseBatchRayCounts.data(),
                    global.raw.rseBatchRayCounts.data(),
                    static_cast<int>(global.raw.rseBatchRayCounts.size()),
                    MPI_UNSIGNED,
                    MPI_SUM,
                    policy.m_communicator);
                MPI_Allreduce(
                    localRaw.totalRays.data(),
                    global.raw.totalRays.data(),
                    static_cast<int>(global.raw.totalRays.size()),
                    MPI_UNSIGNED,
                    MPI_SUM,
                    policy.m_communicator);
                MPI_Allreduce(
                    localRaw.droppedRays.data(),
                    global.raw.droppedRays.data(),
                    static_cast<int>(global.raw.droppedRays.size()),
                    MPI_UNSIGNED,
                    MPI_SUM,
                    policy.m_communicator);
                MPI_Allreduce(
                    &localRaw.rayCount,
                    &global.raw.rayCount,
                    1,
                    MPI_UNSIGNED,
                    MPI_SUM,
                    policy.m_communicator);
                MPI_Allreduce(&localRuntime, &global.runtime, 1, MPI_FLOAT, MPI_SUM, policy.m_communicator);

                std::array<unsigned, 5u> localStatus{
                    localPresence == 0u ? 0u : static_cast<unsigned>(localRaw.srmStatus),
                    localRaw.srmPasses,
                    localRaw.srmMaxIterations,
                    localRaw.srmDivergenceStreak,
                    localPresence};
                std::array<unsigned, 5u> globalStatus{};
                MPI_Allreduce(
                    localStatus.data(),
                    globalStatus.data(),
                    static_cast<int>(globalStatus.size()),
                    MPI_UNSIGNED,
                    MPI_SUM,
                    policy.m_communicator);
                global.raw.srmStatus = static_cast<SrmStatus>(globalStatus[0u]);
                global.raw.srmPasses = globalStatus[1u];
                global.raw.srmMaxIterations = globalStatus[2u];
                global.raw.srmDivergenceStreak = globalStatus[3u];
                MPI_Allreduce(
                    &localRaw.srmRemainingFraction,
                    &global.raw.srmRemainingFraction,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    policy.m_communicator);
                globalResults.emplace_back(std::move(global));
            }
            return gathered;
        }

        template<typename T_Value, typename T_Reduction>
        [[nodiscard]] static T_Value reduce(T_Policy& policy, T_Value value, T_Reduction reduction)
        {
            static_assert(std::is_trivially_copyable_v<T_Value>, "MPI reduction values must be trivially copyable");
            std::vector<T_Value> values(policy.m_workerCount);
            MPI_Allgather(
                std::addressof(value),
                static_cast<int>(sizeof(T_Value)),
                MPI_BYTE,
                values.data(),
                static_cast<int>(sizeof(T_Value)),
                MPI_BYTE,
                policy.m_communicator);
            T_Value result = values.front();
            for(auto const& item : values | std::views::drop(1u))
                result = reduction(std::move(result), item);
            return result;
        }
    };

    /** @brief Execute one complete forward-ray batch on the rank-owned device. */
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    struct HaseWorkItemDispatch<MPIRank<T_Device, T_Exec>, ForwardRayBatch>
    {
        using T_Policy = MPIRank<T_Device, T_Exec>;

        [[nodiscard]] static ForwardRayBatchResult run(T_Policy& policy, ForwardRayBatch const& batch)
        {
            ForwardRayBatchResult result;
            result.index = batch.index;
            policy.m_deviceContext.evaluate(
                policy.m_mesh,
                result.raw,
                result.runtime,
                batch.rayCount,
                batch.rngSeed,
                batch.index,
                policy.m_betaVolumeTotal,
                policy.m_experiment);
            return result;
        }
    };

    /** @brief Finalize gathered batches on the rank-owned device. */
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    struct HaseWorkItemDispatch<MPIRank<T_Device, T_Exec>, FinalizeForwardAse>
    {
        using T_Policy = MPIRank<T_Device, T_Exec>;

        [[nodiscard]] static Result run(T_Policy& policy, FinalizeForwardAse const& item)
        {
            policy.m_deviceContext.uploadAndFinalize(
                policy.m_mesh,
                item.raw,
                policy.m_betaVolumeTotal,
                item.fluorescenceRate,
                item.sigmaA,
                item.sigmaE);
            auto result = policy.m_deviceContext.downloadFinalizedResult(true, true, true, true);
            result.dndtAse = policy.m_deviceContext.downloadVolumeDndtAse();
            result.srmStatus = item.raw.srmStatus;
            result.srmPasses = item.raw.srmPasses;
            result.srmRemainingFraction = item.raw.srmRemainingFraction;
            result.srmMaxIterations = item.raw.srmMaxIterations;
            result.srmDivergenceStreak = item.raw.srmDivergenceStreak;
            return result;
        }
    };

    /** @brief Return the process-local device index assigned to the current MPI rank. */
    inline unsigned mpiRankDeviceIndex(unsigned const localDeviceCount)
    {
        if(localDeviceCount == 0u)
            throw std::runtime_error("MPI forward ASE requires at least one local device");
        MPI_Comm nodeCommunicator = MPI_COMM_NULL;
        int worldRank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, worldRank, MPI_INFO_NULL, &nodeCommunicator);
        int localRank = 0;
        MPI_Comm_rank(nodeCommunicator, &localRank);
        MPI_Comm_free(&nodeCommunicator);
        return static_cast<unsigned>(localRank) % localDeviceCount;
    }

    /** @brief Describe the active one-rank/one-device MPI worker topology. */
    inline RuntimeTopology mpiWorkerTopology()
    {
        int worldSize = 1;
        int worldRank = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
        MPI_Comm nodeCommunicator = MPI_COMM_NULL;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, worldRank, MPI_INFO_NULL, &nodeCommunicator);
        int localRank = 0;
        int localSize = 1;
        MPI_Comm_rank(nodeCommunicator, &localRank);
        MPI_Comm_size(nodeCommunicator, &localSize);
        MPI_Comm_free(&nodeCommunicator);

        int const nodeContribution = localRank == 0 ? 1 : 0;
        int nodeCount = 1;
        int minRanksPerNode = std::numeric_limits<int>::max();
        int maxRanksPerNode = 0;
        int const minContribution = localRank == 0 ? localSize : std::numeric_limits<int>::max();
        int const maxContribution = localRank == 0 ? localSize : 0;
        MPI_Allreduce(&nodeContribution, &nodeCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&minContribution, &minRanksPerNode, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&maxContribution, &maxRanksPerNode, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        RuntimeTopology topology;
        topology.activeNodes = static_cast<unsigned>(nodeCount);
        topology.activeRanks = static_cast<unsigned>(worldSize);
        topology.avgActiveRanksPerNode = static_cast<double>(worldSize) / static_cast<double>(nodeCount);
        topology.minActiveRanksPerNode = static_cast<unsigned>(minRanksPerNode);
        topology.maxActiveRanksPerNode = static_cast<unsigned>(maxRanksPerNode);
        topology.activeGpus = static_cast<unsigned>(worldSize);
        topology.avgGpusPerRank = 1.0;
        topology.avgGpusPerNode = static_cast<double>(worldSize) / static_cast<double>(nodeCount);
        topology.minGpusPerNode = static_cast<unsigned>(minRanksPerNode);
        topology.maxGpusPerNode = static_cast<unsigned>(maxRanksPerNode);
        return topology;
    }
} // namespace hase::core

#endif
