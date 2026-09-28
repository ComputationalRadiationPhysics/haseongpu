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

#    include <core/Runtime.hpp>
#    include <core/calcPhiAseThreaded.hpp>
#    include <mpi.h>

#    include <algorithm>
#    include <array>
#    include <concepts>
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
     * Each rank executes zero or more scheduler-owned domain/batch assignments.
     * Device state, asynchronous queues, accumulators, and boundary buffers
     * remain rank-private. Only batch results and generic
     * collective values cross rank boundaries.
     *
     * @tparam T_Device Alpaka device type.
     * @tparam T_Exec Alpaka executor type.
     */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    class MPIRank
    {
    public:
        /** @brief Bind the current rank to exactly one persistent local device context. */
        MPIRank(
            MPI_Comm const communicator,
            hase::data::TraceView const mesh,
            ForwardPhiAseDeviceContext<T_Device, T_Exec>& deviceContext,
            AseTraceControls const& experiment,
            double const betaVolumeTotal,
            unsigned const volumeCount,
            unsigned const vertexCount,
            unsigned const numIndependentRayPopulations,
            hase::data::AseDomainInterfaceView const interfaceMap,
            hase::data::AseDomainSourceView const domainSources)
            : m_communicator(communicator)
            , m_mesh(mesh)
            , m_deviceContext(deviceContext)
            , m_experiment(experiment)
            , m_betaVolumeTotal(betaVolumeTotal)
            , m_volumeCount(volumeCount)
            , m_vertexCount(vertexCount)
            , m_numIndependentRayPopulations(numIndependentRayPopulations)
            , m_interfaceMap(interfaceMap)
            , m_domainSources(domainSources)
        {
            int rank = 0;
            int size = 0;
            MPI_Comm_rank(m_communicator, &rank);
            MPI_Comm_size(m_communicator, &size);
            m_workerIndex = static_cast<unsigned>(rank);
            m_workerCount = static_cast<unsigned>(size);
            MPI_Comm nodeCommunicator = MPI_COMM_NULL;
            MPI_Comm_split_type(m_communicator, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &nodeCommunicator);
            int localRank = 0;
            MPI_Comm_rank(nodeCommunicator, &localRank);
            int nodeLeader = rank;
            MPI_Allreduce(&rank, &nodeLeader, 1, MPI_INT, MPI_MIN, nodeCommunicator);
            MPI_Comm_free(&nodeCommunicator);
            m_node = static_cast<NodeId>(nodeLeader);
            m_localWorker = static_cast<std::uint32_t>(localRank);
        }

    private:
        MPI_Comm m_communicator;
        unsigned m_workerIndex = 0u;
        unsigned m_workerCount = 1u;
        NodeId m_node = 0u;
        std::uint32_t m_localWorker = 0u;
        hase::data::TraceView m_mesh;
        ForwardPhiAseDeviceContext<T_Device, T_Exec>& m_deviceContext;
        AseTraceControls const& m_experiment;
        double m_betaVolumeTotal;
        unsigned m_volumeCount;
        unsigned m_vertexCount;
        unsigned m_numIndependentRayPopulations;
        hase::data::AseDomainInterfaceView m_interfaceMap;
        hase::data::AseDomainSourceView m_domainSources;

        friend struct HaseWorkerDispatch<MPIRank<T_Device, T_Exec>>;
        friend struct HaseWorkItemDispatch<MPIRank<T_Device, T_Exec>, ForwardRayPopulationGroup>;
        friend struct HaseWorkItemDispatch<MPIRank<T_Device, T_Exec>, FinalizeForwardAse>;
        template<typename, typename>
        friend struct HaseWorkItemDispatch;
    };

    /** @brief Identity and collective dispatch for one-rank/one-device workers. */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
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

        [[nodiscard]] static WorkerDescriptor descriptor(T_Policy const& policy)
        {
            return WorkerDescriptor{
                static_cast<WorkerId>(policy.m_workerIndex),
                policy.m_node,
                policy.m_localWorker,
                1.0,
                std::numeric_limits<std::uint64_t>::max()};
        }

        [[nodiscard]] static bool requiresFinalizedDeviceState(T_Policy const&)
        {
            return true;
        }

        template<typename T_Value>
        [[nodiscard]] static T_Value scatter(T_Policy& policy, T_Value value)
        {
            static_assert(std::is_trivially_copyable_v<T_Value>, "MPI scatter values must be trivially copyable");
            MPI_Bcast(std::addressof(value), static_cast<int>(sizeof(T_Value)), MPI_BYTE, 0, policy.m_communicator);
            return value;
        }

        template<typename T_Value>
        requires std::is_trivially_copyable_v<T_Value>
        [[nodiscard]] static std::vector<T_Value> scatter(T_Policy& policy, std::vector<T_Value> value)
        {
            std::uint64_t count = value.size();
            MPI_Bcast(&count, 1, MPI_UINT64_T, 0, policy.m_communicator);
            if(count > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) / sizeof(T_Value))
                throw std::overflow_error("MPI population broadcast exceeds the byte-count range");
            value.resize(static_cast<std::size_t>(count));
            MPI_Bcast(value.data(), static_cast<int>(count * sizeof(T_Value)), MPI_BYTE, 0, policy.m_communicator);
            return value;
        }

        template<typename T_Value>
        requires std::is_trivially_copyable_v<T_Value>
        [[nodiscard]] static std::vector<T_Value> distribute(T_Policy& policy, std::vector<std::vector<T_Value>> value)
        {
            std::vector<std::uint64_t> counts(policy.m_workerCount, 0u);
            int valid = 1;
            if(policy.m_workerIndex == 0u)
            {
                valid = value.size() == policy.m_workerCount;
                if(valid)
                    for(unsigned rank = 0u; rank < policy.m_workerCount; ++rank)
                    {
                        counts[rank] = value[rank].size();
                        valid &= counts[rank]
                                 <= static_cast<std::uint64_t>(std::numeric_limits<int>::max()) / sizeof(T_Value);
                    }
            }
            MPI_Bcast(&valid, 1, MPI_INT, 0, policy.m_communicator);
            if(!valid)
                throw std::invalid_argument("MPI population distribution has invalid worker counts");
            std::uint64_t localCount = 0u;
            MPI_Scatter(counts.data(), 1, MPI_UINT64_T, &localCount, 1, MPI_UINT64_T, 0, policy.m_communicator);
            if(policy.m_workerIndex == 0u)
            {
                for(unsigned rank = 1u; rank < policy.m_workerCount; ++rank)
                    MPI_Send(
                        value[rank].data(),
                        static_cast<int>(value[rank].size() * sizeof(T_Value)),
                        MPI_BYTE,
                        static_cast<int>(rank),
                        0,
                        policy.m_communicator);
                return std::move(value.front());
            }
            std::vector<T_Value> local(static_cast<std::size_t>(localCount));
            MPI_Recv(
                local.data(),
                static_cast<int>(localCount * sizeof(T_Value)),
                MPI_BYTE,
                0,
                0,
                policy.m_communicator,
                MPI_STATUS_IGNORE);
            return local;
        }

        template<typename T_Value>
        requires std::is_trivially_copyable_v<T_Value>
        [[nodiscard]] static std::shared_ptr<std::vector<std::vector<T_Value>> const> gather(
            T_Policy& policy,
            std::vector<T_Value> value)
        {
            auto result = std::make_shared<std::vector<std::vector<T_Value>>>(policy.m_workerCount);
            // Only the combing owner needs ray payloads. Every rank sees counts so the
            // collective order stays fixed, but payloads travel only to rank zero.
            for(unsigned rank = 0u; rank < policy.m_workerCount; ++rank)
            {
                auto& part = result->at(rank);
                if(rank == policy.m_workerIndex)
                    part = std::move(value);
                std::uint64_t count = part.size();
                MPI_Bcast(&count, 1, MPI_UINT64_T, static_cast<int>(rank), policy.m_communicator);
                if(count > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) / sizeof(T_Value))
                    throw std::overflow_error("MPI population gather exceeds the byte-count range");
                if(rank != 0u)
                {
                    if(policy.m_workerIndex == rank)
                        MPI_Send(
                            part.data(),
                            static_cast<int>(count * sizeof(T_Value)),
                            MPI_BYTE,
                            0,
                            0,
                            policy.m_communicator);
                    else if(policy.m_workerIndex == 0u)
                    {
                        part.resize(static_cast<std::size_t>(count));
                        MPI_Recv(
                            part.data(),
                            static_cast<int>(count * sizeof(T_Value)),
                            MPI_BYTE,
                            static_cast<int>(rank),
                            0,
                            policy.m_communicator,
                            MPI_STATUS_IGNORE);
                    }
                }
            }
            return result;
        }

        [[nodiscard]] static std::shared_ptr<std::vector<ForwardWorkerResult> const> gather(
            T_Policy& policy,
            ForwardWorkerResult local)
        {
            auto gathered = std::make_shared<std::vector<ForwardWorkerResult>>(1u);
            auto& global = gathered->front();
            global.raw = makeForwardRawResult(
                policy.m_volumeCount,
                policy.m_vertexCount,
                policy.m_numIndependentRayPopulations);
            MPI_Allreduce(
                local.raw.vertexPopulationScoreSum.data(),
                global.raw.vertexPopulationScoreSum.data(),
                static_cast<int>(global.raw.vertexPopulationScoreSum.size()),
                MPI_DOUBLE,
                MPI_SUM,
                policy.m_communicator);
            MPI_Allreduce(
                local.raw.rayPopulationRayCounts.data(),
                global.raw.rayPopulationRayCounts.data(),
                static_cast<int>(global.raw.rayPopulationRayCounts.size()),
                MPI_UNSIGNED,
                MPI_SUM,
                policy.m_communicator);
            MPI_Allreduce(
                local.raw.totalRays.data(),
                global.raw.totalRays.data(),
                static_cast<int>(global.raw.totalRays.size()),
                MPI_UNSIGNED,
                MPI_SUM,
                policy.m_communicator);
            MPI_Allreduce(
                local.raw.droppedRays.data(),
                global.raw.droppedRays.data(),
                static_cast<int>(global.raw.droppedRays.size()),
                MPI_UNSIGNED,
                MPI_SUM,
                policy.m_communicator);
            MPI_Allreduce(&local.raw.rayCount, &global.raw.rayCount, 1, MPI_UNSIGNED, MPI_SUM, policy.m_communicator);

            std::array localStatus{
                boundaryStatusPriority(local.raw.boundaryStatus),
                local.raw.boundaryPasses,
                local.raw.boundaryMaxPasses,
                local.raw.boundaryDivergenceStreak};
            std::array<unsigned, 4u> globalStatus{};
            MPI_Allreduce(
                localStatus.data(),
                globalStatus.data(),
                static_cast<int>(globalStatus.size()),
                MPI_UNSIGNED,
                MPI_MAX,
                policy.m_communicator);
            global.raw.boundaryStatus = boundaryStatusFromPriority(globalStatus[0u]);
            global.raw.boundaryPasses = globalStatus[1u];
            global.raw.boundaryMaxPasses = globalStatus[2u];
            global.raw.boundaryDivergenceStreak = globalStatus[3u];
            MPI_Allreduce(
                &local.raw.boundaryRemainingFraction,
                &global.raw.boundaryRemainingFraction,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                policy.m_communicator);
            std::array localTail{
                local.raw.boundaryGamma,
                local.raw.boundaryGammaStandardError,
                local.raw.boundaryTailFactor,
                local.raw.boundaryTailClosure};
            std::array<double, 4u> globalTail{};
            MPI_Allreduce(
                localTail.data(),
                globalTail.data(),
                static_cast<int>(globalTail.size()),
                MPI_DOUBLE,
                MPI_MAX,
                policy.m_communicator);
            global.raw.boundaryGamma = globalTail[0u];
            global.raw.boundaryGammaStandardError = globalTail[1u];
            global.raw.boundaryTailFactor = globalTail[2u];
            global.raw.boundaryTailClosure = globalTail[3u];
            return gathered;
        }

        template<typename T_Value>
        requires std::is_trivially_copyable_v<T_Value>
        [[nodiscard]] static std::shared_ptr<std::vector<T_Value> const> gather(T_Policy& policy, T_Value value)
        {
            auto gathered = std::make_shared<std::vector<T_Value>>(policy.m_workerCount);
            MPI_Allgather(
                std::addressof(value),
                static_cast<int>(sizeof(T_Value)),
                MPI_BYTE,
                gathered->data(),
                static_cast<int>(sizeof(T_Value)),
                MPI_BYTE,
                policy.m_communicator);
            return gathered;
        }

        template<typename T_Value>
        [[nodiscard]] static T_Value reduce(
            T_Policy& policy,
            T_Value value,
            std::invocable<T_Value, T_Value const&> auto reduction)
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

    template<
        alpaka::onHost::concepts::Device T_Device,
        alpaka::concepts::Executor T_Exec,
        std::derived_from<ForwardPopulationOperation> T_Work>
    struct HaseWorkItemDispatch<MPIRank<T_Device, T_Exec>, T_Work>
    {
        static auto run(MPIRank<T_Device, T_Exec>& policy, T_Work work)
        {
            return work.execute(
                policy.m_deviceContext,
                policy.m_mesh,
                policy.m_domainSources,
                policy.m_interfaceMap,
                policy.m_experiment);
        }
    };

    /** @brief Enqueue all rank-local rayPopulations and download their shared accumulator once. */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    struct HaseWorkItemDispatch<MPIRank<T_Device, T_Exec>, ForwardRayPopulationGroup>
    {
        using T_Policy = MPIRank<T_Device, T_Exec>;

        [[nodiscard]] static ForwardWorkerResult run(T_Policy& policy, ForwardRayPopulationGroup const& group)
        {
            ForwardWorkerResult result;
            if(group.rayPopulations.empty())
            {
                result.raw = policy.m_deviceContext.makeEmptyRawResult();
                return result;
            }

            auto const started = std::chrono::steady_clock::now();
            bool resetAccumulators = true;
            for(auto const& batch : group.rayPopulations)
            {
                policy.m_deviceContext.begin(
                    policy.m_mesh,
                    batch.rayCount,
                    batch.rngSeed,
                    batch.index,
                    policy.m_betaVolumeTotal,
                    policy.m_experiment,
                    policy.m_interfaceMap,
                    policy.m_domainSources,
                    batch.domainRayCounts,
                    batch.domainSourceWeights,
                    batch.domainPopulationCounts,
                    resetAccumulators);
                resetAccumulators = false;
            }
            float ignoredRuntime = 0.0f;
            policy.m_deviceContext.finish(result.raw, ignoredRuntime);
            result.runtime = static_cast<float>(
                std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count());
            return result;
        }
    };

    /** @brief Finalize gathered rayPopulations on the rank-owned device. */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    struct HaseWorkItemDispatch<MPIRank<T_Device, T_Exec>, FinalizeForwardAse>
    {
        using T_Policy = MPIRank<T_Device, T_Exec>;

        [[nodiscard]] static data::PhiAseResult run(T_Policy& policy, FinalizeForwardAse const& item)
        {
            if(item.upload)
                policy.m_deviceContext.uploadAndFinalize(policy.m_mesh, item.raw, policy.m_betaVolumeTotal);
            auto result = policy.m_deviceContext.downloadFinalizedResult(
                item.downloadFullResult,
                item.downloadFullResult,
                true,
                item.downloadFullResult);
            if(item.downloadFullResult)
                result.dndtAse = policy.m_deviceContext.downloadVolumeDndtAse();
            result.boundaryStatus = item.raw.boundaryStatus;
            result.boundaryPasses = item.raw.boundaryPasses;
            result.boundaryRemainingFraction = item.raw.boundaryRemainingFraction;
            result.boundaryMaxPasses = item.raw.boundaryMaxPasses;
            result.boundaryDivergenceStreak = item.raw.boundaryDivergenceStreak;
            result.boundaryGamma = item.raw.boundaryGamma;
            result.boundaryGammaStandardError = item.raw.boundaryGammaStandardError;
            result.boundaryTailFactor = item.raw.boundaryTailFactor;
            result.boundaryTailClosure = item.raw.boundaryTailClosure;
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
