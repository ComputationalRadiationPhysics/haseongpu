/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/calcForwardPhiAse.hpp>
#include <core/cancellableBarrier.hpp>
#include <core/haseWorker.hpp>

#include <algorithm>
#include <any>
#include <chrono>
#include <concepts>
#include <limits>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace hase::core
{
    struct ForwardPopulationOperation
    {
    };

    struct PrepareRayPopulationWork : ForwardPopulationOperation
    {
        std::span<ForwardPopulationBatch const> work;
        std::uint32_t seed;

        bool execute(
            alpaka::concepts::SpecializationOf<ForwardPhiAseDeviceContext> auto& context,
            data::TraceView mesh,
            data::AseDomainSourceView sources,
            data::AseDomainInterfaceView,
            AseTraceControls const&) const
        {
            context.prepareRayPopulations(mesh, sources, work, seed);
            return true;
        }
    };

    struct TraceRayPopulationWork : ForwardPopulationOperation
    {
        std::uint32_t rayPopulationId;
        std::span<ForwardPopulationRay const> incoming;
        bool primary;

        auto execute(
            alpaka::concepts::SpecializationOf<ForwardPhiAseDeviceContext> auto& context,
            data::TraceView mesh,
            data::AseDomainSourceView,
            data::AseDomainInterfaceView interfaces,
            AseTraceControls const& controls) const
        {
            return context.traceRayPopulation(mesh, interfaces, controls, rayPopulationId, incoming, primary);
        }
    };

    struct TraceResidentRayPopulationWork : ForwardPopulationOperation
    {
        std::uint32_t rayPopulationId;
        bool primary;

        double execute(
            alpaka::concepts::SpecializationOf<ForwardPhiAseDeviceContext> auto& context,
            data::TraceView mesh,
            data::AseDomainSourceView,
            data::AseDomainInterfaceView interfaces,
            AseTraceControls const& controls) const
        {
            return context.traceResidentRayPopulation(mesh, interfaces, controls, rayPopulationId, primary);
        }
    };

    struct SelectResidentRayPopulationWork : ForwardPopulationOperation
    {
        std::uint32_t seed;
        std::uint32_t pass;

        void execute(
            alpaka::concepts::SpecializationOf<ForwardPhiAseDeviceContext> auto& context,
            data::TraceView,
            data::AseDomainSourceView,
            data::AseDomainInterfaceView,
            AseTraceControls const& controls) const
        {
            context.selectResidentRayPopulation(forwardPopulationDomainCount(controls.domainCount), seed, pass);
        }
    };

    struct SelectRayPopulationWork : ForwardPopulationOperation
    {
        std::vector<ForwardPopulationRay> candidates;
        std::uint32_t seed;
        std::uint32_t pass;

        auto execute(
            alpaka::concepts::SpecializationOf<ForwardPhiAseDeviceContext> auto& context,
            data::TraceView,
            data::AseDomainSourceView,
            data::AseDomainInterfaceView,
            AseTraceControls const& controls)
        {
            return context.selectRayPopulation(
                std::move(candidates),
                forwardPopulationDomainCount(controls.domainCount),
                seed,
                pass);
        }
    };

    struct CollectRayPopulationWork : ForwardPopulationOperation
    {
        auto execute(
            alpaka::concepts::SpecializationOf<ForwardPhiAseDeviceContext> auto& context,
            data::TraceView,
            data::AseDomainSourceView,
            data::AseDomainInterfaceView,
            AseTraceControls const&) const
        {
            ForwardPhiAseRawResult raw;
            float runtime = 0.0f;
            context.finish(raw, runtime);
            return raw;
        }
    };

    struct TraceLogicalSrmBatches : ForwardPopulationOperation
    {
        std::uint32_t seed;

        bool execute(
            alpaka::concepts::SpecializationOf<ForwardPhiAseDeviceContext> auto& context,
            data::TraceView mesh,
            data::AseDomainSourceView,
            data::AseDomainInterfaceView interfaces,
            AseTraceControls const& controls) const
        {
            context.traceLogicalSrmBatches(mesh, interfaces, controls, seed);
            return true;
        }
    };

    /** @brief One statistical batch containing this worker's assigned domain work. */
    struct ForwardRayPopulation
    {
        unsigned index = 0u; //!< Statistical batch index.
        unsigned rayCount = 0u; //!< Complete number of histories in this batch.
        unsigned rngSeed = 0u; //!< Seed shared by all rayPopulations in one adaptive launch.
        std::vector<std::uint32_t> domainRayCounts;
        std::vector<double> domainSourceWeights;
        std::vector<std::uint32_t> domainPopulationCounts;
    };

    /** @brief Domain/batch assignments grouped for one worker and adaptive launch. */
    struct ForwardRayPopulationGroup
    {
        std::vector<ForwardRayPopulation> rayPopulations;
    };

    /** @brief One worker's device-resident accumulation downloaded once per adaptive launch. */
    struct ForwardWorkerResult
    {
        ForwardPhiAseRawResult raw;
        float runtime = 0.0f;
    };

    /** @brief Device-finalization work item created after batch gathering. */
    struct FinalizeForwardAse
    {
        ForwardPhiAseRawResult const& raw; //!< Gathered, unnormalized batch accumulators.
        bool upload = true; //!< Upload `raw` before finalization.
        bool downloadFullResult = true; //!< Download all public arrays rather than convergence only.
    };

    namespace detail
    {
        /** @brief Reusable host-thread collective storage for a single worker group. */
        class ThreadWorkerGroup
        {
        public:
            explicit ThreadWorkerGroup(unsigned const workerCount) : m_barrier(workerCount), m_values(workerCount)
            {
                if(workerCount == 0u)
                    throw std::invalid_argument("a HASE thread worker group cannot be empty");
            }

            void cancel(std::exception_ptr failure)
            {
                m_barrier.cancel(std::move(failure));
            }

            template<typename T_Value>
            [[nodiscard]] std::shared_ptr<std::vector<T_Value> const> gather(unsigned const workerIndex, T_Value value)
            {
                m_values.at(workerIndex) = std::move(value);
                m_barrier.arriveAndWait();
                if(workerIndex == 0u)
                {
                    auto gathered = std::make_shared<std::vector<T_Value>>();
                    gathered->reserve(m_values.size());
                    for(auto& item : m_values)
                        gathered->emplace_back(std::any_cast<T_Value>(std::move(item)));
                    m_collectiveValue = std::move(gathered);
                }
                m_barrier.arriveAndWait();
                auto result = std::any_cast<std::shared_ptr<std::vector<T_Value>>>(m_collectiveValue);
                m_barrier.arriveAndWait();
                return result;
            }

            template<typename T_Value>
            [[nodiscard]] T_Value scatter(unsigned const workerIndex, T_Value value)
            {
                if(workerIndex == 0u)
                    m_collectiveValue = std::move(value);
                m_barrier.arriveAndWait();
                auto result = std::any_cast<T_Value>(m_collectiveValue);
                m_barrier.arriveAndWait();
                return result;
            }

            template<typename T_Value>
            [[nodiscard]] T_Value reduce(
                unsigned const workerIndex,
                T_Value value,
                std::invocable<T_Value, T_Value const&> auto reduction)
            {
                auto values = gather(workerIndex, std::move(value));
                T_Value result{};
                if(workerIndex == 0u)
                {
                    result = values->front();
                    for(auto const& item : *values | std::views::drop(1u))
                        result = reduction(std::move(result), item);
                }
                return scatter(workerIndex, std::move(result));
            }

        private:
            CancellableBarrier m_barrier;
            std::vector<std::any> m_values;
            std::any m_collectiveValue;
        };
    } // namespace detail

    /**
     * @brief Worker policy representing one host thread that owns exactly one device.
     *
     * A worker executes zero or more scheduler-owned domain/batch assignments.
     * It never repartitions an assignment and never shares its device queue,
     * accumulators, or boundary buffers with another worker.
     *
     * @tparam T_Device Alpaka device type.
     * @tparam T_Exec Alpaka executor type.
     */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    class ThreadOwnedDevices
    {
    public:
        /** @brief Bind one worker identity to one persistent device context. */
        ThreadOwnedDevices(
            unsigned const workerIndex,
            unsigned const workerCount,
            detail::ThreadWorkerGroup& group,
            hase::data::TraceView const mesh,
            ForwardPhiAseDeviceContext<T_Device, T_Exec>& deviceContext,
            AseTraceControls const& experiment,
            double const betaVolumeTotal,
            hase::data::AseDomainInterfaceView const interfaceMap,
            hase::data::AseDomainSourceView const domainSources)
            : m_workerIndex(workerIndex)
            , m_workerCount(workerCount)
            , m_group(group)
            , m_mesh(mesh)
            , m_deviceContext(deviceContext)
            , m_experiment(experiment)
            , m_betaVolumeTotal(betaVolumeTotal)
            , m_interfaceMap(interfaceMap)
            , m_domainSources(domainSources)
        {
            if(workerIndex >= workerCount)
                throw std::out_of_range("thread worker index exceeds worker count");
        }

    private:
        unsigned m_workerIndex;
        unsigned m_workerCount;
        detail::ThreadWorkerGroup& m_group;
        hase::data::TraceView m_mesh;
        ForwardPhiAseDeviceContext<T_Device, T_Exec>& m_deviceContext;
        AseTraceControls const& m_experiment;
        double m_betaVolumeTotal;
        hase::data::AseDomainInterfaceView m_interfaceMap;
        hase::data::AseDomainSourceView m_domainSources;

        friend struct HaseWorkerDispatch<ThreadOwnedDevices<T_Device, T_Exec>>;
        friend struct HaseWorkItemDispatch<ThreadOwnedDevices<T_Device, T_Exec>, ForwardRayPopulationGroup>;
        friend struct HaseWorkItemDispatch<ThreadOwnedDevices<T_Device, T_Exec>, FinalizeForwardAse>;
        template<typename, typename>
        friend struct HaseWorkItemDispatch;
    };

    /** @brief Identity and collective dispatch for one-thread/one-device workers. */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    struct HaseWorkerDispatch<ThreadOwnedDevices<T_Device, T_Exec>>
    {
        using T_Policy = ThreadOwnedDevices<T_Device, T_Exec>;

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
                NodeId{0u},
                static_cast<std::uint32_t>(policy.m_workerIndex),
                1.0,
                std::numeric_limits<std::uint64_t>::max()};
        }

        [[nodiscard]] static bool requiresFinalizedDeviceState(T_Policy const& policy)
        {
            return isRoot(policy);
        }

        template<typename T_Value>
        [[nodiscard]] static auto scatter(T_Policy& policy, T_Value&& value)
        {
            using T = std::remove_cvref_t<T_Value>;
            return policy.m_group.scatter(policy.m_workerIndex, T(std::forward<T_Value>(value)));
        }

        template<typename T_Value>
        [[nodiscard]] static std::vector<T_Value> distribute(T_Policy& policy, std::vector<std::vector<T_Value>> value)
        {
            auto rootParts = policy.m_workerIndex == 0u
                                 ? std::make_shared<std::vector<std::vector<T_Value>>>(std::move(value))
                                 : std::shared_ptr<std::vector<std::vector<T_Value>>>{};
            auto parts = policy.m_group.scatter(policy.m_workerIndex, std::move(rootParts));
            return std::move(parts->at(policy.m_workerIndex));
        }

        template<typename T_Value>
        [[nodiscard]] static auto gather(T_Policy& policy, T_Value&& value)
        {
            using T = std::remove_cvref_t<T_Value>;
            return policy.m_group.gather(policy.m_workerIndex, T(std::forward<T_Value>(value)));
        }

        template<typename T_Value>
        [[nodiscard]] static auto reduce(
            T_Policy& policy,
            T_Value&& value,
            std::invocable<std::remove_cvref_t<T_Value>, std::remove_cvref_t<T_Value> const&> auto reduction)
        {
            using T = std::remove_cvref_t<T_Value>;
            return policy.m_group.reduce(policy.m_workerIndex, T(std::forward<T_Value>(value)), std::move(reduction));
        }
    };

    template<
        alpaka::onHost::concepts::Device T_Device,
        alpaka::concepts::Executor T_Exec,
        std::derived_from<ForwardPopulationOperation> T_Work>
    struct HaseWorkItemDispatch<ThreadOwnedDevices<T_Device, T_Exec>, T_Work>
    {
        static auto run(ThreadOwnedDevices<T_Device, T_Exec>& policy, T_Work work)
        {
            return work.execute(
                policy.m_deviceContext,
                policy.m_mesh,
                policy.m_domainSources,
                policy.m_interfaceMap,
                policy.m_experiment);
        }
    };

    /** @brief Enqueue all locally assigned rayPopulations and download their shared accumulator once. */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    struct HaseWorkItemDispatch<ThreadOwnedDevices<T_Device, T_Exec>, ForwardRayPopulationGroup>
    {
        using T_Policy = ThreadOwnedDevices<T_Device, T_Exec>;

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

    /** @brief Finalize gathered rayPopulations on one thread-owned device. */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    struct HaseWorkItemDispatch<ThreadOwnedDevices<T_Device, T_Exec>, FinalizeForwardAse>
    {
        using T_Policy = ThreadOwnedDevices<T_Device, T_Exec>;

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
} // namespace hase::core
