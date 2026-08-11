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
#include <core/haseWorker.hpp>

#include <algorithm>
#include <any>
#include <barrier>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace hase::core
{
    /** @brief One complete, indivisible forward-ray batch. */
    struct ForwardRayBatch
    {
        unsigned index = 0u; //!< Statistical batch index.
        unsigned rayCount = 0u; //!< Complete number of histories in this batch.
        unsigned rngSeed = 0u; //!< Seed shared by all batches in one adaptive launch.
    };

    /** @brief Result of one complete forward-ray batch, retaining its statistical identity. */
    struct ForwardRayBatchResult
    {
        unsigned index = 0u; //!< Statistical batch index.
        ForwardPhiAseRawResult raw; //!< Unnormalized batch accumulation.
        float runtime = 0.0f; //!< Device execution time in seconds.
    };

    /** @brief Batch-preserving result collection used by worker gather operations. */
    using ForwardRayBatchResults = std::vector<ForwardRayBatchResult>;

    /** @brief Device-finalization work item created after batch gathering. */
    struct FinalizeForwardAse
    {
        ForwardPhiAseRawResult const& raw; //!< Gathered, unnormalized batch accumulators.
        double fluorescenceRate; //!< Fluorescence normalization applied on the device.
        double sigmaA; //!< Absorption cross section used for ASE depletion.
        double sigmaE; //!< Emission cross section used for ASE depletion.
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

            template<typename T_Value>
            [[nodiscard]] std::shared_ptr<std::vector<T_Value> const> gather(unsigned const workerIndex, T_Value value)
            {
                m_values.at(workerIndex) = std::move(value);
                m_barrier.arrive_and_wait();
                if(workerIndex == 0u)
                {
                    auto gathered = std::make_shared<std::vector<T_Value>>();
                    gathered->reserve(m_values.size());
                    for(auto& item : m_values)
                        gathered->emplace_back(std::any_cast<T_Value>(std::move(item)));
                    m_collectiveValue = std::move(gathered);
                }
                m_barrier.arrive_and_wait();
                auto result = std::any_cast<std::shared_ptr<std::vector<T_Value>>>(m_collectiveValue);
                m_barrier.arrive_and_wait();
                return result;
            }

            template<typename T_Value>
            [[nodiscard]] T_Value scatter(unsigned const workerIndex, T_Value value)
            {
                if(workerIndex == 0u)
                    m_collectiveValue = std::move(value);
                m_barrier.arrive_and_wait();
                auto result = std::any_cast<T_Value>(m_collectiveValue);
                m_barrier.arrive_and_wait();
                return result;
            }

            template<typename T_Value, typename T_Reduction>
            [[nodiscard]] T_Value reduce(unsigned const workerIndex, T_Value value, T_Reduction reduction)
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
            std::barrier<> m_barrier;
            std::vector<std::any> m_values;
            std::any m_collectiveValue;
        };
    } // namespace detail

    /**
     * @brief Worker policy representing one host thread that owns exactly one device.
     *
     * A worker executes zero or more complete batches selected by `hase::mapIdx`.
     * It never repartitions a selected batch and never shares its device queue,
     * accumulators, or reflection reservoir with another worker.
     *
     * @tparam T_Device Alpaka device type.
     * @tparam T_Exec Alpaka executor type.
     */
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    class ThreadOwnedDevices
    {
    public:
        /** @brief Bind one worker identity to one persistent device context. */
        ThreadOwnedDevices(
            unsigned const workerIndex,
            unsigned const workerCount,
            detail::ThreadWorkerGroup& group,
            DeviceMeshView const mesh,
            ForwardPhiAseDeviceContext<T_Device, T_Exec>& deviceContext,
            ExperimentParameters const& experiment,
            double const betaVolumeTotal)
            : m_workerIndex(workerIndex)
            , m_workerCount(workerCount)
            , m_group(group)
            , m_mesh(mesh)
            , m_deviceContext(deviceContext)
            , m_experiment(experiment)
            , m_betaVolumeTotal(betaVolumeTotal)
        {
            if(workerIndex >= workerCount)
                throw std::out_of_range("thread worker index exceeds worker count");
        }

    private:
        unsigned m_workerIndex;
        unsigned m_workerCount;
        detail::ThreadWorkerGroup& m_group;
        DeviceMeshView m_mesh;
        ForwardPhiAseDeviceContext<T_Device, T_Exec>& m_deviceContext;
        ExperimentParameters const& m_experiment;
        double m_betaVolumeTotal;

        friend struct HaseWorkerDispatch<ThreadOwnedDevices<T_Device, T_Exec>>;
        friend struct HaseWorkItemDispatch<ThreadOwnedDevices<T_Device, T_Exec>, ForwardRayBatch>;
        friend struct HaseWorkItemDispatch<ThreadOwnedDevices<T_Device, T_Exec>, FinalizeForwardAse>;
    };

    /** @brief Identity and collective dispatch for one-thread/one-device workers. */
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
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

        template<typename T_Value>
        [[nodiscard]] static auto scatter(T_Policy& policy, T_Value&& value)
        {
            using T = std::remove_cvref_t<T_Value>;
            return policy.m_group.scatter(policy.m_workerIndex, T(std::forward<T_Value>(value)));
        }

        template<typename T_Value>
        [[nodiscard]] static auto gather(T_Policy& policy, T_Value&& value)
        {
            using T = std::remove_cvref_t<T_Value>;
            return policy.m_group.gather(policy.m_workerIndex, T(std::forward<T_Value>(value)));
        }

        template<typename T_Value, typename T_Reduction>
        [[nodiscard]] static auto reduce(T_Policy& policy, T_Value&& value, T_Reduction reduction)
        {
            using T = std::remove_cvref_t<T_Value>;
            return policy.m_group.reduce(policy.m_workerIndex, T(std::forward<T_Value>(value)), std::move(reduction));
        }
    };

    /** @brief Execute one complete forward-ray batch on one thread-owned device. */
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    struct HaseWorkItemDispatch<ThreadOwnedDevices<T_Device, T_Exec>, ForwardRayBatch>
    {
        using T_Policy = ThreadOwnedDevices<T_Device, T_Exec>;

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

    /** @brief Finalize gathered batches on one thread-owned device. */
    template<alpaka::onHost::concepts::Device T_Device, typename T_Exec>
    struct HaseWorkItemDispatch<ThreadOwnedDevices<T_Device, T_Exec>, FinalizeForwardAse>
    {
        using T_Policy = ThreadOwnedDevices<T_Device, T_Exec>;

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
} // namespace hase::core
