/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace hase
{
    /**
     * @brief Map a one-dimensional Alpaka index range onto one worker.
     *
     * Worker `i` of `n` receives the first range element at offset `i`, then
     * follows the original stride multiplied by `n`. Across all workers every
     * index is produced exactly once. A mapped work item is indivisible and
     * must not be repartitioned by the receiving worker.
     */
    [[nodiscard]] constexpr auto mapIdx(auto const& worker, auto const& range)
    {
        static_assert(std::remove_cvref_t<decltype(range)>::dim() == 1u, "HASE workers currently map 1D ranges");
        unsigned const workerCount = worker.workerCount();
        unsigned const workerIndex = worker.workerIndex();
        if(workerCount == 0u || workerIndex >= workerCount)
            throw std::invalid_argument("invalid HASE worker group identity");
        auto first = range.m_begin + range.m_stride * workerIndex;
        first[0u] = std::min(first[0u], range.m_end[0u]);
        return alpaka::IdxRange{first, range.m_end, range.m_stride * workerCount};
    }
} // namespace hase

namespace hase::core
{
    /**
     * @brief Policy customization point for worker identity and collectives.
     *
     * A worker policy specialization supplies `workerIndex`, `workerCount`,
     * `isRoot`, `scatter`, `gather`, and `reduce`. Collective operations must be
     * called in identical order by all workers in the group and outside mapped
     * loops, because workers may own different numbers of work items.
     */
    template<typename T_WorkerPolicy>
    struct HaseWorkerDispatch;

    /**
     * @brief Customization point for executing one indivisible work item.
     *
     * Physics features specialize this template for their work-item type.
     * Adding pump or ASE batches therefore does not extend an overload chain in
     * HaseWorker itself.
     */
    template<typename T_WorkerPolicy, typename T_WorkItem>
    struct HaseWorkItemDispatch;

    /**
     * @brief Uniform worker facade used by the policy-independent simulation loop.
     *
     * @tparam T_WorkerPolicy Stateful ownership and communication policy.
     */
    template<typename T_WorkerPolicy>
    class HaseWorker
    {
    public:
        /** @brief Construct a worker by taking ownership of its policy state. */
        explicit HaseWorker(T_WorkerPolicy policy) : m_policy(std::move(policy))
        {
        }

        /** @brief Return this worker's zero-based index in its worker group. */
        [[nodiscard]] unsigned workerIndex() const
        {
            return HaseWorkerDispatch<T_WorkerPolicy>::workerIndex(m_policy);
        }

        /** @brief Return the number of workers participating in the group. */
        [[nodiscard]] unsigned workerCount() const
        {
            return HaseWorkerDispatch<T_WorkerPolicy>::workerCount(m_policy);
        }

        /** @brief Return whether this worker owns root-only host responsibilities. */
        [[nodiscard]] bool isRoot() const
        {
            return HaseWorkerDispatch<T_WorkerPolicy>::isRoot(m_policy);
        }

        /** @brief Execute one complete work item through its dispatch specialization. */
        template<typename T_WorkItem>
        auto operator()(T_WorkItem&& workItem)
        {
            using T_Item = std::remove_cvref_t<T_WorkItem>;
            return HaseWorkItemDispatch<T_WorkerPolicy, T_Item>::run(m_policy, std::forward<T_WorkItem>(workItem));
        }

        /** @brief Scatter root-owned worker inputs according to the policy. */
        template<typename T_Value>
        auto scatter(T_Value&& value)
        {
            return HaseWorkerDispatch<T_WorkerPolicy>::scatter(m_policy, std::forward<T_Value>(value));
        }

        /** @brief Gather indexed results without replacing batch identity by worker identity. */
        template<typename T_Value>
        auto gather(T_Value&& value)
        {
            return HaseWorkerDispatch<T_WorkerPolicy>::gather(m_policy, std::forward<T_Value>(value));
        }

        /** @brief Reduce a genuinely additive or ordered diagnostic quantity. */
        template<typename T_Value, typename T_Reduction>
        auto reduce(T_Value&& value, T_Reduction reduction)
        {
            return HaseWorkerDispatch<T_WorkerPolicy>::reduce(
                m_policy,
                std::forward<T_Value>(value),
                std::move(reduction));
        }

        /** @brief Access policy state for worker-group lifecycle management. */
        [[nodiscard]] T_WorkerPolicy& policy()
        {
            return m_policy;
        }

    private:
        T_WorkerPolicy m_policy;
    };

    template<typename T_WorkerPolicy>
    HaseWorker(T_WorkerPolicy) -> HaseWorker<T_WorkerPolicy>;
} // namespace hase::core
