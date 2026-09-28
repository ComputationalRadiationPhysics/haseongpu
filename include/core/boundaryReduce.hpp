/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <tuple>
#include <utility>

namespace hase::core::detail
{
    // Workaround class to avoid a GCC 12 constraint cache error in the main translation unit.
    template<typename T_Queue, typename T_Output, typename T_Input>
    struct ReduceWrapper
    {
        static void run(
            T_Queue const& queue,
            alpaka::GetValueType_t<T_Output> neutral,
            T_Output& output,
            T_Input const& input)
        {
            alpaka::onHost::reduce(queue, neutral, output, std::plus{}, input);
        }
    };

    template<typename T_Spec, typename T_Value>
    struct EagerReduceWrapper
    {
        using Selector = ALPAKA_TYPEOF(alpaka::onHost::makeDeviceSelector(std::declval<T_Spec>()));
        using Device = ALPAKA_TYPEOF(std::declval<Selector&>().makeDevice(0u));
        using Queue = ALPAKA_TYPEOF(std::declval<Device&>().makeQueue(alpaka::queueKind::nonBlocking));
        using Buffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<T_Value>(std::declval<Device&>(), std::size_t{1u}));
        using View = ALPAKA_TYPEOF(std::declval<Buffer&>().getView().getSubView(alpaka::Vec{std::size_t{1u}}));
        static constexpr auto fn = &ReduceWrapper<Queue, Buffer, View>::run;
    };

    // Instantiate before the nested simulation dispatch to avoid GCC 12's constraint-cache crash.
    inline constexpr auto eagerReduceWrappers = []<typename... T_Specs>(std::tuple<T_Specs...>)
    {
        return std::tuple{
            EagerReduceWrapper<T_Specs, double>::fn...,
            EagerReduceWrapper<T_Specs, std::uint32_t>::fn...};
    }(alpaka::onHost::enabledDeviceSpecs);

    template<typename T_Queue, typename T_Output, typename T_Input>
    void reduce(T_Queue const& queue, alpaka::GetValueType_t<T_Output> neutral, T_Output& output, T_Input const& input)
    {
        using Function = decltype(&ReduceWrapper<T_Queue, T_Output, T_Input>::run);
        std::get<Function>(eagerReduceWrappers)(queue, neutral, output, input);
    }
} // namespace hase::core::detail
