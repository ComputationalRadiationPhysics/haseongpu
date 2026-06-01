/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once
#include <alpaka/alpaka.hpp>

#include <concepts/concepts.hpp>

#include <vector>

namespace hase::alpakaUtils
{

    struct ToDevice
    {
        template<typename T>
        struct Op
        {
            auto operator()(concepts::Queue auto const& queue, T const& inputView)
            {
                auto deviceBuffer = alpaka::onHost::allocLike(queue.getDevice(), inputView);
                alpaka::onHost::memcpy(queue, deviceBuffer, inputView);
                alpaka::onHost::wait(queue);
                return deviceBuffer;
            }
        };
    };

    // A fully variadic partial specialization such as alpaka::View<TArgs...>
    // triggers an nvcc segmentation fault during compilation (cuda 12.9).
    // Keeping at least one explicit template parameter, e.g. alpaka::View<A, TArgs...>,
    // avoids the crash.
    template<typename A, typename... TArgs>
    struct ToDevice::Op<alpaka::View<A, TArgs...>>
    {
        auto operator()(concepts::Queue auto const& queue, alpaka::View<A, TArgs...> const& inputView)
        {
            auto deviceBuffer = alpaka::onHost::allocLike(queue.getDevice(), inputView);
            alpaka::onHost::memcpy(queue, deviceBuffer, inputView);
            alpaka::onHost::wait(queue);
            return deviceBuffer;
        }
    };

    template<typename A>
    struct ToDevice::Op<std::vector<A>>
    {
        auto operator()(concepts::Queue auto const& queue, std::vector<A> const& inputBuffer)
        {
            auto deviceBuffer = alpaka::onHost::allocLike(queue.getDevice(), inputBuffer);
            alpaka::onHost::memcpy(queue, deviceBuffer, inputBuffer);
            alpaka::onHost::wait(queue);
            return deviceBuffer;
        }
    };

    auto toDevice(concepts::Queue auto const& queue, auto const& inputBuffer)
    {
        return hase::alpakaUtils::ToDevice::Op<ALPAKA_TYPEOF(inputBuffer)>{}(queue, inputBuffer);
    }
} // namespace hase::alpakaUtils
