/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/utils.hpp>
#include <concepts/concepts.hpp>
#include <core/mesh.hpp>

#include <cmath>

namespace hase::kernels
{
    struct BuildActivePointMask
    {
        ALPAKA_FN_ACC void operator()(auto const& acc, hase::core::DeviceMeshView const mesh, auto activeMask) const
        {
            for(auto [point] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfPoints}))
            {
                activeMask[point] = 1u;
            }
        }
    };

    void enqueueBuildActivePointMask(
        auto& devBundle,
        hase::concepts::Queue auto const& queue,
        auto const& mesh,
        auto& activeMask)
    {
        alpaka::onHost::fill(queue, activeMask, 0u);
        auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfPoints});
        queue.enqueue(frameSpec, alpaka::KernelBundle{BuildActivePointMask{}, mesh, activeMask});
    }

} // namespace hase::kernels
