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
                unsigned active = 0u;
                for(unsigned triangle = 0u; triangle < mesh.numberOfTriangles && active == 0u; ++triangle)
                {
                    if(mesh.claddingCellTypes[triangle] == mesh.claddingNumber)
                    {
                        continue;
                    }

                    for(unsigned vertex = 0u; vertex < 3u; ++vertex)
                    {
                        if(mesh.trianglePointIndices[triangle + vertex * mesh.numberOfTriangles] == point)
                        {
                            active = 1u;
                        }
                    }
                }
                activeMask[point] = active;
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
            alpaka::Vec{mesh.numberOfTriangles});
        queue.enqueue(frameSpec, alpaka::KernelBundle{BuildActivePointMask{}, mesh, activeMask});
    }

} // namespace hase::kernels
