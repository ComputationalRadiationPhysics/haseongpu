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
    struct MapPointBetaToPrismBeta
    {
        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            auto betaCells,
            auto betaVolume) const
        {
            for(auto [prism] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfPrisms}))
            {
                if(!mesh.samplePointsAreMeshPoints)
                {
                    betaVolume[prism] = betaCells[prism];
                    continue;
                }

                double sum = 0.0;
                for(unsigned vertex = 0u; vertex < mesh.numberOfCellVertices; ++vertex)
                {
                    unsigned const point = mesh.cellPointIndices[prism * mesh.numberOfCellVertices + vertex];
                    sum += betaCells[point];
                }
                betaVolume[prism] = sum / static_cast<double>(mesh.numberOfCellVertices);
            }
        }
    };

    void enqueueMapPointBetaToPrismBeta(
        auto& devBundle,
        hase::concepts::Queue auto const& queue,
        auto const& mesh,
        auto& betaCells,
        auto& betaVolume)
    {
        auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfPrisms});
        queue.enqueue(frameSpec, alpaka::KernelBundle{MapPointBetaToPrismBeta{}, mesh, betaCells, betaVolume});
    }

} // namespace hase::kernels
