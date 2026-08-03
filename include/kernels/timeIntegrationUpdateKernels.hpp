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
    struct AddScaled
    {
        double scale = 1.0;

        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            auto base,
            auto slope,
            auto out) const
        {
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                out[cell] = base[cell] + scale * slope[cell];
            }
        }
    };

    struct CombineHeun
    {
        double timeStep = 0.0;

        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            auto base,
            auto first,
            auto second,
            auto out) const
        {
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                out[cell] = base[cell] + 0.5 * timeStep * (first[cell] + second[cell]);
            }
        }
    };

    struct CombineRungeKutta4
    {
        double timeStep = 0.0;

        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            auto base,
            auto k1,
            auto k2,
            auto k3,
            auto k4,
            auto out) const
        {
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                out[cell] = base[cell] + (timeStep / 6.0) * (k1[cell] + 2.0 * k2[cell] + 2.0 * k3[cell] + k4[cell]);
            }
        }
    };

    struct ExponentialEulerUpdate
    {
        double timeStep = 0.0;
        double tau = 1.0;

        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            auto betaVolume,
            auto dndtPump,
            auto dndtAse,
            auto out) const
        {
            double const decay = alpaka::math::exp(-timeStep / tau);
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                double const source = dndtPump[cell] - dndtAse[cell];
                out[cell] = tau * source * (1.0 - decay) + betaVolume[cell] * decay;
            }
        }
    };

    struct ClipBeta
    {
        ALPAKA_FN_ACC void operator()(auto const& acc, hase::core::DeviceMeshView const mesh, auto betaVolume) const
        {
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                betaVolume[cell] = alpaka::math::min(1.0, alpaka::math::max(0.0, betaVolume[cell]));
            }
        }
    };

} // namespace hase::kernels
