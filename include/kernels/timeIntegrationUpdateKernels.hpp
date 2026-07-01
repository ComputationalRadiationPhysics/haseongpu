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

        ALPAKA_FN_ACC void operator()(auto const& acc, hase::core::DeviceMeshView const mesh, auto base, auto slope, auto out)
            const
        {
            for(auto [sample] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfSamples}))
            {
                out[sample] = base[sample] + scale * slope[sample];
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
            for(auto [sample] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfSamples}))
            {
                out[sample] = base[sample] + 0.5 * timeStep * (first[sample] + second[sample]);
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
            for(auto [sample] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfSamples}))
            {
                out[sample] = base[sample] + (timeStep / 6.0) * (k1[sample] + 2.0 * k2[sample] + 2.0 * k3[sample] + k4[sample]);
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
            auto betaCells,
            auto dndtPump,
            auto dndtAse,
            auto out) const
        {
            double const decay = alpaka::math::exp(-timeStep / tau);
            for(auto [sample] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfSamples}))
            {
                double const source = dndtPump[sample] - dndtAse[sample];
                out[sample] = tau * source * (1.0 - decay) + betaCells[sample] * decay;
            }
        }
    };

    struct ClipBeta
    {
        ALPAKA_FN_ACC void operator()(auto const& acc, hase::core::DeviceMeshView const mesh, auto betaCells) const
        {
            for(auto [sample] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfSamples}))
            {
                betaCells[sample] = alpaka::math::min(1.0, alpaka::math::max(0.0, betaCells[sample]));
            }
        }
    };

} // namespace hase::kernels
