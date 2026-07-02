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
    template<typename T>
    concept ComposeDerivativeBufferHandle = requires(T buffers) {
        buffers.betaCells;
        buffers.pumpedBeta;
        buffers.phiAse;
        buffers.activeMask;
        buffers.dndtPump;
        buffers.dndtAse;
        buffers.derivative;
    };

    struct ComposeDerivative
    {
        double sigmaAbsorption = 0.0;
        double sigmaEmission = 0.0;
        double tau = 1.0;
        double pumpDuration = 1.0;
        bool pumpEnabled = true;

        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            auto betaCells,
            auto pumpedBeta,
            auto phiAse,
            auto activeMask,
            auto dndtPump,
            auto dndtAse,
            auto derivative) const
        {
            for(auto [sample] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfSamples}))
            {
                unsigned const point = sample % mesh.numberOfPoints;
                double const beta = betaCells[sample];
                double const pumpTerm = pumpEnabled ? (pumpedBeta[sample] - beta) / pumpDuration : 0.0;
                double const gainPerDensity = beta * (sigmaEmission + sigmaAbsorption) - sigmaAbsorption;
                double const aseTerm
                    = activeMask[point] != 0u ? gainPerDensity * static_cast<double>(phiAse[sample]) : 0.0;

                dndtPump[sample] = pumpTerm;
                dndtAse[sample] = aseTerm;
                derivative[sample] = pumpTerm - aseTerm - beta / tau;
            }
        }
    };

    void enqueueComposeDerivative(
        auto& devBundle,
        hase::concepts::Queue auto const& queue,
        auto const& mesh,
        double sigmaAbsorption,
        double sigmaEmission,
        double tau,
        double pumpDuration,
        bool pumpEnabled,
        ComposeDerivativeBufferHandle auto& buffers)
    {
        auto frameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfSamples});
        queue.enqueue(
            frameSpec,
            alpaka::KernelBundle{
                ComposeDerivative{sigmaAbsorption, sigmaEmission, tau, pumpDuration, pumpEnabled},
                mesh,
                buffers.betaCells,
                buffers.pumpedBeta,
                buffers.phiAse,
                buffers.activeMask,
                buffers.dndtPump,
                buffers.dndtAse,
                buffers.derivative});
    }

} // namespace hase::kernels
