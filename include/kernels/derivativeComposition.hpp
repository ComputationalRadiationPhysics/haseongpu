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
#include <concepts/concepts.hpp>
#include <core/mesh.hpp>

namespace hase::kernels
{
    template<typename T>
    concept ComposeDerivativeBufferHandle = requires(T buffers) {
        buffers.betaVolume;
        buffers.phiAse;
        buffers.dndtPump;
        buffers.dndtAse;
        buffers.derivative;
    };

    struct ComposeDerivative
    {
        double sigmaAbsorption = 0.0;
        double sigmaEmission = 0.0;
        double tau = 1.0;
        bool pumpEnabled = true;

        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            core::DeviceMeshView const mesh,
            auto betaVolume,
            auto phiAse,
            auto dndtPump,
            auto dndtAse,
            auto derivative) const
        {
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                if(mesh.getCellType(cell) == mesh.claddingNumber)
                {
                    dndtPump[cell] = 0.0;
                    dndtAse[cell] = 0.0;
                    derivative[cell] = 0.0;
                    continue;
                }

                double const pumpTerm = pumpEnabled ? dndtPump[cell] : 0.0;
                double const gainPerDensity = betaVolume[cell] * (sigmaEmission + sigmaAbsorption) - sigmaAbsorption;
                double const aseTerm = gainPerDensity * static_cast<double>(phiAse[cell]);
                if(!pumpEnabled)
                    dndtPump[cell] = 0.0;
                dndtAse[cell] = aseTerm;
                derivative[cell] = pumpTerm - aseTerm - betaVolume[cell] / tau;
            }
        }
    };

    void enqueueComposeDerivative(
        auto& devBundle,
        concepts::Queue auto const& queue,
        auto const& mesh,
        double sigmaAbsorption,
        double sigmaEmission,
        double tau,
        bool pumpEnabled,
        ComposeDerivativeBufferHandle auto& buffers)
    {
        auto frameSpec = alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfCells});
        queue.enqueue(
            frameSpec,
            alpaka::KernelBundle{
                ComposeDerivative{sigmaAbsorption, sigmaEmission, tau, pumpEnabled},
                mesh,
                buffers.betaVolume,
                buffers.phiAse,
                buffers.dndtPump,
                buffers.dndtAse,
                buffers.derivative});
    }

} // namespace hase::kernels
