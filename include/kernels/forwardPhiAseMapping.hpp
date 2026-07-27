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
#include <core/mesh.hpp>

#include <limits>

namespace hase::kernels
{
    struct BuildBetaVolumePrefix
    {
        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            core::DeviceMeshView const mesh,
            auto betaVolume,
            auto betaVolumePrefix,
            auto betaVolumeTotal) const
        {
            for(auto [worker] :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{1u}))
            {
                double total = 0.0;
                for(unsigned cell = 0u; cell < mesh.numberOfCells; ++cell)
                {
                    total += betaVolume[cell] * static_cast<double>(mesh.cellVolumes[cell]);
                    betaVolumePrefix[cell] = total;
                }
                betaVolumeTotal[worker] = total;
            }
        }
    };

    struct FinalizeForwardVolumePhiAse
    {
        unsigned rayCount;
        double betaVolumeTotal;
        double fluorescenceRate;
        double sigmaA;
        double sigmaE;

        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            core::DeviceMeshView const mesh,
            auto scoreSum,
            auto scoreSquareSum,
            auto droppedRays,
            auto volumePhiAse,
            auto standardError,
            auto relativeStandardError,
            auto volumeDndtAse) const
        {
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                double const volume = static_cast<double>(mesh.cellVolumes[cell]);
                double const maximum = std::numeric_limits<double>::max();
                double relativeError = maximum;
                double absoluteError = maximum;
                double estimate = 0.0;
                if(rayCount > 0u && volume > 0.0)
                {
                    estimate = scoreSum[cell] * betaVolumeTotal / (static_cast<double>(rayCount) * volume);
                    if(droppedRays[cell] == 0u && rayCount >= 2u && alpaka::math::isfinite(scoreSum[cell])
                       && alpaka::math::isfinite(scoreSquareSum[cell]))
                    {
                        if(scoreSum[cell] == 0.0)
                        {
                            relativeError = std::numeric_limits<double>::quiet_NaN();
                            absoluteError = 0.0;
                        }
                        else
                        {
                            double const n = static_cast<double>(rayCount);
                            double const relativeVariance
                                = (n * scoreSquareSum[cell] / (scoreSum[cell] * scoreSum[cell]) - 1.0) / n;
                            relativeError = alpaka::math::sqrt(alpaka::math::max(0.0, relativeVariance));
                            absoluteError = relativeError * alpaka::math::abs(estimate) * fluorescenceRate;
                        }
                    }
                }
                float const phiAse = static_cast<float>(estimate * fluorescenceRate);
                volumePhiAse[cell] = phiAse;
                standardError[cell] = absoluteError;
                relativeStandardError[cell] = relativeError;
                double const gainPerDensity = mesh.betaVolume[cell] * (sigmaE + sigmaA) - sigmaA;
                volumeDndtAse[cell] = gainPerDensity * static_cast<double>(phiAse);
            }
        }
    };

    void enqueueBuildBetaVolumePrefix(
        auto& devBundle,
        auto const& queue,
        core::DeviceMeshView const mesh,
        auto const& betaVolume,
        auto& betaVolumePrefix,
        auto& betaVolumeTotal)
    {
        auto frameSpec
            = hase::alpakaUtils::getFrameSpec<uint32_t>(devBundle.device, devBundle.executor, alpaka::Vec{1u});
        queue.enqueue(
            frameSpec,
            alpaka::KernelBundle{BuildBetaVolumePrefix{}, mesh, betaVolume, betaVolumePrefix, betaVolumeTotal});
    }

    void enqueueFinalizeForwardCellPhiAse(
        auto& devBundle,
        auto const& queue,
        core::DeviceMeshView const mesh,
        auto const& scoreSum,
        auto const& scoreSquareSum,
        auto const& droppedRays,
        auto& volumePhiAse,
        auto& standardError,
        auto& relativeStandardError,
        auto& volumeDndtAse,
        unsigned rayCount,
        double betaVolumeTotal,
        double fluorescenceRate,
        double sigmaA,
        double sigmaE)
    {
        auto cellFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfCells});
        queue.enqueue(
            cellFrameSpec,
            alpaka::KernelBundle{
                FinalizeForwardVolumePhiAse{rayCount, betaVolumeTotal, fluorescenceRate, sigmaA, sigmaE},
                mesh,
                scoreSum,
                scoreSquareSum,
                droppedRays,
                volumePhiAse,
                standardError,
                relativeStandardError,
                volumeDndtAse});
    }

} // namespace hase::kernels
