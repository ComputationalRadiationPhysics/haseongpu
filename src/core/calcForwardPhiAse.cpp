/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <core/calcForwardPhiAse.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace hase::core
{
    ForwardPhiAseRawResult makeForwardRawResult(unsigned const volumeCount)
    {
        return ForwardPhiAseRawResult{
            std::vector<double>(volumeCount, 0.0),
            std::vector<double>(volumeCount, 0.0),
            std::vector<unsigned>(volumeCount, 0u),
            std::vector<unsigned>(volumeCount, 0u),
            0u,
            SrmStatus::DISABLED,
            0u,
            0.0,
            0u,
            0u};
    }

    double calcForwardBetaVolumeTotal(HostMesh const& hostMesh)
    {
        double total = 0.0;
        unsigned const count = std::min(
            static_cast<unsigned>(hostMesh.betaVolume.size()),
            static_cast<unsigned>(hostMesh.cellVolumes.size()));
        for(unsigned volume = 0u; volume < count; ++volume)
        {
            total += hostMesh.betaVolume.at(volume) * static_cast<double>(hostMesh.cellVolumes.at(volume));
        }
        return total;
    }

    void mergeForwardRawResult(ForwardPhiAseRawResult& target, ForwardPhiAseRawResult const& source)
    {
        if(target.scoreSum.empty())
        {
            target = source;
            return;
        }

        target.rayCount += source.rayCount;
        if(srmStatusPriority(source.srmStatus) > srmStatusPriority(target.srmStatus))
            target.srmStatus = source.srmStatus;
        target.srmPasses = std::max(target.srmPasses, source.srmPasses);
        target.srmRemainingFraction = std::max(target.srmRemainingFraction, source.srmRemainingFraction);
        target.srmMaxIterations = std::max(target.srmMaxIterations, source.srmMaxIterations);
        target.srmDivergenceStreak = std::max(target.srmDivergenceStreak, source.srmDivergenceStreak);
        for(unsigned volume = 0u; volume < target.scoreSum.size(); ++volume)
        {
            target.scoreSum.at(volume) += source.scoreSum.at(volume);
            target.scoreSquareSum.at(volume) += source.scoreSquareSum.at(volume);
            target.totalRays.at(volume) += source.totalRays.at(volume);
            target.droppedRays.at(volume) += source.droppedRays.at(volume);
        }
    }

    double calcForwardStandardError(
        double const scoreSum,
        double const scoreSquareSum,
        unsigned const rayCount,
        double const normalizationVolume,
        double const volumeSize)
    {
        if(rayCount < 2u || volumeSize <= 0.0 || normalizationVolume < 0.0 || !std::isfinite(scoreSum)
           || !std::isfinite(scoreSquareSum))
        {
            return std::numeric_limits<double>::max();
        }

        double const n = static_cast<double>(rayCount);
        double const varianceOfMean = (scoreSquareSum - scoreSum * scoreSum / n) / (n * (n - 1.0));
        double const volumeScale = normalizationVolume / volumeSize;
        return std::sqrt(std::max(0.0, varianceOfMean)) * volumeScale;
    }

    void finalizeForwardPhiAse(
        HostMesh const& hostMesh,
        ForwardPhiAseRawResult const& rawResult,
        Result& result)
    {
        unsigned const volumeCount = static_cast<unsigned>(rawResult.scoreSum.size());
        double const betaVolumeTotal = calcForwardBetaVolumeTotal(hostMesh);

        result = Result(
            std::vector(volumeCount, 0.0f),
            std::vector(volumeCount, 0.0),
            rawResult.totalRays,
            std::vector(volumeCount, 0.0),
            rawResult.droppedRays,
            rawResult.srmStatus,
            rawResult.srmPasses,
            rawResult.srmRemainingFraction,
            rawResult.srmMaxIterations,
            rawResult.srmDivergenceStreak);
        for(unsigned volume = 0u; volume < volumeCount; ++volume)
        {
            double const volumeSize = hostMesh.cellVolumes.at(volume);
            double const scoreSum = rawResult.scoreSum.at(volume);
            if(volumeSize > 0.0 && rawResult.rayCount > 0u)
            {
                double const estimate
                    = scoreSum * betaVolumeTotal / (static_cast<double>(rawResult.rayCount) * volumeSize);
                result.phiAse.at(volume) = static_cast<float>(estimate);
                result.mse.at(volume) = result.droppedRays[volume] == 0u ? calcForwardStandardError(
                                                                               scoreSum,
                                                                               rawResult.scoreSquareSum.at(volume),
                                                                               rawResult.rayCount,
                                                                               betaVolumeTotal,
                                                                               volumeSize)
                                                                         : std::numeric_limits<double>::max();
            }
            else
            {
                result.phiAse.at(volume) = 0.0f;
                result.mse.at(volume) = std::numeric_limits<double>::max();
            }
        }
    }
} // namespace hase::core
