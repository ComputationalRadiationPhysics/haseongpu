/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <alpaka/math.hpp>

#include <core/calcForwardPhiAse.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace hase::core
{
    ForwardPhiAseRawResult makeForwardRawResult(unsigned const volumeCount, unsigned const vertexCount)
    {
        return ForwardPhiAseRawResult{
            std::vector<double>(hase::kernels::forward::forwardRseBatchCount * 2u * vertexCount, 0.0),
            {},
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
        if(target.vertexBatchScoreSum.empty())
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
        if(target.vertexBatchScoreSum.size() != source.vertexBatchScoreSum.size())
            throw std::runtime_error("cannot merge forward ASE results with different vertex counts");
        for(unsigned vertex = 0u; vertex < target.vertexBatchScoreSum.size(); ++vertex)
            target.vertexBatchScoreSum.at(vertex) += source.vertexBatchScoreSum.at(vertex);
        for(unsigned batch = 0u; batch < hase::kernels::forward::forwardRseBatchCount; ++batch)
            target.rseBatchRayCounts.at(batch) += source.rseBatchRayCounts.at(batch);
        for(unsigned volume = 0u; volume < target.totalRays.size(); ++volume)
        {
            target.totalRays.at(volume) += source.totalRays.at(volume);
            target.droppedRays.at(volume) += source.droppedRays.at(volume);
        }
    }

    double calcForwardRelativeStandardError(
        double const scoreSum,
        double const scoreSquareSum,
        unsigned const rayCount)
    {
        if(rayCount < 2u || !alpaka::math::isfinite(scoreSum) || !alpaka::math::isfinite(scoreSquareSum))
        {
            return std::numeric_limits<double>::max();
        }
        if(scoreSum == 0.0)
        {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double const n = static_cast<double>(rayCount);
        double const relativeVariance = (n * scoreSquareSum / (scoreSum * scoreSum) - 1.0) / n;
        return std::sqrt(std::max(0.0, relativeVariance));
    }

    double calcForwardStandardError(
        double const scoreSum,
        double const scoreSquareSum,
        unsigned const rayCount,
        double const normalizationVolume,
        double const volumeSize)
    {
        if(rayCount < 2u || volumeSize <= 0.0 || normalizationVolume < 0.0 || !alpaka::math::isfinite(scoreSum)
           || !alpaka::math::isfinite(scoreSquareSum))
        {
            return std::numeric_limits<double>::max();
        }

        double const relativeStandardError = calcForwardRelativeStandardError(scoreSum, scoreSquareSum, rayCount);
        if(alpaka::math::isnan(relativeStandardError))
        {
            return 0.0;
        }
        if(!alpaka::math::isfinite(relativeStandardError))
        {
            return std::numeric_limits<double>::max();
        }

        double const volumeScale = normalizationVolume / volumeSize;
        double const estimate = scoreSum * volumeScale / rayCount;
        return relativeStandardError * std::abs(estimate);
    }

    void finalizeForwardPhiAse(HostMesh const& hostMesh, ForwardPhiAseRawResult const& rawResult, Result& result)
    {
        finalizeForwardPhiAse(hostMesh, rawResult, calcForwardBetaVolumeTotal(hostMesh), result);
    }

    void finalizeForwardPhiAse(
        HostMesh const& hostMesh,
        ForwardPhiAseRawResult const& rawResult,
        double const betaVolumeTotal,
        Result& result)
    {
        unsigned const volumeCount = hostMesh.numberOfCells;
        unsigned const materialVertexCount = 2u * hostMesh.numberOfMeshPoints;
        if(rawResult.vertexBatchScoreSum.size() != hase::kernels::forward::forwardRseBatchCount * materialVertexCount)
            throw std::runtime_error("forward ASE vertex score count does not match the mesh");
        std::array<std::vector<double>, hase::kernels::forward::forwardRseBatchCount> cellBatchScoreDensity;
        for(unsigned batch = 0u; batch < hase::kernels::forward::forwardRseBatchCount; ++batch)
        {
            auto const begin = rawResult.vertexBatchScoreSum.cbegin() + batch * materialVertexCount;
            std::vector<double> const vertexBatch(begin, begin + materialVertexCount);
            cellBatchScoreDensity.at(batch)
                = hase::kernels::accumulateMaterialVertexIntegralsToCellDensities(hostMesh, vertexBatch);
        }
        result = Result(
            std::vector(volumeCount, 0.0f),
            std::vector(volumeCount, 0.0),
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
            if(volumeSize > 0.0 && rawResult.rayCount > 0u)
            {
                double scoreSum = 0.0;
                double batchMeanSum = 0.0;
                double batchMeanSquareSum = 0.0;
                unsigned activeBatchCount = 0u;
                for(unsigned batch = 0u; batch < hase::kernels::forward::forwardRseBatchCount; ++batch)
                {
                    unsigned const batchRayCount = rawResult.rseBatchRayCounts.at(batch);
                    double const batchScore = cellBatchScoreDensity.at(batch).at(volume) * volumeSize;
                    scoreSum += batchScore;
                    if(batchRayCount == 0u)
                        continue;
                    double const batchMean = batchScore / static_cast<double>(batchRayCount);
                    batchMeanSum += batchMean;
                    batchMeanSquareSum += batchMean * batchMean;
                    ++activeBatchCount;
                }
                double const estimate
                    = scoreSum * betaVolumeTotal / (static_cast<double>(rawResult.rayCount) * volumeSize);
                result.phiAse.at(volume) = static_cast<float>(estimate);
                double relativeError = std::numeric_limits<double>::max();
                if(result.droppedRays[volume] == 0u && activeBatchCount >= 2u)
                {
                    double const count = static_cast<double>(activeBatchCount);
                    double const batchMean = batchMeanSum / count;
                    if(batchMean == 0.0)
                        relativeError = std::numeric_limits<double>::quiet_NaN();
                    else
                    {
                        double const sampleVariance = std::max(
                            0.0,
                            (batchMeanSquareSum - batchMeanSum * batchMeanSum / count) / (count - 1.0));
                        relativeError = std::sqrt(sampleVariance / count) / std::abs(batchMean);
                    }
                }
                result.relativeStandardError.at(volume) = relativeError;
                result.standardError.at(volume)
                    = result.droppedRays[volume] != 0u || !std::isfinite(relativeError)
                          ? (std::isnan(relativeError) ? 0.0 : std::numeric_limits<double>::max())
                          : relativeError * std::abs(estimate);
            }
            else
            {
                result.phiAse.at(volume) = 0.0f;
                result.standardError.at(volume) = std::numeric_limits<double>::max();
                result.relativeStandardError.at(volume) = std::numeric_limits<double>::max();
            }
        }
    }
} // namespace hase::core
