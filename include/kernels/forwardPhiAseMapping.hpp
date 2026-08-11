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
#include <kernels/forward/accumulation.hpp>

#include <limits>

namespace hase::kernels
{
    struct BuildBetaVolumeWeights
    {
        ALPAKA_FN_ACC void operator()(
            auto const& acc,
            core::DeviceMeshView const mesh,
            auto betaVolume,
            auto betaVolumeWeights) const
        {
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                betaVolumeWeights[cell] = betaVolume[cell] * static_cast<double>(mesh.cellVolumes[cell]);
            }
        }
    };

    struct CaptureBetaVolumeTotal
    {
        ALPAKA_FN_ACC void operator()(
            auto const&,
            unsigned const numberOfCells,
            auto betaVolumePrefix,
            auto betaVolumeTotal) const
        {
            betaVolumeTotal[0u] = numberOfCells == 0u ? 0.0 : betaVolumePrefix[numberOfCells - 1u];
        }
    };

    struct FinalizeForwardVolumePhiAse
    {
        unsigned rayCount;
        unsigned batchCount;
        double betaVolumeTotal;
        double fluorescenceRate;
        double sigmaA;
        double sigmaE;

        template<
            typename T_Acc,
            typename T_VertexBatchScoreSum,
            typename T_RseBatchRayCounts,
            typename T_LumpedMaterialVertexVolume,
            typename T_DroppedRays,
            typename T_VolumePhiAse,
            typename T_StandardError,
            typename T_RelativeStandardError,
            typename T_VolumeDndtAse>
        ALPAKA_FN_ACC void operator()(
            T_Acc const& acc,
            core::DeviceMeshView const mesh,
            T_VertexBatchScoreSum vertexBatchScoreSum,
            T_RseBatchRayCounts rseBatchRayCounts,
            T_LumpedMaterialVertexVolume lumpedMaterialVertexVolume,
            T_DroppedRays droppedRays,
            T_VolumePhiAse volumePhiAse,
            T_StandardError standardError,
            T_RelativeStandardError relativeStandardError,
            T_VolumeDndtAse volumeDndtAse) const
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
                    unsigned const materialVertexOffset
                        = mesh.getCellType(cell) == mesh.claddingNumber ? mesh.numberOfMeshPoints : 0u;
                    double scoreSum = 0.0;
                    double batchMeanSum = 0.0;
                    double batchMeanSquareSum = 0.0;
                    unsigned activeBatchCount = 0u;
                    for(unsigned batch = 0u; batch < batchCount; ++batch)
                    {
                        double batchScoreDensity = 0.0;
                        for(unsigned localVertex = 0u; localVertex < mesh.numberOfCellVertices; ++localVertex)
                        {
                            unsigned const materialVertex
                                = materialVertexOffset
                                  + mesh.cellPointIndices[cell * mesh.numberOfCellVertices + localVertex];
                            double const vertexVolume = lumpedMaterialVertexVolume[materialVertex];
                            unsigned const vertex = batch * (2u * mesh.numberOfMeshPoints) + materialVertex;
                            batchScoreDensity += vertexVolume > 0.0 ? vertexBatchScoreSum[vertex] / vertexVolume : 0.0;
                        }
                        batchScoreDensity /= static_cast<double>(mesh.numberOfCellVertices);
                        double const batchScore = batchScoreDensity * volume;
                        scoreSum += batchScore;
                        if(rseBatchRayCounts[batch] == 0u)
                            continue;
                        double const batchMean = batchScore / static_cast<double>(rseBatchRayCounts[batch]);
                        batchMeanSum += batchMean;
                        batchMeanSquareSum += batchMean * batchMean;
                        ++activeBatchCount;
                    }
                    estimate = scoreSum * betaVolumeTotal / (static_cast<double>(rayCount) * volume);
                    if(droppedRays[cell] == 0u && activeBatchCount >= 2u)
                    {
                        double const count = static_cast<double>(activeBatchCount);
                        double const batchMean = batchMeanSum / count;
                        if(batchMean == 0.0)
                        {
                            relativeError = std::numeric_limits<double>::quiet_NaN();
                            absoluteError = 0.0;
                        }
                        else
                        {
                            double const sampleVariance = alpaka::math::max(
                                0.0,
                                (batchMeanSquareSum - batchMeanSum * batchMeanSum / count) / (count - 1.0));
                            relativeError = alpaka::math::sqrt(sampleVariance / count) / alpaka::math::abs(batchMean);
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

    template<
        typename T_DevBundle,
        typename T_Queue,
        typename T_VertexBatchScoreSum,
        typename T_RseBatchRayCounts,
        typename T_LumpedMaterialVertexVolume,
        typename T_DroppedRays,
        typename T_VolumePhiAse,
        typename T_StandardError,
        typename T_RelativeStandardError,
        typename T_VolumeDndtAse>
    void enqueueFinalizeForwardCellPhiAse(
        T_DevBundle& devBundle,
        T_Queue const& queue,
        core::DeviceMeshView const mesh,
        T_VertexBatchScoreSum const& vertexBatchScoreSum,
        T_RseBatchRayCounts const& rseBatchRayCounts,
        T_LumpedMaterialVertexVolume const& lumpedMaterialVertexVolume,
        T_DroppedRays const& droppedRays,
        T_VolumePhiAse& volumePhiAse,
        T_StandardError& standardError,
        T_RelativeStandardError& relativeStandardError,
        T_VolumeDndtAse& volumeDndtAse,
        unsigned rayCount,
        unsigned batchCount,
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
                FinalizeForwardVolumePhiAse{rayCount, batchCount, betaVolumeTotal, fluorescenceRate, sigmaA, sigmaE},
                mesh,
                vertexBatchScoreSum,
                rseBatchRayCounts,
                lumpedMaterialVertexVolume,
                droppedRays,
                volumePhiAse,
                standardError,
                relativeStandardError,
                volumeDndtAse});
    }

} // namespace hase::kernels
