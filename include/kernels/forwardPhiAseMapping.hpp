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
        alpaka::Vec<unsigned, hase::kernels::forward::forwardRseBatchCount> rseBatchRayCounts;
        double betaVolumeTotal;
        double fluorescenceRate;
        double sigmaA;
        double sigmaE;

        template<
            typename T_Acc,
            typename T_VertexBatchScoreSum,
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
                    for(unsigned batch = 0u; batch < hase::kernels::forward::forwardRseBatchCount; ++batch)
                    {
                        double batchScoreDensity = 0.0;
                        for(unsigned localVertex = 0u; localVertex < mesh.numberOfCellVertices; ++localVertex)
                        {
                            unsigned const materialVertex
                                = materialVertexOffset
                                  + mesh.cellPointIndices[cell * mesh.numberOfCellVertices + localVertex];
                            double const vertexVolume = lumpedMaterialVertexVolume[materialVertex];
                            unsigned const vertex
                                = batch * (2u * mesh.numberOfMeshPoints) + materialVertex;
                            batchScoreDensity
                                += vertexVolume > 0.0 ? vertexBatchScoreSum[vertex] / vertexVolume : 0.0;
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

    template<
        typename T_DevBundle,
        typename T_Queue,
        typename T_VertexBatchScoreSum,
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
        T_LumpedMaterialVertexVolume const& lumpedMaterialVertexVolume,
        T_DroppedRays const& droppedRays,
        T_VolumePhiAse& volumePhiAse,
        T_StandardError& standardError,
        T_RelativeStandardError& relativeStandardError,
        T_VolumeDndtAse& volumeDndtAse,
        unsigned rayCount,
        alpaka::Vec<unsigned, hase::kernels::forward::forwardRseBatchCount> rseBatchRayCounts,
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
                FinalizeForwardVolumePhiAse{
                    rayCount,
                    rseBatchRayCounts,
                    betaVolumeTotal,
                    fluorescenceRate,
                    sigmaA,
                    sigmaE},
                mesh,
                vertexBatchScoreSum,
                lumpedMaterialVertexVolume,
                droppedRays,
                volumePhiAse,
                standardError,
                relativeStandardError,
                volumeDndtAse});
    }

} // namespace hase::kernels
