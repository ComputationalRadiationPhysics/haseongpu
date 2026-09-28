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
#include <data/TraceData.hpp>
#include <kernels/forward/accumulation.hpp>
#include <kernels/forward/rayPopulationStatistics.hpp>

#include <concepts>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace hase::kernels
{
    /** @brief Device operation building per-cell spontaneous-source weights. */
    struct BuildSourceStrengthWeights
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            data::TraceView const mesh,
            alpaka::concepts::IView<double> auto sourceStrengthWeights) const
        {
            for(auto [cell] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{mesh.numberOfCells}))
            {
                sourceStrengthWeights[cell] = mesh.isActive(cell)
                                                  ? mesh.getBetaVolume(cell) * mesh.getCellVolume(cell)
                                                        * mesh.activeIonDensity(cell) / mesh.fluorescenceLifetime(cell)
                                                  : 0.0;
            }
        }
    };

    /** @brief Transform domain-grouped global cell ids into spontaneous-source weights. */
    struct DomainSourceStrength
    {
        data::TraceView mesh;

        ALPAKA_FN_ACC auto operator()(alpaka::concepts::Simd auto const& globalCells) const
        {
            using T_Cells = std::remove_cvref_t<decltype(globalCells)>;
            static_assert(std::same_as<alpaka::GetValueType_t<T_Cells>, std::uint32_t>);
            return alpaka::Simd<double, T_Cells::width()>{
                [&](auto const lane)
                {
                    auto const cell = globalCells[decltype(lane)::value];
                    return mesh.isActive(cell) ? mesh.getBetaVolume(cell) * mesh.getCellVolume(cell)
                                                     * mesh.activeIonDensity(cell) / mesh.fluorescenceLifetime(cell)
                                               : 0.0;
                }};
        }
    };

    /** @brief Read each domain's total directly from its local CDF endpoint. */
    struct ComputeDomainSourceStrengthTotals
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            alpaka::concepts::IView<std::uint32_t> auto const offsets,
            alpaka::concepts::IView<double> auto const sourceStrengthPrefix,
            alpaka::concepts::IView<double> auto sourceStrengthTotals,
            std::uint32_t const domainCount) const
        {
            for(auto [domain] :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{domainCount}))
            {
                auto const begin = offsets[domain];
                auto const end = offsets[domain + 1u];
                sourceStrengthTotals[domain] = begin == end ? 0.0 : sourceStrengthPrefix[end - 1u];
            }
        }
    };

    /** @brief Device operation copying the final prefix sum into scalar storage. */
    struct CaptureSourceStrengthTotal
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const&,
            unsigned const numberOfCells,
            alpaka::concepts::IView<double> auto sourceStrengthPrefix,
            alpaka::concepts::IView<double> auto sourceStrengthTotal) const
        {
            sourceStrengthTotal[0u] = numberOfCells == 0u ? 0.0 : sourceStrengthPrefix[numberOfCells - 1u];
        }
    };

    /** @brief Normalize forward scores and derive cell PhiASE, RSE, and ASE rate. */
    struct FinalizeForwardVolumePhiAse
    {
        unsigned rayCount;
        unsigned numIndependentRayPopulations;
        double sourceStrengthTotal;

        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            data::TraceView const mesh,
            alpaka::concepts::IView<double const> auto vertexPopulationScoreSum,
            alpaka::concepts::IView<std::uint32_t> auto rayPopulationRayCounts,
            alpaka::concepts::IView<std::uint32_t const> auto droppedRays,
            alpaka::concepts::IView<float> auto volumePhiAse,
            alpaka::concepts::IView<double> auto standardError,
            alpaka::concepts::IView<double> auto relativeStandardError,
            alpaka::concepts::IView<double> auto volumeDndtAse) const
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
                    unsigned const materialVertexOffset = mesh.getMaterialId(cell) * mesh.numberOfMeshPoints;
                    forward::RayPopulationStatistics statistics;
                    for(unsigned batch = 0u; batch < numIndependentRayPopulations; ++batch)
                    {
                        double batchScoreDensity = 0.0;
                        for(unsigned localVertex = 0u; localVertex < mesh.numberOfCellVertices; ++localVertex)
                        {
                            unsigned const materialVertex
                                = materialVertexOffset
                                  + mesh.cellPointIndices[cell * mesh.numberOfCellVertices + localVertex];
                            double const vertexVolume = mesh.lumpedMaterialVertexVolumes[materialVertex];
                            std::size_t const vertex
                                = static_cast<std::size_t>(batch) * mesh.numberOfMaterials * mesh.numberOfMeshPoints
                                  + materialVertex;
                            batchScoreDensity
                                += vertexVolume > 0.0 ? vertexPopulationScoreSum[vertex] / vertexVolume : 0.0;
                        }
                        batchScoreDensity /= static_cast<double>(mesh.numberOfCellVertices);
                        statistics.add(batchScoreDensity, rayPopulationRayCounts[batch]);
                    }
                    auto const summary = statistics.finalize(sourceStrengthTotal, droppedRays[cell] != 0u);
                    estimate = summary.value;
                    relativeError = summary.relativeStandardError;
                    absoluteError = summary.standardError;
                }
                float const phiAse = static_cast<float>(estimate);
                volumePhiAse[cell] = phiAse;
                standardError[cell] = absoluteError;
                relativeStandardError[cell] = relativeError;
                unsigned const material = mesh.getMaterialId(cell);
                double const gainPerDensity
                    = mesh.isActive(cell)
                          ? mesh.betaVolume[cell]
                                    * (mesh.materialPeakEmission[material] + mesh.materialPeakAbsorption[material])
                                - mesh.materialPeakAbsorption[material]
                          : 0.0;
                volumeDndtAse[cell] = gainPerDensity * estimate;
            }
        }
    };

    /**
     * @brief Enqueue final normalization of accumulated forward-ray statistics.
     * @param devBundle Device and executor used to build the cell launch.
     * @param queue Queue receiving the finalization kernel.
     * @param mesh Device-resident trace view and lumped material-vertex volumes.
     * @param vertexPopulationScoreSum Raw score sums indexed by batch and material vertex.
     * @param rayPopulationRayCounts Number of histories contributing to each batch.
     * @param droppedRays Per-cell traversal-failure counts.
     * @param volumePhiAse Cell PhiASE output.
     * @param standardError Cell absolute standard-error output.
     * @param relativeStandardError Cell relative standard-error output.
     * @param volumeDndtAse Cell ASE population-rate output.
     * @param rayCount Total histories represented by all batches.
     * @param numIndependentRayPopulations Number of statistical batches.
     * @param sourceStrengthTotal Integral used to normalize the sampled source.
     */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Executor>
    void enqueueFinalizeForwardCellPhiAse(
        alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
        concepts::Queue auto const& queue,
        data::TraceView const mesh,
        alpaka::concepts::IBuffer<double> auto const& vertexPopulationScoreSum,
        alpaka::concepts::IView<std::uint32_t> auto const& rayPopulationRayCounts,
        alpaka::concepts::IBuffer<std::uint32_t> auto const& droppedRays,
        alpaka::concepts::IBuffer<float> auto& volumePhiAse,
        alpaka::concepts::IBuffer<double> auto& standardError,
        alpaka::concepts::IBuffer<double> auto& relativeStandardError,
        alpaka::concepts::IBuffer<double> auto& volumeDndtAse,
        unsigned rayCount,
        unsigned numIndependentRayPopulations,
        double sourceStrengthTotal)
    {
        auto cellFrameSpec = hase::alpakaUtils::getFrameSpec<uint32_t>(
            devBundle.device,
            devBundle.executor,
            alpaka::Vec{mesh.numberOfCells});
        queue.enqueue(
            cellFrameSpec,
            alpaka::KernelBundle{
                FinalizeForwardVolumePhiAse{rayCount, numIndependentRayPopulations, sourceStrengthTotal},
                mesh,
                vertexPopulationScoreSum,
                rayPopulationRayCounts,
                droppedRays,
                volumePhiAse,
                standardError,
                relativeStandardError,
                volumeDndtAse});
    }

} // namespace hase::kernels
