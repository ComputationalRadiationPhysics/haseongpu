/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/mesh.hpp>
#include <kernels/calcSampleGainSum.hpp>
#include <kernels/forward/rayWalk.hpp>
#include <kernels/forward/volumeSampling.hpp>

#include <cassert>
#include <limits>

namespace hase::kernels::forward
{
    struct AccumulateForwardPhiAse
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            unsigned const forwardRayCount,
            double const forwardRayLength,
            double const betaVolumeTotal,
            alpaka::concepts::IMdSpan auto phiAccumulator,
            alpaka::concepts::IMdSpan auto phiSquareAccumulator,
            alpaka::concepts::IMdSpan auto volumeRayVisits,
            alpaka::concepts::IMdSpan auto droppedRays,
            alpaka::concepts::IMdSpan auto const sigmaA,
            alpaka::concepts::IMdSpan auto const sigmaE,
            unsigned const lambdaResolution,
            unsigned const threadLocalStridingRNG) const
        {
            auto const tIdx = hase::alpakaUtils::getLinGlobalIdx(acc);
            auto rndEngine = alpaka::rand::engine::Philox4x32x10{threadLocalStridingRNG + tIdx};
            for(auto rayNumber : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                (void) rayNumber;
                unsigned tet = sampleVolumeByBetaVolume(mesh, betaVolumeTotal, rndEngine);
                double const sourceWeight = betaVolumeTotal > 0.0 ? 1.0 : 0.0;
                hase::core::Point origin = samplePointInVolume(mesh, tet, rndEngine);
                hase::core::Point const direction = sampleIsotropicDirection(rndEngine);
                unsigned const sigmaIndex = GenRndSigmas{}(lambdaResolution, rndEngine);
                walkRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    forwardRayLength,
                    sourceWeight,
                    sigmaA[sigmaIndex],
                    sigmaE[sigmaIndex],
                    phiAccumulator,
                    phiSquareAccumulator,
                    volumeRayVisits,
                    droppedRays);
            }
        }

        ALPAKA_FN_HOST_ACC void walkRay(
            auto const& acc,
            hase::core::DeviceMeshView const& mesh,
            unsigned tet,
            hase::core::Point origin,
            hase::core::Point const direction,
            double remaining,
            double const sourceWeight,
            double const sigmaA,
            double const sigmaE,
            alpaka::concepts::IMdSpan auto phiAccumulator,
            alpaka::concepts::IMdSpan auto phiSquareAccumulator,
            alpaka::concepts::IMdSpan auto volumeRayVisits,
            alpaka::concepts::IMdSpan auto droppedRays) const
        {
            int forbiddenFace = -1;
            double accumulatedGain = 1.0;
            constexpr double nudgeFactor = 64.0 * std::numeric_limits<double>::epsilon();
            while(remaining > SMALL)
            {
                assert(tet < mesh.numberOfCells);
                double segmentLength = remaining;
                int const nextFace = nextFaceIntersection(mesh, tet, origin, direction, forbiddenFace, segmentLength);
                double const segmentGain = localSegmentGain(mesh, tet, segmentLength, sigmaA, sigmaE);
                double contribution = sourceWeight * accumulatedGain;
                contribution *= localSegmentTrackLengthIntegral(mesh, tet, segmentLength, sigmaA, sigmaE);
                if(alpaka::math::isfinite(contribution))
                {
                    alpaka::onAcc::atomicAdd(acc, &phiAccumulator[tet], contribution);
                    alpaka::onAcc::atomicAdd(acc, &phiSquareAccumulator[tet], contribution * contribution);
                    alpaka::onAcc::atomicAdd(acc, &volumeRayVisits[tet], 1u);
                }
                else
                {
                    alpaka::onAcc::atomicAdd(acc, &droppedRays[tet], 1u);
                }

                accumulatedGain *= segmentGain;
                origin = advance(origin, direction, segmentLength);
                remaining -= segmentLength;
                if(nextFace < 0)
                {
                    break;
                }
                int const neighbor = mesh.getCellNeighbor(tet, static_cast<unsigned>(nextFace));
                if(neighbor < 0)
                {
                    break;
                }
                forbiddenFace = mesh.getCellNeighborLocalFace(tet, static_cast<unsigned>(nextFace));
                tet = static_cast<unsigned>(neighbor);
                double const nudge = nudgeFactor * forwardLengthScale(remaining, segmentLength);
                if(nudge > 0.0 && remaining > nudge)
                {
                    origin = advance(origin, direction, nudge);
                    remaining -= nudge;
                }
            }
        }

        [[nodiscard]] ALPAKA_FN_HOST_ACC double forwardLengthScale(double const remaining, double const segmentLength)
            const
        {
            return alpaka::math::max(remaining, segmentLength);
        }
    };

    [[nodiscard]] ALPAKA_FN_HOST_ACC double boundaryReflectance(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        unsigned const localFace,
        hase::core::Point const direction,
        hase::core::Point const outwardNormal)
    {
        int const boundary = mesh.cellFaceBoundaries[tet * mesh.numberOfFacesPerCell + localFace];
        if(boundary <= 0)
        {
            return 0.0;
        }
        double const nInside = static_cast<double>(mesh.getSurfaceRefractiveIndexInside(tet, localFace));
        double const nOutside = static_cast<double>(mesh.getSurfaceRefractiveIndexOutside(tet, localFace));
        if(nInside > 0.0 && nOutside > 0.0)
        {
            double const cosIncident = alpaka::math::min(
                1.0,
                alpaka::math::abs(hase::core::dot(normalize(direction), outwardNormal)));
            double const sin2Incident = alpaka::math::max(0.0, 1.0 - cosIncident * cosIncident);
            double const ratio = nInside / nOutside;
            if(ratio * ratio * sin2Incident > 1.0)
            {
                return 1.0;
            }
        }
        return alpaka::math::max(0.0, static_cast<double>(mesh.getSurfaceReflectivity(tet, localFace)));
    }

    struct AccumulateForwardPhiAseReservoir
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            unsigned const forwardRayCount,
            double const forwardRayLength,
            double const betaVolumeTotal,
            alpaka::concepts::IMdSpan auto phiAccumulator,
            alpaka::concepts::IMdSpan auto phiSquareAccumulator,
            alpaka::concepts::IMdSpan auto volumeRayVisits,
            alpaka::concepts::IMdSpan auto droppedRays,
            alpaka::concepts::IMdSpan auto reservoirCounts,
            alpaka::concepts::IMdSpan auto reservoirDirX,
            alpaka::concepts::IMdSpan auto reservoirDirY,
            alpaka::concepts::IMdSpan auto reservoirDirZ,
            alpaka::concepts::IMdSpan auto reservoirWeights,
            alpaka::concepts::IMdSpan auto reservoirSigmaIndices,
            alpaka::concepts::IMdSpan auto reservoirTotalWeight,
            unsigned const surfaceReservoirSize,
            alpaka::concepts::IMdSpan auto const sigmaA,
            alpaka::concepts::IMdSpan auto const sigmaE,
            unsigned const lambdaResolution,
            unsigned const threadLocalStridingRNG) const
        {
            auto const tIdx = hase::alpakaUtils::getLinGlobalIdx(acc);
            auto rndEngine = alpaka::rand::engine::Philox4x32x10{threadLocalStridingRNG + tIdx};
            for(auto rayNumber : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                (void) rayNumber;
                unsigned tet = sampleVolumeByBetaVolume(mesh, betaVolumeTotal, rndEngine);
                double const sourceWeight = betaVolumeTotal > 0.0 ? 1.0 : 0.0;
                hase::core::Point origin = samplePointInVolume(mesh, tet, rndEngine);
                hase::core::Point const direction = sampleIsotropicDirection(rndEngine);
                unsigned const sigmaIndex = GenRndSigmas{}(lambdaResolution, rndEngine);
                walkRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    forwardRayLength,
                    sourceWeight,
                    sigmaA[sigmaIndex],
                    sigmaE[sigmaIndex],
                    sigmaIndex,
                    phiAccumulator,
                    phiSquareAccumulator,
                    volumeRayVisits,
                    droppedRays,
                    reservoirCounts,
                    reservoirDirX,
                    reservoirDirY,
                    reservoirDirZ,
                    reservoirWeights,
                    reservoirSigmaIndices,
                    reservoirTotalWeight,
                    surfaceReservoirSize);
            }
        }

        ALPAKA_FN_HOST_ACC void walkRay(
            auto const& acc,
            hase::core::DeviceMeshView const& mesh,
            unsigned tet,
            hase::core::Point origin,
            hase::core::Point const direction,
            double remaining,
            double const sourceWeight,
            double const sigmaA,
            double const sigmaE,
            unsigned const sigmaIndex,
            alpaka::concepts::IMdSpan auto phiAccumulator,
            alpaka::concepts::IMdSpan auto phiSquareAccumulator,
            alpaka::concepts::IMdSpan auto volumeRayVisits,
            alpaka::concepts::IMdSpan auto droppedRays,
            alpaka::concepts::IMdSpan auto reservoirCounts,
            alpaka::concepts::IMdSpan auto reservoirDirX,
            alpaka::concepts::IMdSpan auto reservoirDirY,
            alpaka::concepts::IMdSpan auto reservoirDirZ,
            alpaka::concepts::IMdSpan auto reservoirWeights,
            alpaka::concepts::IMdSpan auto reservoirSigmaIndices,
            alpaka::concepts::IMdSpan auto reservoirTotalWeight,
            unsigned const surfaceReservoirSize) const
        {
            int forbiddenFace = -1;
            double accumulatedGain = 1.0;
            while(remaining > SMALL)
            {
                assert(tet < mesh.numberOfCells);
                double segmentLength = remaining;
                int const nextFace = nextFaceIntersection(mesh, tet, origin, direction, forbiddenFace, segmentLength);
                double const segmentGain = localSegmentGain(mesh, tet, segmentLength, sigmaA, sigmaE);
                double contribution = sourceWeight * accumulatedGain;
                contribution *= localSegmentTrackLengthIntegral(mesh, tet, segmentLength, sigmaA, sigmaE);
                if(alpaka::math::isfinite(contribution))
                {
                    alpaka::onAcc::atomicAdd(acc, &phiAccumulator[tet], contribution);
                    alpaka::onAcc::atomicAdd(acc, &phiSquareAccumulator[tet], contribution * contribution);
                    alpaka::onAcc::atomicAdd(acc, &volumeRayVisits[tet], 1u);
                }
                else
                {
                    alpaka::onAcc::atomicAdd(acc, &droppedRays[tet], 1u);
                }

                accumulatedGain *= segmentGain;
                origin = advance(origin, direction, segmentLength);
                remaining -= segmentLength;
                if(nextFace < 0)
                {
                    break;
                }
                int const neighbor = mesh.getCellNeighbor(tet, static_cast<unsigned>(nextFace));
                if(neighbor < 0)
                {
                    depositReflection(
                        acc,
                        mesh,
                        tet,
                        static_cast<unsigned>(nextFace),
                        direction,
                        sourceWeight * accumulatedGain,
                        sigmaIndex,
                        reservoirCounts,
                        reservoirDirX,
                        reservoirDirY,
                        reservoirDirZ,
                        reservoirWeights,
                        reservoirSigmaIndices,
                        reservoirTotalWeight,
                        surfaceReservoirSize);
                    break;
                }
                forbiddenFace = mesh.getCellNeighborLocalFace(tet, static_cast<unsigned>(nextFace));
                tet = static_cast<unsigned>(neighbor);
                constexpr double nudgeFactor = 64.0 * std::numeric_limits<double>::epsilon();
                double const nudge = nudgeFactor * alpaka::math::max(remaining, segmentLength);
                if(nudge > 0.0 && remaining > nudge)
                {
                    origin = advance(origin, direction, nudge);
                    remaining -= nudge;
                }
            }
        }

        ALPAKA_FN_HOST_ACC void depositReflection(
            auto const& acc,
            hase::core::DeviceMeshView const& mesh,
            unsigned const tet,
            unsigned const localFace,
            hase::core::Point const direction,
            double const incidentWeight,
            unsigned const sigmaIndex,
            alpaka::concepts::IMdSpan auto reservoirCounts,
            alpaka::concepts::IMdSpan auto reservoirDirX,
            alpaka::concepts::IMdSpan auto reservoirDirY,
            alpaka::concepts::IMdSpan auto reservoirDirZ,
            alpaka::concepts::IMdSpan auto reservoirWeights,
            alpaka::concepts::IMdSpan auto reservoirSigmaIndices,
            alpaka::concepts::IMdSpan auto reservoirTotalWeight,
            unsigned const surfaceReservoirSize) const
        {
            if(surfaceReservoirSize == 0u || incidentWeight <= 0.0 || !alpaka::math::isfinite(incidentWeight))
            {
                return;
            }
            hase::core::Point const normal = outwardFaceNormal(mesh, tet, localFace);
            double const reflectance = boundaryReflectance(mesh, tet, localFace, direction, normal);
            double const reflectedWeight = incidentWeight * reflectance;
            if(reflectedWeight <= 0.0 || !alpaka::math::isfinite(reflectedWeight))
            {
                return;
            }
            unsigned const faceId = tet * mesh.numberOfFacesPerCell + localFace;
            unsigned const slot = alpaka::onAcc::atomicAdd(acc, &reservoirCounts[faceId], 1u);
            alpaka::onAcc::atomicAdd(acc, &reservoirTotalWeight[0], reflectedWeight);
            if(slot >= surfaceReservoirSize)
            {
                return;
            }
            unsigned const index = faceId * surfaceReservoirSize + slot;
            hase::core::Point const reflected = reflectedDirection(direction, normal);
            reservoirDirX[index] = reflected.x;
            reservoirDirY[index] = reflected.y;
            reservoirDirZ[index] = reflected.z;
            reservoirWeights[index] = reflectedWeight;
            reservoirSigmaIndices[index] = sigmaIndex;
        }
    };

    struct AccumulateReflectedForwardPhiAse
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            unsigned const totalSlots,
            double const forwardRayLength,
            alpaka::concepts::IMdSpan auto phiAccumulator,
            alpaka::concepts::IMdSpan auto phiSquareAccumulator,
            alpaka::concepts::IMdSpan auto volumeRayVisits,
            alpaka::concepts::IMdSpan auto droppedRays,
            alpaka::concepts::IMdSpan auto const inReservoirCounts,
            alpaka::concepts::IMdSpan auto const inReservoirDirX,
            alpaka::concepts::IMdSpan auto const inReservoirDirY,
            alpaka::concepts::IMdSpan auto const inReservoirDirZ,
            alpaka::concepts::IMdSpan auto const inReservoirWeights,
            alpaka::concepts::IMdSpan auto const inReservoirSigmaIndices,
            alpaka::concepts::IMdSpan auto outReservoirCounts,
            alpaka::concepts::IMdSpan auto outReservoirDirX,
            alpaka::concepts::IMdSpan auto outReservoirDirY,
            alpaka::concepts::IMdSpan auto outReservoirDirZ,
            alpaka::concepts::IMdSpan auto outReservoirWeights,
            alpaka::concepts::IMdSpan auto outReservoirSigmaIndices,
            alpaka::concepts::IMdSpan auto outReservoirTotalWeight,
            unsigned const surfaceReservoirSize,
            alpaka::concepts::IMdSpan auto const sigmaA,
            alpaka::concepts::IMdSpan auto const sigmaE) const
        {
            AccumulateForwardPhiAseReservoir walker;
            for(auto [slotIndex] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{totalSlots}))
            {
                if(surfaceReservoirSize == 0u)
                {
                    continue;
                }
                unsigned const faceId = slotIndex / surfaceReservoirSize;
                unsigned const localSlot = slotIndex - faceId * surfaceReservoirSize;
                unsigned const count = inReservoirCounts[faceId];
                if(localSlot >= count || localSlot >= surfaceReservoirSize)
                {
                    continue;
                }
                unsigned const tet = faceId / mesh.numberOfFacesPerCell;
                unsigned const localFace = faceId - tet * mesh.numberOfFacesPerCell;
                hase::core::Point const direction = normalize(hase::core::Point{
                    inReservoirDirX[slotIndex],
                    inReservoirDirY[slotIndex],
                    inReservoirDirZ[slotIndex]});
                hase::core::Point origin = faceCentroid(mesh, tet, localFace);
                double const nudge = 64.0 * std::numeric_limits<double>::epsilon() * forwardRayLength;
                origin = advance(origin, direction, alpaka::math::max(nudge, 1.0e-12));
                unsigned const sigmaIndex = inReservoirSigmaIndices[slotIndex];
                walker.walkRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    forwardRayLength,
                    inReservoirWeights[slotIndex],
                    sigmaA[sigmaIndex],
                    sigmaE[sigmaIndex],
                    sigmaIndex,
                    phiAccumulator,
                    phiSquareAccumulator,
                    volumeRayVisits,
                    droppedRays,
                    outReservoirCounts,
                    outReservoirDirX,
                    outReservoirDirY,
                    outReservoirDirZ,
                    outReservoirWeights,
                    outReservoirSigmaIndices,
                    outReservoirTotalWeight,
                    surfaceReservoirSize);
            }
        }
    };
} // namespace hase::kernels::forward
