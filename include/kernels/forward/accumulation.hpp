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
#include <type_traits>

namespace hase::kernels::forward
{
    template<
        alpaka::concepts::IMdSpan TPhi,
        alpaka::concepts::IMdSpan TPhiSquare,
        alpaka::concepts::IMdSpan TVolumeRayVisits,
        alpaka::concepts::IMdSpan TDroppedRays>
    struct ForwardAccumulationSpans
    {
        TPhi phi;
        TPhiSquare phiSquare;
        TVolumeRayVisits volumeRayVisits;
        TDroppedRays droppedRays;
    };

    template<alpaka::concepts::IMdSpan TSigmaA, alpaka::concepts::IMdSpan TSigmaE>
    struct ForwardSpectrumSpans
    {
        TSigmaA sigmaA;
        TSigmaE sigmaE;
        unsigned lambdaResolution;
    };

    template<
        alpaka::concepts::IMdSpan TCounts,
        alpaka::concepts::IMdSpan TDirX,
        alpaka::concepts::IMdSpan TDirY,
        alpaka::concepts::IMdSpan TDirZ,
        alpaka::concepts::IMdSpan TWeights,
        alpaka::concepts::IMdSpan TSigmaIndices,
        alpaka::concepts::IMdSpan TTotalWeight>
    struct SurfaceReservoirSpans
    {
        TCounts counts;
        TDirX dirX;
        TDirY dirY;
        TDirZ dirZ;
        TWeights weights;
        TSigmaIndices sigmaIndices;
        TTotalWeight totalWeight;
        unsigned slotsPerFace;
    };
} // namespace hase::kernels::forward

namespace alpaka::onHost
{
    template<
        alpaka::concepts::IMdSpan TPhi,
        alpaka::concepts::IMdSpan TPhiSquare,
        alpaka::concepts::IMdSpan TVolumeRayVisits,
        alpaka::concepts::IMdSpan TDroppedRays>
    struct MakeAccessibleOnAcc::Op<
        hase::kernels::forward::ForwardAccumulationSpans<TPhi, TPhiSquare, TVolumeRayVisits, TDroppedRays>>
    {
        auto operator()(
            hase::kernels::forward::ForwardAccumulationSpans<TPhi, TPhiSquare, TVolumeRayVisits, TDroppedRays>& spans)
            const
        {
            return hase::kernels::forward::ForwardAccumulationSpans{
                makeAccessibleOnAcc(spans.phi),
                makeAccessibleOnAcc(spans.phiSquare),
                makeAccessibleOnAcc(spans.volumeRayVisits),
                makeAccessibleOnAcc(spans.droppedRays)};
        }

        auto operator()(
            hase::kernels::forward::ForwardAccumulationSpans<TPhi, TPhiSquare, TVolumeRayVisits, TDroppedRays> const&
                spans) const
        {
            return hase::kernels::forward::ForwardAccumulationSpans{
                makeAccessibleOnAcc(spans.phi),
                makeAccessibleOnAcc(spans.phiSquare),
                makeAccessibleOnAcc(spans.volumeRayVisits),
                makeAccessibleOnAcc(spans.droppedRays)};
        }
    };

    template<alpaka::concepts::IMdSpan TSigmaA, alpaka::concepts::IMdSpan TSigmaE>
    struct MakeAccessibleOnAcc::Op<hase::kernels::forward::ForwardSpectrumSpans<TSigmaA, TSigmaE>>
    {
        auto operator()(hase::kernels::forward::ForwardSpectrumSpans<TSigmaA, TSigmaE>& spans) const
        {
            return hase::kernels::forward::ForwardSpectrumSpans{
                makeAccessibleOnAcc(spans.sigmaA),
                makeAccessibleOnAcc(spans.sigmaE),
                spans.lambdaResolution};
        }

        auto operator()(hase::kernels::forward::ForwardSpectrumSpans<TSigmaA, TSigmaE> const& spans) const
        {
            return hase::kernels::forward::ForwardSpectrumSpans{
                makeAccessibleOnAcc(spans.sigmaA),
                makeAccessibleOnAcc(spans.sigmaE),
                spans.lambdaResolution};
        }
    };

    template<
        alpaka::concepts::IMdSpan TCounts,
        alpaka::concepts::IMdSpan TDirX,
        alpaka::concepts::IMdSpan TDirY,
        alpaka::concepts::IMdSpan TDirZ,
        alpaka::concepts::IMdSpan TWeights,
        alpaka::concepts::IMdSpan TSigmaIndices,
        alpaka::concepts::IMdSpan TTotalWeight>
    struct MakeAccessibleOnAcc::Op<
        hase::kernels::forward::
            SurfaceReservoirSpans<TCounts, TDirX, TDirY, TDirZ, TWeights, TSigmaIndices, TTotalWeight>>
    {
        auto operator()(hase::kernels::forward::
                            SurfaceReservoirSpans<TCounts, TDirX, TDirY, TDirZ, TWeights, TSigmaIndices, TTotalWeight>&
                                spans) const
        {
            return hase::kernels::forward::SurfaceReservoirSpans{
                makeAccessibleOnAcc(spans.counts),
                makeAccessibleOnAcc(spans.dirX),
                makeAccessibleOnAcc(spans.dirY),
                makeAccessibleOnAcc(spans.dirZ),
                makeAccessibleOnAcc(spans.weights),
                makeAccessibleOnAcc(spans.sigmaIndices),
                makeAccessibleOnAcc(spans.totalWeight),
                spans.slotsPerFace};
        }

        auto operator()(
            hase::kernels::forward::
                SurfaceReservoirSpans<TCounts, TDirX, TDirY, TDirZ, TWeights, TSigmaIndices, TTotalWeight> const&
                    spans) const
        {
            return hase::kernels::forward::SurfaceReservoirSpans{
                makeAccessibleOnAcc(spans.counts),
                makeAccessibleOnAcc(spans.dirX),
                makeAccessibleOnAcc(spans.dirY),
                makeAccessibleOnAcc(spans.dirZ),
                makeAccessibleOnAcc(spans.weights),
                makeAccessibleOnAcc(spans.sigmaIndices),
                makeAccessibleOnAcc(spans.totalWeight),
                spans.slotsPerFace};
        }
    };
} // namespace alpaka::onHost

namespace alpaka::trait
{
    template<
        alpaka::concepts::IMdSpan TPhi,
        alpaka::concepts::IMdSpan TPhiSquare,
        alpaka::concepts::IMdSpan TVolumeRayVisits,
        alpaka::concepts::IMdSpan TDroppedRays>
    struct IsKernelArgumentTriviallyCopyable<
        hase::kernels::forward::ForwardAccumulationSpans<TPhi, TPhiSquare, TVolumeRayVisits, TDroppedRays>>
        : std::bool_constant<
              IsKernelArgumentTriviallyCopyable<TPhi>::value && IsKernelArgumentTriviallyCopyable<TPhiSquare>::value
              && IsKernelArgumentTriviallyCopyable<TVolumeRayVisits>::value
              && IsKernelArgumentTriviallyCopyable<TDroppedRays>::value>
    {
    };

    template<alpaka::concepts::IMdSpan TSigmaA, alpaka::concepts::IMdSpan TSigmaE>
    struct IsKernelArgumentTriviallyCopyable<hase::kernels::forward::ForwardSpectrumSpans<TSigmaA, TSigmaE>>
        : std::bool_constant<
              IsKernelArgumentTriviallyCopyable<TSigmaA>::value && IsKernelArgumentTriviallyCopyable<TSigmaE>::value>
    {
    };

    template<
        alpaka::concepts::IMdSpan TCounts,
        alpaka::concepts::IMdSpan TDirX,
        alpaka::concepts::IMdSpan TDirY,
        alpaka::concepts::IMdSpan TDirZ,
        alpaka::concepts::IMdSpan TWeights,
        alpaka::concepts::IMdSpan TSigmaIndices,
        alpaka::concepts::IMdSpan TTotalWeight>
    struct IsKernelArgumentTriviallyCopyable<
        hase::kernels::forward::
            SurfaceReservoirSpans<TCounts, TDirX, TDirY, TDirZ, TWeights, TSigmaIndices, TTotalWeight>>
        : std::bool_constant<
              IsKernelArgumentTriviallyCopyable<TCounts>::value && IsKernelArgumentTriviallyCopyable<TDirX>::value
              && IsKernelArgumentTriviallyCopyable<TDirY>::value && IsKernelArgumentTriviallyCopyable<TDirZ>::value
              && IsKernelArgumentTriviallyCopyable<TWeights>::value
              && IsKernelArgumentTriviallyCopyable<TSigmaIndices>::value
              && IsKernelArgumentTriviallyCopyable<TTotalWeight>::value>
    {
    };
} // namespace alpaka::trait

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
            auto accumulation,
            auto spectrum,
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
                unsigned const sigmaIndex = GenRndSigmas{}(spectrum.lambdaResolution, rndEngine);
                walkRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    forwardRayLength,
                    sourceWeight,
                    spectrum.sigmaA[sigmaIndex],
                    spectrum.sigmaE[sigmaIndex],
                    accumulation);
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
            auto accumulation) const
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
                    alpaka::onAcc::atomicAdd(acc, &accumulation.phi[tet], contribution);
                    alpaka::onAcc::atomicAdd(acc, &accumulation.phiSquare[tet], contribution * contribution);
                    alpaka::onAcc::atomicAdd(acc, &accumulation.volumeRayVisits[tet], 1u);
                }
                else
                {
                    alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[tet], 1u);
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

    [[nodiscard]] inline ALPAKA_FN_HOST_ACC double boundaryReflectance(
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
            double const cosIncident
                = alpaka::math::min(1.0, alpaka::math::abs(hase::core::dot(normalize(direction), outwardNormal)));
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
        ALPAKA_FN_HOST_ACC void depositReflection(
            auto const& acc,
            hase::core::DeviceMeshView const& mesh,
            unsigned const tet,
            unsigned const localFace,
            hase::core::Point const direction,
            double const incidentWeight,
            unsigned const sigmaIndex,
            auto reservoir) const
        {
            if(reservoir.slotsPerFace == 0u || incidentWeight <= 0.0 || !alpaka::math::isfinite(incidentWeight))
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
            unsigned const slot = alpaka::onAcc::atomicAdd(acc, &reservoir.counts[faceId], 1u);
            alpaka::onAcc::atomicAdd(acc, &reservoir.totalWeight[0], reflectedWeight);
            if(slot >= reservoir.slotsPerFace)
            {
                return;
            }
            unsigned const index = faceId * reservoir.slotsPerFace + slot;
            hase::core::Point const reflected = reflectedDirection(direction, normal);
            reservoir.dirX[index] = reflected.x;
            reservoir.dirY[index] = reflected.y;
            reservoir.dirZ[index] = reflected.z;
            reservoir.weights[index] = reflectedWeight;
            reservoir.sigmaIndices[index] = sigmaIndex;
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
            auto accumulation,
            auto reservoir) const
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
                    alpaka::onAcc::atomicAdd(acc, &accumulation.phi[tet], contribution);
                    alpaka::onAcc::atomicAdd(acc, &accumulation.phiSquare[tet], contribution * contribution);
                    alpaka::onAcc::atomicAdd(acc, &accumulation.volumeRayVisits[tet], 1u);
                }
                else
                {
                    alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[tet], 1u);
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
                        reservoir);
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

        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            unsigned const forwardRayCount,
            double const forwardRayLength,
            double const betaVolumeTotal,
            auto accumulation,
            auto reservoir,
            auto spectrum,
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
                unsigned const sampledSigmaIndex = GenRndSigmas{}(spectrum.lambdaResolution, rndEngine);
                walkRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    forwardRayLength,
                    sourceWeight,
                    spectrum.sigmaA[sampledSigmaIndex],
                    spectrum.sigmaE[sampledSigmaIndex],
                    sampledSigmaIndex,
                    accumulation,
                    reservoir);
            }
        }
    };

    struct AccumulateReflectedForwardPhiAse
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            unsigned const totalSlots,
            double const forwardRayLength,
            auto accumulation,
            auto inReservoir,
            auto outReservoir,
            auto spectrum) const
        {
            AccumulateForwardPhiAseReservoir walker;
            for(auto [slotIndex] :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{totalSlots}))
            {
                if(inReservoir.slotsPerFace == 0u)
                {
                    continue;
                }
                unsigned const faceId = slotIndex / inReservoir.slotsPerFace;
                unsigned const localSlot = slotIndex - faceId * inReservoir.slotsPerFace;
                unsigned const count = inReservoir.counts[faceId];
                if(localSlot >= count || localSlot >= inReservoir.slotsPerFace)
                {
                    continue;
                }
                unsigned const tet = faceId / mesh.numberOfFacesPerCell;
                unsigned const localFace = faceId - tet * mesh.numberOfFacesPerCell;
                hase::core::Point const direction = normalize(
                    hase::core::Point{
                        inReservoir.dirX[slotIndex],
                        inReservoir.dirY[slotIndex],
                        inReservoir.dirZ[slotIndex]});
                hase::core::Point origin = faceCentroid(mesh, tet, localFace);
                double const nudge = 64.0 * std::numeric_limits<double>::epsilon() * forwardRayLength;
                origin = advance(origin, direction, alpaka::math::max(nudge, 1.0e-12));
                unsigned const sigmaIndex = inReservoir.sigmaIndices[slotIndex];
                walker.walkRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    forwardRayLength,
                    inReservoir.weights[slotIndex],
                    spectrum.sigmaA[sigmaIndex],
                    spectrum.sigmaE[sigmaIndex],
                    sigmaIndex,
                    accumulation,
                    outReservoir);
            }
        }
    };
} // namespace hase::kernels::forward
