/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/mesh.hpp>
#include <kernels/forward/rayWalk.hpp>
#include <kernels/forward/volumeSampling.hpp>
#include <kernels/calcSampleGainSum.hpp>

#include <cassert>

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

        [[nodiscard]] ALPAKA_FN_HOST_ACC double forwardLengthScale(double const remaining, double const segmentLength) const
        {
            return alpaka::math::max(remaining, segmentLength);
        }
    };
} // namespace hase::kernels::forward
