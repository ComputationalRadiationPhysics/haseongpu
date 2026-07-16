/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/mesh.hpp>
#include <kernels/forward/rayTransition.hpp>
#include <kernels/importanceSampling.hpp>
#include <kernels/propagateRay.hpp>
#include <kernels/reflection.hpp>

#include <cassert>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace hase::kernels::forward
{
    ALPAKA_FN_HOST_ACC constexpr std::uint64_t rayHistoryId(unsigned const pass, unsigned const rayIndex)
    {
        return (static_cast<std::uint64_t>(pass) << 32u) | rayIndex;
    }

    ALPAKA_FN_HOST_ACC constexpr std::uint64_t surfaceSamplingHistoryId(unsigned const pass)
    {
        return (std::uint64_t{1u} << 63u) | pass;
    }

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
        alpaka::concepts::IMdSpan TFaceWeights>
    struct SurfaceReservoirSpans
    {
        TCounts counts;
        TDirX dirX;
        TDirY dirY;
        TDirZ dirZ;
        TWeights weights;
        TSigmaIndices sigmaIndices;
        TFaceWeights faceWeights;
        unsigned slotsPerFace;
    };

    template<
        alpaka::concepts::IMdSpan TCdf,
        alpaka::concepts::IMdSpan TTotalWeight,
        alpaka::concepts::IMdSpan TRayFaces>
    struct SurfaceReservoirSamplingCdfSpans
    {
        TCdf cdf;
        TTotalWeight totalWeight;
        TRayFaces rayFaces;
        bool useFaceStratification;
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
        alpaka::concepts::IMdSpan TFaceWeights>
    struct MakeAccessibleOnAcc::Op<
        hase::kernels::forward::
            SurfaceReservoirSpans<TCounts, TDirX, TDirY, TDirZ, TWeights, TSigmaIndices, TFaceWeights>>
    {
        auto operator()(hase::kernels::forward::
                            SurfaceReservoirSpans<TCounts, TDirX, TDirY, TDirZ, TWeights, TSigmaIndices, TFaceWeights>&
                                spans) const
        {
            return hase::kernels::forward::SurfaceReservoirSpans{
                makeAccessibleOnAcc(spans.counts),
                makeAccessibleOnAcc(spans.dirX),
                makeAccessibleOnAcc(spans.dirY),
                makeAccessibleOnAcc(spans.dirZ),
                makeAccessibleOnAcc(spans.weights),
                makeAccessibleOnAcc(spans.sigmaIndices),
                makeAccessibleOnAcc(spans.faceWeights),
                spans.slotsPerFace};
        }

        auto operator()(
            hase::kernels::forward::
                SurfaceReservoirSpans<TCounts, TDirX, TDirY, TDirZ, TWeights, TSigmaIndices, TFaceWeights> const&
                    spans) const
        {
            return hase::kernels::forward::SurfaceReservoirSpans{
                makeAccessibleOnAcc(spans.counts),
                makeAccessibleOnAcc(spans.dirX),
                makeAccessibleOnAcc(spans.dirY),
                makeAccessibleOnAcc(spans.dirZ),
                makeAccessibleOnAcc(spans.weights),
                makeAccessibleOnAcc(spans.sigmaIndices),
                makeAccessibleOnAcc(spans.faceWeights),
                spans.slotsPerFace};
        }
    };

    template<
        alpaka::concepts::IMdSpan TCdf,
        alpaka::concepts::IMdSpan TTotalWeight,
        alpaka::concepts::IMdSpan TRayFaces>
    struct MakeAccessibleOnAcc::Op<
        hase::kernels::forward::SurfaceReservoirSamplingCdfSpans<TCdf, TTotalWeight, TRayFaces>>
    {
        auto operator()(
            hase::kernels::forward::SurfaceReservoirSamplingCdfSpans<TCdf, TTotalWeight, TRayFaces>& spans) const
        {
            return hase::kernels::forward::SurfaceReservoirSamplingCdfSpans{
                makeAccessibleOnAcc(spans.cdf),
                makeAccessibleOnAcc(spans.totalWeight),
                makeAccessibleOnAcc(spans.rayFaces),
                spans.useFaceStratification};
        }

        auto operator()(
            hase::kernels::forward::SurfaceReservoirSamplingCdfSpans<TCdf, TTotalWeight, TRayFaces> const& spans) const
        {
            return hase::kernels::forward::SurfaceReservoirSamplingCdfSpans{
                makeAccessibleOnAcc(spans.cdf),
                makeAccessibleOnAcc(spans.totalWeight),
                makeAccessibleOnAcc(spans.rayFaces),
                spans.useFaceStratification};
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
        alpaka::concepts::IMdSpan TFaceWeights>
    struct IsKernelArgumentTriviallyCopyable<
        hase::kernels::forward::
            SurfaceReservoirSpans<TCounts, TDirX, TDirY, TDirZ, TWeights, TSigmaIndices, TFaceWeights>>
        : std::bool_constant<
              IsKernelArgumentTriviallyCopyable<TCounts>::value && IsKernelArgumentTriviallyCopyable<TDirX>::value
              && IsKernelArgumentTriviallyCopyable<TDirY>::value && IsKernelArgumentTriviallyCopyable<TDirZ>::value
              && IsKernelArgumentTriviallyCopyable<TWeights>::value
              && IsKernelArgumentTriviallyCopyable<TSigmaIndices>::value
              && IsKernelArgumentTriviallyCopyable<TFaceWeights>::value>
    {
    };

    template<
        alpaka::concepts::IMdSpan TCdf,
        alpaka::concepts::IMdSpan TTotalWeight,
        alpaka::concepts::IMdSpan TRayFaces>
    struct IsKernelArgumentTriviallyCopyable<
        hase::kernels::forward::SurfaceReservoirSamplingCdfSpans<TCdf, TTotalWeight, TRayFaces>>
        : std::bool_constant<
              IsKernelArgumentTriviallyCopyable<TCdf>::value && IsKernelArgumentTriviallyCopyable<TTotalWeight>::value
              && IsKernelArgumentTriviallyCopyable<TRayFaces>::value>
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
            unsigned const globalRayOffset,
            unsigned const globalRayCount,
            double const sourceStratificationOffset,
            unsigned const spectrumStratificationPhase,
            double const betaVolumeTotal,
            auto accumulation,
            auto spectrum,
            unsigned const rngSeed) const
        {
            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                unsigned const globalRayIndex = globalRayOffset + rayNumber;
                auto rndEngine = alpaka::rand::engine::Philox4x32x10{rngSeed, rayHistoryId(0u, globalRayIndex)};
                unsigned const tet = sampleStratifiedVolumeByBetaVolume(
                    mesh,
                    betaVolumeTotal,
                    globalRayIndex,
                    globalRayCount,
                    sourceStratificationOffset,
                    rndEngine);
                // Beta-volume sampling carries its importance factor in the
                // source probability, so the compensating weight is one.
                double const sourceWeight = betaVolumeTotal > 0.0 ? 1.0 : 0.0;
                hase::core::Point origin = samplePointInVolume(mesh, tet, rndEngine);
                hase::core::Point const direction = sampleIsotropicDirection(rndEngine);
                unsigned const sigmaIndex = stratifiedSpectrumIndex(
                    spectrum.lambdaResolution,
                    globalRayIndex,
                    globalRayCount,
                    spectrumStratificationPhase);
                walkVolumeSeededForwardRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    sourceWeight,
                    spectrum.sigmaA[sigmaIndex],
                    spectrum.sigmaE[sigmaIndex],
                    accumulation);
            }
        }

        ALPAKA_FN_HOST_ACC void walkVolumeSeededForwardRay(
            auto const& acc,
            hase::core::DeviceMeshView const& mesh,
            unsigned tet,
            hase::core::Point origin,
            hase::core::Point const direction,
            double const sourceWeight,
            double const sigmaAbsorption,
            double const sigmaEmission,
            auto accumulation) const
        {
            int forbiddenFace = -1;
            double accumulatedGain = 1.0;
            constexpr unsigned maxTraversalSteps = 10000u;
            for(unsigned step = 0u; step < maxTraversalSteps; ++step)
            {
                assert(tet < mesh.numberOfCells);
                Tet4FaceIntersection const intersection
                    = nextFaceIntersection(mesh, tet, origin, direction, forbiddenFace);
                if(intersection.localFace < 0)
                {
                    int const recoveryFace = isNearTet4Face(mesh, tet, origin)
                                                 ? immediateExitFace(mesh, tet, origin, direction, forbiddenFace)
                                                 : -1;
                    if(recoveryFace >= 0)
                    {
                        Tet4FaceTransition const transition
                            = recoverFaceTransition(mesh, tet, recoveryFace, origin, direction);
                        if(transition.status == Tet4TransitionStatus::failed)
                        {
                            alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[transition.cell], 1u);
                            break;
                        }
                        if(transition.status == Tet4TransitionStatus::reachedBoundary)
                        {
                            break;
                        }
                        tet = transition.cell;
                        forbiddenFace = transition.forbiddenFace;
                        continue;
                    }
                    alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[tet], 1u);
                    break;
                }
                double const segmentLength = intersection.length;
                double const segmentGain = localSegmentGain(mesh, tet, segmentLength, sigmaAbsorption, sigmaEmission);
                double contribution = sourceWeight * accumulatedGain;
                contribution
                    *= localSegmentTrackLengthIntegral(mesh, tet, segmentLength, sigmaAbsorption, sigmaEmission);
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
                Tet4FaceTransition const transition
                    = transitionAcrossIntersection(mesh, tet, intersection, origin, direction);
                if(transition.status == Tet4TransitionStatus::failed)
                {
                    alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[transition.cell], 1u);
                    break;
                }
                if(transition.status == Tet4TransitionStatus::reachedBoundary)
                {
                    break;
                }
                tet = transition.cell;
                forbiddenFace = transition.forbiddenFace;
                if(step + 1u == maxTraversalSteps)
                    alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[tet], 1u);
            }
        }
    };

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
            auto reservoir,
            auto& rndEngine) const
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
            unsigned const seen = alpaka::onAcc::atomicAdd(acc, &reservoir.counts[faceId], 1u);
            double const priorWeight = alpaka::onAcc::atomicAdd(acc, &reservoir.faceWeights[faceId], reflectedWeight);
            unsigned slot = seen;
            if(slot >= reservoir.slotsPerFace)
            {
                double const replacementProbability = alpaka::math::min(
                    1.0,
                    static_cast<double>(reservoir.slotsPerFace) * reflectedWeight / (priorWeight + reflectedWeight));
                double const selection = alpaka::rand::distribution::UniformReal<double>{}(rndEngine);
                if(selection >= replacementProbability)
                {
                    return;
                }
                slot = static_cast<unsigned>(
                    alpaka::rand::distribution::UniformReal<double>{}(rndEngine)
                    * static_cast<double>(reservoir.slotsPerFace));
                if(slot >= reservoir.slotsPerFace)
                    slot = reservoir.slotsPerFace - 1u;
            }
            unsigned const index = faceId * reservoir.slotsPerFace + slot;
            hase::core::Point const reflected = reflectedDirection(direction, normal);
            reservoir.dirX[index] = reflected.x;
            reservoir.dirY[index] = reflected.y;
            reservoir.dirZ[index] = reflected.z;
            reservoir.weights[index] = reflectedWeight;
            reservoir.sigmaIndices[index] = sigmaIndex;
        }

        ALPAKA_FN_HOST_ACC void walkForwardRay(
            auto const& acc,
            hase::core::DeviceMeshView const& mesh,
            unsigned tet,
            hase::core::Point origin,
            hase::core::Point const direction,
            int const initialForbiddenFace,
            double const sourceWeight,
            double const sigmaAbsorption,
            double const sigmaEmission,
            unsigned const sigmaIndex,
            auto accumulation,
            auto reservoir,
            auto& rndEngine) const
        {
            int forbiddenFace = initialForbiddenFace;
            double accumulatedGain = 1.0;
            constexpr unsigned maxTraversalSteps = 10000u;
            for(unsigned step = 0u; step < maxTraversalSteps; ++step)
            {
                assert(tet < mesh.numberOfCells);
                Tet4FaceIntersection const intersection
                    = nextFaceIntersection(mesh, tet, origin, direction, forbiddenFace);
                if(intersection.localFace < 0)
                {
                    int const recoveryFace = isNearTet4Face(mesh, tet, origin)
                                                 ? immediateExitFace(mesh, tet, origin, direction, forbiddenFace)
                                                 : -1;
                    if(recoveryFace >= 0)
                    {
                        Tet4FaceTransition const transition
                            = recoverFaceTransition(mesh, tet, recoveryFace, origin, direction);
                        if(transition.status == Tet4TransitionStatus::failed)
                        {
                            alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[transition.cell], 1u);
                            break;
                        }
                        if(transition.status == Tet4TransitionStatus::reachedBoundary)
                        {
                            depositReflection(
                                acc,
                                mesh,
                                transition.cell,
                                static_cast<unsigned>(transition.boundaryFace),
                                direction,
                                sourceWeight * accumulatedGain,
                                sigmaIndex,
                                reservoir,
                                rndEngine);
                            break;
                        }
                        tet = transition.cell;
                        forbiddenFace = transition.forbiddenFace;
                        continue;
                    }
                    alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[tet], 1u);
                    break;
                }
                double const segmentLength = intersection.length;
                double const segmentGain = localSegmentGain(mesh, tet, segmentLength, sigmaAbsorption, sigmaEmission);
                double contribution = sourceWeight * accumulatedGain;
                contribution
                    *= localSegmentTrackLengthIntegral(mesh, tet, segmentLength, sigmaAbsorption, sigmaEmission);
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
                Tet4FaceTransition const transition
                    = transitionAcrossIntersection(mesh, tet, intersection, origin, direction);
                if(transition.status == Tet4TransitionStatus::failed)
                {
                    alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[transition.cell], 1u);
                    break;
                }
                if(transition.status == Tet4TransitionStatus::reachedBoundary)
                {
                    depositReflection(
                        acc,
                        mesh,
                        transition.cell,
                        static_cast<unsigned>(transition.boundaryFace),
                        direction,
                        sourceWeight * accumulatedGain,
                        sigmaIndex,
                        reservoir,
                        rndEngine);
                    break;
                }
                tet = transition.cell;
                forbiddenFace = transition.forbiddenFace;
                if(step + 1u == maxTraversalSteps)
                    alpaka::onAcc::atomicAdd(acc, &accumulation.droppedRays[tet], 1u);
            }
        }

        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            unsigned const forwardRayCount,
            unsigned const globalRayOffset,
            unsigned const globalRayCount,
            double const sourceStratificationOffset,
            unsigned const spectrumStratificationPhase,
            double const betaVolumeTotal,
            auto accumulation,
            auto reservoir,
            auto spectrum,
            unsigned const rngSeed) const
        {
            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                unsigned const globalRayIndex = globalRayOffset + rayNumber;
                auto rndEngine = alpaka::rand::engine::Philox4x32x10{rngSeed, rayHistoryId(0u, globalRayIndex)};
                unsigned const tet = sampleStratifiedVolumeByBetaVolume(
                    mesh,
                    betaVolumeTotal,
                    globalRayIndex,
                    globalRayCount,
                    sourceStratificationOffset,
                    rndEngine);
                double const sourceWeight = betaVolumeTotal > 0.0 ? 1.0 : 0.0;
                hase::core::Point origin = samplePointInVolume(mesh, tet, rndEngine);
                hase::core::Point const direction = sampleIsotropicDirection(rndEngine);
                unsigned const sampledSigmaIndex = stratifiedSpectrumIndex(
                    spectrum.lambdaResolution,
                    globalRayIndex,
                    globalRayCount,
                    spectrumStratificationPhase);
                walkForwardRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    -1,
                    sourceWeight,
                    spectrum.sigmaA[sampledSigmaIndex],
                    spectrum.sigmaE[sampledSigmaIndex],
                    sampledSigmaIndex,
                    accumulation,
                    reservoir,
                    rndEngine);
            }
        }
    };

    struct NormalizeSurfaceReservoirSamplingCdf
    {
        ALPAKA_FN_HOST_ACC void operator()(auto const& acc, unsigned const faceCount, auto samplingCdf) const
        {
            double const totalWeight = samplingCdf.totalWeight[0u];
            for(auto [face] :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{faceCount}))
            {
                samplingCdf.cdf[face] = totalWeight > 0.0 ? samplingCdf.cdf[face] / totalWeight : 0.0;
            }
        }
    };

    struct CaptureSurfaceReservoirSamplingTotalWeight
    {
        ALPAKA_FN_HOST_ACC void operator()(auto const&, unsigned const faceCount, auto samplingCdf) const
        {
            samplingCdf.totalWeight[0u] = faceCount == 0u ? 0.0 : samplingCdf.cdf[faceCount - 1u];
        }
    };

    struct GenerateSurfaceReservoirSystematicOffset
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const&,
            auto systematicOffset,
            unsigned const rngSeed,
            unsigned const pass) const
        {
            auto rng = alpaka::rand::engine::Philox4x32x10{rngSeed, surfaceSamplingHistoryId(pass)};
            systematicOffset[0u] = alpaka::rand::distribution::UniformReal<double, alpaka::rand::interval::OO>{}(rng);
        }
    };

    struct AssignSurfaceReservoirStratifiedRayCounts
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            unsigned const faceCount,
            unsigned const rayCount,
            auto samplingCdf,
            auto systematicOffset,
            auto rayCounts) const
        {
            double const offset = systematicOffset[0u];
            for(auto [face] :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{faceCount}))
            {
                double const lowerCdf = face == 0u ? 0.0 : samplingCdf.cdf[face - 1u];
                double const scaledLower = static_cast<double>(rayCount) * lowerCdf - offset;
                double const scaledUpper = static_cast<double>(rayCount) * samplingCdf.cdf[face] - offset;
                rayCounts[face]
                    = static_cast<unsigned>(alpaka::math::floor(scaledUpper) - alpaka::math::floor(scaledLower));
            }
        }
    };

    struct ScatterSurfaceReservoirStratifiedRayFaces
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            unsigned const faceCount,
            auto rayCounts,
            auto rayOffsets,
            auto rayFaces) const
        {
            for(auto [face] :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{faceCount}))
            {
                unsigned const firstRay = rayOffsets[face];
                unsigned const endRay = firstRay + rayCounts[face];
                for(unsigned ray = firstRay; ray < endRay; ++ray)
                {
                    rayFaces[ray] = face;
                }
            }
        }
    };

    struct AccumulateReflectedForwardPhiAse
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            unsigned const forwardRayCount,
            double const sourceWeight,
            auto accumulation,
            auto inReservoir,
            auto samplingCdf,
            auto outReservoir,
            auto spectrum,
            unsigned const rngSeed,
            unsigned const reflectionPass) const
        {
            AccumulateForwardPhiAseReservoir walker;
            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                auto rndEngine
                    = alpaka::rand::engine::Philox4x32x10{rngSeed, rayHistoryId(reflectionPass, rayNumber)};
                unsigned const faceCount = mesh.numberOfCells * mesh.numberOfFacesPerCell;
                if(inReservoir.slotsPerFace == 0u || faceCount == 0u || samplingCdf.totalWeight[0u] <= 0.0)
                {
                    continue;
                }

                unsigned faceId = 0u;
                if(samplingCdf.useFaceStratification)
                {
                    faceId = samplingCdf.rayFaces[rayNumber];
                }
                else
                {
                    double const faceTarget = alpaka::rand::distribution::UniformReal<double>{}(rndEngine);
                    unsigned lower = 0u;
                    unsigned upper = faceCount;
                    while(lower < upper)
                    {
                        unsigned const middle = lower + (upper - lower) / 2u;
                        if(samplingCdf.cdf[middle] <= faceTarget)
                        {
                            lower = middle + 1u;
                        }
                        else
                        {
                            upper = middle;
                        }
                    }
                    faceId = lower < faceCount ? lower : faceCount - 1u;
                }
                unsigned const filledSlots = alpaka::math::min(inReservoir.counts[faceId], inReservoir.slotsPerFace);
                if(filledSlots == 0u)
                    continue;
                double slotWeight = 0.0;
                unsigned const offset = faceId * inReservoir.slotsPerFace;
                for(unsigned slot = 0u; slot < filledSlots; ++slot)
                    slotWeight += inReservoir.weights[offset + slot];
                if(slotWeight <= 0.0)
                    continue;
                double const slotTarget = alpaka::rand::distribution::UniformReal<double>{}(rndEngine) *slotWeight;
                double cumulativeSlotWeight = 0.0;
                unsigned localSlot = filledSlots - 1u;
                for(unsigned slot = 0u; slot < filledSlots; ++slot)
                {
                    cumulativeSlotWeight += inReservoir.weights[offset + slot];
                    if(cumulativeSlotWeight >= slotTarget)
                    {
                        localSlot = slot;
                        break;
                    }
                }
                unsigned const slotIndex = offset + localSlot;
                unsigned const tet = faceId / mesh.numberOfFacesPerCell;
                unsigned const localFace = faceId - tet * mesh.numberOfFacesPerCell;
                hase::core::Point const direction = normalize(
                    hase::core::Point{
                        inReservoir.dirX[slotIndex],
                        inReservoir.dirY[slotIndex],
                        inReservoir.dirZ[slotIndex]});
                hase::core::Point const origin = faceCentroid(mesh, tet, localFace);
                unsigned const sigmaIndex = inReservoir.sigmaIndices[slotIndex];
                walker.walkForwardRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    static_cast<int>(localFace),
                    sourceWeight,
                    spectrum.sigmaA[sigmaIndex],
                    spectrum.sigmaE[sigmaIndex],
                    sigmaIndex,
                    accumulation,
                    outReservoir,
                    rndEngine);
            }
        }
    };
} // namespace hase::kernels::forward
