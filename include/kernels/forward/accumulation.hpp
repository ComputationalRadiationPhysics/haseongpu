/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/mesh.hpp>
#include <kernels/forward/barycentric.hpp>
#include <kernels/forward/policyRay.hpp>
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
    inline constexpr unsigned forwardRseBatchCount = 8u;

    ALPAKA_FN_HOST_ACC constexpr unsigned rseBatchRayCount(
        unsigned const globalRayOffset,
        unsigned const rayCount,
        unsigned const batch)
    {
        unsigned const end = globalRayOffset + rayCount;
        unsigned const first
            = globalRayOffset + (batch + forwardRseBatchCount - globalRayOffset % forwardRseBatchCount)
                                    % forwardRseBatchCount;
        return first < end ? 1u + (end - 1u - first) / forwardRseBatchCount : 0u;
    }

    ALPAKA_FN_HOST_ACC constexpr std::uint64_t mixRseBatchSeed64(std::uint64_t value)
    {
        std::uint64_t const multiplier = 0xe9'846a'fb1a'615dull;
        value ^= value >> 32u;
        value *= multiplier;
        value ^= value >> 32u;
        value *= multiplier;
        value ^= value >> 28u;
        return value;
    }

    ALPAKA_FN_HOST_ACC constexpr unsigned rseBatchSeed(unsigned const applicationSeed, unsigned const batch)
    {
        return static_cast<unsigned>(mixRseBatchSeed64(
            static_cast<std::uint64_t>(applicationSeed) + 0x9e37'79b9ull + 0x85eb'ca6bull * (batch + 1u)));
    }

    ALPAKA_FN_HOST_ACC constexpr double rseBatchSourceStratificationOffset(
        unsigned const applicationSeed,
        unsigned const batch)
    {
        return static_cast<double>(rseBatchSeed(rseBatchSeed(applicationSeed, batch), 0x7d3a'9f21u))
               / 4294967296.0;
    }

    ALPAKA_FN_HOST_ACC constexpr unsigned rseBatchSpectrumStratificationPhase(
        unsigned const applicationSeed,
        unsigned const batch,
        unsigned const spectrumSize)
    {
        return spectrumSize == 0u
                   ? 0u
                   : rseBatchSeed(rseBatchSeed(applicationSeed, batch), 0x6ca4'c37du) % spectrumSize;
    }

    ALPAKA_FN_HOST_ACC constexpr unsigned rseBatchRayIndex(unsigned const globalRayIndex)
    {
        return globalRayIndex / forwardRseBatchCount;
    }

    struct ForwardAseRayState
        : ray::TraversalState
        , ray::SrmPositionStorage<typename std::remove_cvref_t<ALPAKA_TYPEOF(ray::aseSrmPolicy)>::PositionPolicy>
    {
        double weight = 0.0;
        double sigmaAbsorption = 0.0;
        double sigmaEmission = 0.0;
        double accumulatedGain = 1.0;
        unsigned sigmaIndex = 0u;
        unsigned rseBatch = 0u;
    };

    static_assert(!std::derived_from<ForwardAseRayState, ray::BarycentricSrmPositionStorage>);

    ALPAKA_FN_HOST_ACC constexpr std::uint64_t rayHistoryId(unsigned const pass, unsigned const rayIndex)
    {
        return (static_cast<std::uint64_t>(pass) << 32u) | rayIndex;
    }

    ALPAKA_FN_HOST_ACC constexpr std::uint64_t surfaceSamplingHistoryId(unsigned const pass)
    {
        return (std::uint64_t{1u} << 63u) | pass;
    }

    template<
        alpaka::concepts::IMdSpan TVertexBatchScoreSum,
        alpaka::concepts::IMdSpan TCellRayVisits,
        alpaka::concepts::IMdSpan TCellDroppedRays>
    struct ForwardAccumulationSpans
    {
        TVertexBatchScoreSum vertexBatchScoreSum;
        TCellRayVisits cellRayVisits;
        TCellDroppedRays cellDroppedRays;
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

    template<typename T_Accumulation>
    struct ForwardAseCellPolicy : ray::behaviourDimension::Cell
    {
        T_Accumulation accumulation;

        ALPAKA_FN_HOST_ACC constexpr explicit ForwardAseCellPolicy(T_Accumulation value) : accumulation{value}
        {
        }

        ALPAKA_FN_ACC bool operator()(
            auto const& acc,
            hase::core::DeviceMeshView const& mesh,
            auto& rayState,
            unsigned const tet,
            Tet4FaceIntersection const intersection)
        {
            double const segmentGain
                = localSegmentGain(mesh, tet, intersection.length, rayState.sigmaAbsorption, rayState.sigmaEmission);
            double contribution = rayState.weight * rayState.accumulatedGain;
            contribution *= localSegmentTrackLengthIntegral(
                mesh,
                tet,
                intersection.length,
                rayState.sigmaAbsorption,
                rayState.sigmaEmission);
            if(alpaka::math::isfinite(contribution))
            {
                alpaka::onAcc::atomicAdd(acc, &accumulation.cellRayVisits[tet], 1u);
                auto const weights = segmentMidpointBarycentricVertexWeights(
                    mesh,
                    tet,
                    rayState.position,
                    rayState.direction,
                    intersection.length);
                unsigned const materialVertexOffset
                    = mesh.getCellType(tet) == mesh.claddingNumber ? mesh.numberOfMeshPoints : 0u;
                for(unsigned localVertex = 0u; localVertex < hase::core::tet4VertexCount; ++localVertex)
                {
                    unsigned const materialVertex = materialVertexOffset
                                                    + mesh.cellPointIndices[
                                                        tet * mesh.numberOfCellVertices + localVertex];
                    unsigned const vertex
                        = rayState.rseBatch * (2u * mesh.numberOfMeshPoints) + materialVertex;
                    double const weight = weights[localVertex];
                    alpaka::onAcc::atomicAdd(
                        acc,
                        &accumulation.vertexBatchScoreSum[vertex],
                        contribution * weight);
                }
            }
            else
            {
                alpaka::onAcc::atomicAdd(acc, &accumulation.cellDroppedRays[tet], 1u);
            }
            rayState.accumulatedGain *= segmentGain;
            return true;
        }
    };

    template<typename T_Accumulation>
    struct CountDroppedForwardRay
    {
        T_Accumulation accumulation;

        ALPAKA_FN_HOST_ACC constexpr explicit CountDroppedForwardRay(T_Accumulation value) : accumulation{value}
        {
        }

        ALPAKA_FN_ACC void operator()(auto const& acc, auto const&, auto& rayState)
        {
            alpaka::onAcc::atomicAdd(acc, &accumulation.cellDroppedRays[rayState.cell], 1u);
        }
    };

    struct AccumulateForwardPhiAse
    {
        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            unsigned const forwardRayCount,
            unsigned const globalRayOffset,
            unsigned const globalRayCount,
            double const,
            unsigned const,
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
                unsigned const batch = globalRayIndex % forwardRseBatchCount;
                unsigned const batchRayIndex = rseBatchRayIndex(globalRayIndex);
                unsigned const batchRayCount = rseBatchRayCount(0u, globalRayCount, batch);
                unsigned const batchSeed = rseBatchSeed(rngSeed, batch);
                auto rndEngine = alpaka::rand::engine::Philox4x32x10{batchSeed, rayHistoryId(0u, batchRayIndex)};
                unsigned const tet = sampleStratifiedVolumeByBetaVolume(
                    mesh,
                    betaVolumeTotal,
                    batchRayIndex,
                    batchRayCount,
                    rseBatchSourceStratificationOffset(rngSeed, batch),
                    rndEngine);
                // Beta-volume sampling carries its importance factor in the
                // source probability, so the compensating weight is one.
                double const sourceWeight = betaVolumeTotal > 0.0 ? 1.0 : 0.0;
                hase::core::Point origin = samplePointInVolume(mesh, tet, rndEngine);
                hase::core::Point const direction = sampleIsotropicDirection(rndEngine);
                unsigned const sigmaIndex = stratifiedSpectrumIndex(
                    spectrum.lambdaResolution,
                    batchRayIndex,
                    batchRayCount,
                    rseBatchSpectrumStratificationPhase(rngSeed, batch, spectrum.lambdaResolution));
                walkVolumeSeededForwardRay(
                    acc,
                    mesh,
                    tet,
                    origin,
                    direction,
                    sourceWeight,
                    spectrum.sigmaA[sigmaIndex],
                    spectrum.sigmaE[sigmaIndex],
                    batch,
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
            unsigned const rseBatch,
            auto accumulation) const
        {
            namespace ray = hase::kernels::forward::ray;
            ForwardAseRayState rayState;
            rayState.position = origin;
            rayState.direction = direction;
            rayState.cell = tet;
            rayState.forbiddenFace = -1;
            rayState.weight = sourceWeight;
            rayState.sigmaAbsorption = sigmaAbsorption;
            rayState.sigmaEmission = sigmaEmission;
            rayState.rseBatch = rseBatch;
            auto const walkResult = ray::walk(
                acc,
                mesh,
                rayState,
                ray::RayWalkBehaviour{
                    ForwardAseCellPolicy<ALPAKA_TYPEOF(accumulation)>{accumulation},
                    ray::BoundaryPolicyEscape{}});
            if(walkResult == ray::WalkResult::failed)
                CountDroppedForwardRay<ALPAKA_TYPEOF(accumulation)>{accumulation}(acc, mesh, rayState);
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

        template<typename T_Reservoir, typename T_Rng>
        struct StoreReflectionBoundary : ray::BoundaryPolicySrm<ray::srmPosition::Centroid>
        {
            T_Reservoir reservoir;
            T_Rng rng;

            ALPAKA_FN_HOST_ACC constexpr StoreReflectionBoundary(T_Reservoir reservoirValue, T_Rng rngValue)
                : reservoir{reservoirValue}
                , rng{rngValue}
            {
            }

            ALPAKA_FN_ACC ray::BoundaryResult operator()(
                auto const& acc,
                hase::core::DeviceMeshView const& mesh,
                auto& rayState,
                unsigned const tet,
                unsigned const localFace)
            {
                AccumulateForwardPhiAseReservoir{}.depositReflection(
                    acc,
                    mesh,
                    tet,
                    localFace,
                    rayState.direction,
                    rayState.weight * rayState.accumulatedGain,
                    rayState.sigmaIndex,
                    reservoir,
                    rng);
                return ray::BoundaryResult::stop;
            }
        };

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
            unsigned const rseBatch,
            auto accumulation,
            auto reservoir,
            auto& rndEngine) const
        {
            namespace ray = hase::kernels::forward::ray;
            ForwardAseRayState rayState;
            rayState.position = origin;
            rayState.direction = direction;
            rayState.cell = tet;
            rayState.forbiddenFace = initialForbiddenFace;
            rayState.weight = sourceWeight;
            rayState.sigmaAbsorption = sigmaAbsorption;
            rayState.sigmaEmission = sigmaEmission;
            rayState.sigmaIndex = sigmaIndex;
            rayState.rseBatch = rseBatch;
            auto const walkResult = ray::walk(
                acc,
                mesh,
                rayState,
                ray::RayWalkBehaviour{
                    ForwardAseCellPolicy<ALPAKA_TYPEOF(accumulation)>{accumulation},
                    StoreReflectionBoundary<ALPAKA_TYPEOF(reservoir), ALPAKA_TYPEOF(rndEngine)>{
                        reservoir,
                        rndEngine}});
            if(walkResult == ray::WalkResult::failed)
                CountDroppedForwardRay<ALPAKA_TYPEOF(accumulation)>{accumulation}(acc, mesh, rayState);
        }

        ALPAKA_FN_HOST_ACC void operator()(
            auto const& acc,
            hase::core::DeviceMeshView const mesh,
            unsigned const forwardRayCount,
            unsigned const globalRayOffset,
            unsigned const globalRayCount,
            double const,
            unsigned const,
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
                unsigned const batch = globalRayIndex % forwardRseBatchCount;
                unsigned const batchRayIndex = rseBatchRayIndex(globalRayIndex);
                unsigned const batchRayCount = rseBatchRayCount(0u, globalRayCount, batch);
                unsigned const batchSeed = rseBatchSeed(rngSeed, batch);
                auto rndEngine = alpaka::rand::engine::Philox4x32x10{batchSeed, rayHistoryId(0u, batchRayIndex)};
                unsigned const tet = sampleStratifiedVolumeByBetaVolume(
                    mesh,
                    betaVolumeTotal,
                    batchRayIndex,
                    batchRayCount,
                    rseBatchSourceStratificationOffset(rngSeed, batch),
                    rndEngine);
                double const sourceWeight = betaVolumeTotal > 0.0 ? 1.0 : 0.0;
                hase::core::Point origin = samplePointInVolume(mesh, tet, rndEngine);
                hase::core::Point const direction = sampleIsotropicDirection(rndEngine);
                unsigned const sampledSigmaIndex = stratifiedSpectrumIndex(
                    spectrum.lambdaResolution,
                    batchRayIndex,
                    batchRayCount,
                    rseBatchSpectrumStratificationPhase(rngSeed, batch, spectrum.lambdaResolution));
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
                    batch,
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
            unsigned const globalRayOffset,
            unsigned const globalRayCount,
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
                unsigned const globalRayIndex = globalRayOffset + rayNumber;
                unsigned const batch = globalRayIndex % forwardRseBatchCount;
                unsigned const batchRayIndex = rseBatchRayIndex(globalRayIndex);
                unsigned const batchRayCount = rseBatchRayCount(0u, globalRayCount, batch);
                unsigned const batchSeed = rseBatchSeed(rngSeed, batch);
                auto rndEngine
                    = alpaka::rand::engine::Philox4x32x10{batchSeed, rayHistoryId(reflectionPass, batchRayIndex)};
                unsigned const faceCount = mesh.numberOfCells * mesh.numberOfFacesPerCell;
                if(inReservoir.slotsPerFace == 0u || faceCount == 0u || samplingCdf.totalWeight[0u] <= 0.0)
                {
                    continue;
                }

                unsigned faceId = 0u;
                if(samplingCdf.useFaceStratification)
                {
                    double const faceTarget = stratifiedUnitInterval(
                        batchRayIndex,
                        batchRayCount,
                        rseBatchSourceStratificationOffset(rngSeed ^ reflectionPass, batch));
                    unsigned lower = 0u;
                    unsigned upper = faceCount;
                    while(lower < upper)
                    {
                        unsigned const middle = lower + (upper - lower) / 2u;
                        if(samplingCdf.cdf[middle] <= faceTarget)
                            lower = middle + 1u;
                        else
                            upper = middle;
                    }
                    faceId = lower < faceCount ? lower : faceCount - 1u;
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
                hase::core::Point const origin
                    = ray::restoreSrmPosition(ray::aseSrmPolicy.positionPolicy, mesh, tet, localFace);
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
                    batch,
                    accumulation,
                    outReservoir,
                    rndEngine);
            }
        }
    };
} // namespace hase::kernels::forward
