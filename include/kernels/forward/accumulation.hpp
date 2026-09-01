/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/geometry.hpp>
#include <data/TraceData.hpp>
#include <kernels/forward/barycentric.hpp>
#include <kernels/forward/policyRay.hpp>
#include <kernels/forward/randomHistory.hpp>
#include <kernels/forward/rayTransition.hpp>
#include <kernels/forward/rayWalk.hpp>
#include <kernels/forward/reflectionResampling.hpp>
#include <kernels/forward/surfaceReservoir.hpp>
#include <kernels/forward/tracePolicyList.hpp>
#include <kernels/forward/volumeSampling.hpp>
#include <kernels/reflection.hpp>
#include <random/randomEngine.hpp>

#include <cassert>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace hase::kernels::forward
{
    inline constexpr unsigned defaultForwardRseBatchCount = 8u;

    /** @brief Resolve the statistical batch count for a worker group. */
    /**
     * @param workerCount Number of workers participating in the trace.
     * @return Statistical batch count, never smaller than the default or worker count.
     */
    ALPAKA_FN_HOST_ACC constexpr unsigned forwardRseBatchCount(unsigned const workerCount)
    {
        return workerCount > defaultForwardRseBatchCount ? workerCount : defaultForwardRseBatchCount;
    }

    /**
     * @param globalRayOffset First global ray index in the launch.
     * @param rayCount Number of launched histories.
     * @param batch Statistical batch index.
     * @param batchCount Total number of interleaved batches.
     * @return Histories in the launch assigned to `batch`.
     */
    ALPAKA_FN_HOST_ACC constexpr unsigned rseBatchRayCount(
        unsigned const globalRayOffset,
        unsigned const rayCount,
        unsigned const batch,
        unsigned const batchCount = defaultForwardRseBatchCount)
    {
        unsigned const end = globalRayOffset + rayCount;
        unsigned const first = globalRayOffset + (batch + batchCount - globalRayOffset % batchCount) % batchCount;
        return first < end ? 1u + (end - 1u - first) / batchCount : 0u;
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

    /**
     * @param applicationSeed Seed shared by the complete launch.
     * @param batch Statistical batch index.
     * @return Deterministically separated seed for the batch.
     */
    ALPAKA_FN_HOST_ACC constexpr unsigned rseBatchSeed(unsigned const applicationSeed, unsigned const batch)
    {
        return static_cast<unsigned>(mixRseBatchSeed64(
            static_cast<std::uint64_t>(applicationSeed) + 0x9e37'79b9ull + 0x85eb'ca6bull * (batch + 1u)));
    }

    /**
     * @param applicationSeed Seed shared by the complete launch.
     * @param batch Statistical batch index.
     * @return Deterministic systematic source offset in `[0, 1)`.
     */
    ALPAKA_FN_HOST_ACC constexpr double rseBatchSourceStratificationOffset(
        unsigned const applicationSeed,
        unsigned const batch)
    {
        return static_cast<double>(rseBatchSeed(rseBatchSeed(applicationSeed, batch), 0x7d3a'9f21u)) / 4294967296.0;
    }

    /**
     * @param applicationSeed Seed shared by the complete launch.
     * @param batch Statistical batch index.
     * @param spectrumSize Number of discrete wavelength bins.
     * @return Batch-specific cyclic spectrum phase.
     */
    ALPAKA_FN_HOST_ACC constexpr unsigned rseBatchSpectrumStratificationPhase(
        unsigned const applicationSeed,
        unsigned const batch,
        unsigned const spectrumSize)
    {
        return spectrumSize == 0u ? 0u
                                  : rseBatchSeed(rseBatchSeed(applicationSeed, batch), 0x6ca4'c37du) % spectrumSize;
    }

    /**
     * @param globalRayIndex Interleaved history index.
     * @param batchCount Total number of statistical batches.
     * @return Index of the history within its batch.
     */
    ALPAKA_FN_HOST_ACC constexpr unsigned rseBatchRayIndex(
        unsigned const globalRayIndex,
        unsigned const batchCount = defaultForwardRseBatchCount)
    {
        return globalRayIndex / batchCount;
    }

    /** @brief Mutable physical and statistical state for one forward ASE history. */
    struct ForwardAseRayState
        : ray::TraversalState
        , ray::SrmPositionStorage<typename std::remove_cvref_t<ALPAKA_TYPEOF(ray::aseSrmPolicy)>::PositionPolicy>
    {
        double weight = 0.0;
        double wavelength = 0.0;
        double accumulatedGain = 1.0;
        unsigned rseBatch = 0u;
    };

    static_assert(!std::derived_from<ForwardAseRayState, ray::BarycentricSrmPositionStorage>);

    /** @brief Device views receiving material-vertex scores and per-cell diagnostics. */
    template<
        alpaka::concepts::IMdSpan<double> TVertexBatchScoreSum,
        alpaka::concepts::IMdSpan<std::uint32_t> TCellRayVisits,
        alpaka::concepts::IMdSpan<std::uint32_t> TCellDroppedRays>
    struct ForwardAccumulationSpans
    {
        TVertexBatchScoreSum vertexBatchScoreSum;
        TCellRayVisits cellRayVisits;
        TCellDroppedRays cellDroppedRays;
    };

} // namespace hase::kernels::forward

namespace hase::kernels::forward
{

    /** @brief Policy-selected per-cell ray-visit diagnostic. */
    struct RecordCellRayVisit
    {
        ALPAKA_FN_ACC void operator()(
            tracePolicy::diagnostics::None,
            alpaka::onAcc::concepts::Acc auto const&,
            alpaka::concepts::IMdSpan<std::uint32_t> auto,
            unsigned) const
        {
        }

        ALPAKA_FN_ACC void operator()(
            tracePolicy::diagnostics::CellRayVisits,
            alpaka::onAcc::concepts::Acc auto const& acc,
            alpaka::concepts::IMdSpan<std::uint32_t> auto cellRayVisits,
            unsigned const tet) const
        {
            alpaka::onAcc::atomicAdd(acc, &cellRayVisits[tet], 1u);
        }
    };

    /** @brief Cell behavior accumulating track-length scores for one ASE ray. */
    template<
        alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> T_Accumulation,
        concepts::TracePolicyList T_TracePolicies>
    struct ForwardAseCellPolicy : ray::behaviourDimension::Cell
    {
        T_Accumulation accumulation;

        ALPAKA_FN_HOST_ACC constexpr ForwardAseCellPolicy(T_Accumulation value, T_TracePolicies) : accumulation{value}
        {
        }

        ALPAKA_FN_ACC bool operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            hase::data::TraceView const& mesh,
            ray::State auto& rayState,
            unsigned const tet,
            Tet4FaceIntersection const intersection)
        {
            double const segmentGain = localSegmentGain(mesh, tet, intersection.length, rayState.wavelength);
            double contribution = rayState.weight * rayState.accumulatedGain;
            contribution *= localSegmentTrackLengthIntegral(mesh, tet, intersection.length, rayState.wavelength);
            if(alpaka::math::isfinite(contribution))
            {
                RecordCellRayVisit{}(T_TracePolicies::getDiagnostics(), acc, accumulation.cellRayVisits, tet);
                auto const weights = segmentMidpointBarycentricVertexWeights(
                    mesh,
                    tet,
                    rayState.position,
                    rayState.direction,
                    intersection.length);
                unsigned const materialVertexOffset = mesh.getMaterialId(tet) * mesh.numberOfMeshPoints;
                for(unsigned localVertex = 0u; localVertex < hase::data::tet4VertexCount; ++localVertex)
                {
                    unsigned const materialVertex
                        = materialVertexOffset + mesh.cellPointIndices[tet * mesh.numberOfCellVertices + localVertex];
                    unsigned const vertex
                        = rayState.rseBatch * (mesh.numberOfMaterials * mesh.numberOfMeshPoints) + materialVertex;
                    double const weight = weights[localVertex];
                    alpaka::onAcc::atomicAdd(acc, &accumulation.vertexBatchScoreSum[vertex], contribution * weight);
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

    /** @brief Failure handler incrementing the dropped-ray count for the last cell. */
    template<alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> T_Accumulation>
    struct CountDroppedForwardRay
    {
        T_Accumulation accumulation;

        ALPAKA_FN_HOST_ACC constexpr explicit CountDroppedForwardRay(T_Accumulation value) : accumulation{value}
        {
        }

        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            hase::data::TraceView const&,
            ray::State auto& rayState)
        {
            alpaka::onAcc::atomicAdd(acc, &accumulation.cellDroppedRays[rayState.cell], 1u);
        }
    };

    /** @brief Kernel launching direct, source-stratified forward ASE histories. */
    struct AccumulateForwardPhiAse
    {
        ALPAKA_FN_HOST_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto const tracePolicies,
            hase::data::TraceView const mesh,
            unsigned const forwardRayCount,
            unsigned const batch,
            double const sourceStrengthTotal,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation,
            unsigned const rngSeed) const
        {
            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                unsigned const batchRayIndex = rayNumber;
                unsigned const batchRayCount = forwardRayCount;
                unsigned const batchSeed = rseBatchSeed(rngSeed, batch);
                auto rndEngine = hase::random::makeRandomEngine(batchSeed, rayHistoryId(0u, batchRayIndex));
                unsigned const tet = sampleStratifiedVolumeBySourceStrength(
                    mesh,
                    sourceStrengthTotal,
                    batchRayIndex,
                    batchRayCount,
                    rseBatchSourceStratificationOffset(rngSeed, batch),
                    rndEngine);
                // Beta-volume sampling carries its importance factor in the
                // source probability, so the compensating weight is one.
                double const sourceWeight = sourceStrengthTotal > 0.0 ? 1.0 : 0.0;
                hase::core::Point origin = samplePointInVolume(mesh, tet, rndEngine);
                hase::core::Point const direction = sampleIsotropicDirection(rndEngine);
                unsigned const material = mesh.getMaterialId(tet);
                unsigned const spectrumSize = mesh.crossSectionCount(material);
                unsigned const spectrumIndex = stratifiedSpectrumIndex(
                    spectrumSize,
                    batchRayIndex,
                    batchRayCount,
                    rseBatchSpectrumStratificationPhase(rngSeed, batch, spectrumSize));
                walkVolumeSeededForwardRay(
                    acc,
                    tracePolicies,
                    mesh,
                    tet,
                    origin,
                    direction,
                    sourceWeight,
                    spectrumSize == 0u ? 0.0 : mesh.emissionWavelength(material, spectrumIndex),
                    batch,
                    accumulation);
            }
        }

        ALPAKA_FN_HOST_ACC void walkVolumeSeededForwardRay(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto const tracePolicies,
            hase::data::TraceView const& mesh,
            unsigned tet,
            hase::core::Point origin,
            hase::core::Point const direction,
            double const sourceWeight,
            double const wavelength,
            unsigned const rseBatch,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation) const
        {
            namespace ray = hase::kernels::forward::ray;
            ForwardAseRayState rayState;
            rayState.position = origin;
            rayState.direction = direction;
            rayState.cell = tet;
            rayState.forbiddenFace = -1;
            rayState.weight = sourceWeight;
            rayState.wavelength = wavelength;
            rayState.rseBatch = rseBatch;
            auto const walkResult = ray::walk(
                acc,
                mesh,
                rayState,
                ray::RayWalkBehaviour{ForwardAseCellPolicy{accumulation, tracePolicies}, ray::BoundaryPolicyEscape{}});
            if(walkResult == ray::WalkResult::failed)
                CountDroppedForwardRay<ALPAKA_TYPEOF(accumulation)>{accumulation}(acc, mesh, rayState);
        }
    };

    /** @brief Kernel launching direct histories while capturing reflected boundary candidates. */
    struct AccumulateForwardPhiAseReflections
    {
        template<alpaka::concepts::SpecializationOf<ReflectionCandidateSpans> T_Candidates>
        struct StoreReflectionBoundary : ray::BoundaryPolicySrm<ray::srmPosition::Centroid>
        {
            T_Candidates candidates;
            std::uint32_t candidateIndex;

            ALPAKA_FN_HOST_ACC constexpr StoreReflectionBoundary(
                T_Candidates candidatesValue,
                std::uint32_t const candidateIndexValue)
                : candidates{candidatesValue}
                , candidateIndex{candidateIndexValue}
            {
            }

            ALPAKA_FN_ACC ray::BoundaryResult operator()(
                alpaka::onAcc::concepts::Acc auto const&,
                hase::data::TraceView const& mesh,
                ray::State auto& rayState,
                unsigned const tet,
                unsigned const localFace)
            {
                core::Direction const normal = outwardFaceNormal(mesh, tet, localFace);
                double const reflectance = boundaryReflectance(mesh, tet, localFace, rayState.direction, normal);
                return storeReflectionCandidate(
                    mesh,
                    rayState,
                    tet,
                    localFace,
                    ReflectionBoundarySample{
                        reflectedDirection(rayState.direction, normal),
                        rayState.weight * rayState.accumulatedGain * reflectance,
                        rayState.wavelength},
                    candidates,
                    candidateIndex);
            }
        };

        ALPAKA_FN_HOST_ACC void walkForwardRay(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto const tracePolicies,
            data::TraceView const& mesh,
            unsigned tet,
            core::Point origin,
            core::Point const direction,
            int const initialForbiddenFace,
            double const sourceWeight,
            double const wavelength,
            unsigned const rseBatch,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation,
            alpaka::concepts::SpecializationOf<ReflectionCandidateSpans> auto candidates,
            std::uint32_t const candidateIndex) const
        {
            ForwardAseRayState rayState;
            rayState.position = origin;
            rayState.direction = direction;
            rayState.cell = tet;
            rayState.forbiddenFace = initialForbiddenFace;
            rayState.weight = sourceWeight;
            rayState.wavelength = wavelength;
            rayState.rseBatch = rseBatch;
            auto const walkResult = ray::walk(
                acc,
                mesh,
                rayState,
                ray::RayWalkBehaviour{
                    ForwardAseCellPolicy{accumulation, tracePolicies},
                    StoreReflectionBoundary<ALPAKA_TYPEOF(candidates)>{candidates, candidateIndex}});
            if(walkResult == ray::WalkResult::failed)
                CountDroppedForwardRay<ALPAKA_TYPEOF(accumulation)>{accumulation}(acc, mesh, rayState);
        }

        ALPAKA_FN_HOST_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto const tracePolicies,
            hase::data::TraceView const mesh,
            unsigned const forwardRayCount,
            unsigned const batch,
            double const sourceStrengthTotal,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation,
            alpaka::concepts::SpecializationOf<ReflectionCandidateSpans> auto candidates,
            unsigned const rngSeed) const
        {
            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                unsigned const batchRayIndex = rayNumber;
                unsigned const batchRayCount = forwardRayCount;
                unsigned const batchSeed = rseBatchSeed(rngSeed, batch);
                auto rndEngine = hase::random::makeRandomEngine(batchSeed, rayHistoryId(0u, batchRayIndex));
                unsigned const tet = sampleStratifiedVolumeBySourceStrength(
                    mesh,
                    sourceStrengthTotal,
                    batchRayIndex,
                    batchRayCount,
                    rseBatchSourceStratificationOffset(rngSeed, batch),
                    rndEngine);
                double const sourceWeight = sourceStrengthTotal > 0.0 ? 1.0 : 0.0;
                core::Point origin = samplePointInVolume(mesh, tet, rndEngine);
                core::Point const direction = sampleIsotropicDirection(rndEngine);
                unsigned const material = mesh.getMaterialId(tet);
                unsigned const spectrumSize = mesh.crossSectionCount(material);
                unsigned const spectrumIndex = stratifiedSpectrumIndex(
                    spectrumSize,
                    batchRayIndex,
                    batchRayCount,
                    rseBatchSpectrumStratificationPhase(rngSeed, batch, spectrumSize));
                walkForwardRay(
                    acc,
                    tracePolicies,
                    mesh,
                    tet,
                    origin,
                    direction,
                    -1,
                    sourceWeight,
                    spectrumSize == 0u ? 0.0 : mesh.emissionWavelength(material, spectrumIndex),
                    batch,
                    accumulation,
                    candidates,
                    batchRayIndex);
            }
        }
    };

    /** @brief Kernel relaunching histories from exact weighted boundary candidates. */
    struct AccumulateReflectedForwardPhiAse
    {
        ALPAKA_FN_HOST_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto const tracePolicies,
            hase::data::TraceView const mesh,
            unsigned const forwardRayCount,
            unsigned const batch,
            double const sourceWeight,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation,
            alpaka::concepts::SpecializationOf<ReflectionCandidateSpans> auto inputCandidates,
            alpaka::concepts::SpecializationOf<ReflectionSamplingSpans> auto sampling,
            alpaka::concepts::SpecializationOf<ReflectionCandidateSpans> auto outputCandidates,
            unsigned const rngSeed,
            unsigned const reflectionPass) const
        {
            AccumulateForwardPhiAseReflections walker;
            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                unsigned const batchRayIndex = rayNumber;
                unsigned const batchSeed = rseBatchSeed(rngSeed, batch);
                ReflectionCandidateSample const sample = sampleReflectionCandidate(
                    inputCandidates,
                    sampling,
                    forwardRayCount,
                    batchRayIndex,
                    forwardRayCount,
                    reflectionResamplingOffset(batchSeed, reflectionPass));
                if(!sample.valid)
                    continue;
                std::uint32_t const candidateIndex = sample.candidateIndex;
                unsigned const faceId = inputCandidates.faceIds[candidateIndex];
                unsigned const tet = faceId / mesh.numberOfFacesPerCell;
                unsigned const localFace = faceId - tet * mesh.numberOfFacesPerCell;
                core::Point const direction = normalize(inputCandidates.directions.at(candidateIndex));
                core::Point const origin = inputCandidates.positions.at(candidateIndex);
                double const wavelength = inputCandidates.wavelengths[candidateIndex];
                walker.walkForwardRay(
                    acc,
                    tracePolicies,
                    mesh,
                    tet,
                    origin,
                    direction,
                    static_cast<int>(localFace),
                    sourceWeight,
                    wavelength,
                    batch,
                    accumulation,
                    outputCandidates,
                    batchRayIndex);
            }
        }
    };

    /** @brief Forward ASE walker that deposits reflected states into a bounded per-face reservoir. */
    struct AccumulateForwardPhiAseSurfaceReservoir
    {
        template<
            alpaka::concepts::SpecializationOf<SurfaceReservoirSpans> T_Reservoir,
            alpaka::rand::concepts::UniformRandomEngine T_Rng>
        struct StoreReflectionBoundary : ray::BoundaryPolicySrm<ray::srmPosition::Centroid>
        {
            T_Reservoir reservoir;
            T_Rng rng;
            std::uint32_t candidateIndex;

            ALPAKA_FN_HOST_ACC constexpr StoreReflectionBoundary(
                T_Reservoir reservoirValue,
                T_Rng rngValue,
                std::uint32_t const candidateIndexValue)
                : reservoir{reservoirValue}
                , rng{rngValue}
                , candidateIndex{candidateIndexValue}
            {
            }

            ALPAKA_FN_ACC ray::BoundaryResult operator()(
                alpaka::onAcc::concepts::Acc auto const& acc,
                hase::data::TraceView const& mesh,
                ray::State auto& rayState,
                unsigned const tet,
                unsigned const localFace)
            {
                core::Direction const normal = outwardFaceNormal(mesh, tet, localFace);
                double const reflectance = boundaryReflectance(mesh, tet, localFace, rayState.direction, normal);
                return storeSurfaceReservoirBoundarySample(
                    acc,
                    mesh,
                    rayState,
                    tet,
                    localFace,
                    SurfaceReservoirBoundarySample{
                        reflectedDirection(rayState.direction, normal),
                        rayState.weight * rayState.accumulatedGain * reflectance,
                        rayState.wavelength},
                    reservoir,
                    candidateIndex,
                    rng);
            }
        };

        ALPAKA_FN_HOST_ACC void walkForwardRay(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto const tracePolicies,
            data::TraceView const& mesh,
            unsigned const tet,
            core::Point const origin,
            core::Point const direction,
            std::int32_t const initialForbiddenFace,
            double const sourceWeight,
            double const wavelength,
            unsigned const rseBatch,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation,
            alpaka::concepts::SpecializationOf<SurfaceReservoirSpans> auto reservoir,
            std::uint32_t const candidateIndex,
            alpaka::rand::concepts::UniformRandomEngine auto rng) const
        {
            ForwardAseRayState rayState;
            rayState.position = origin;
            rayState.direction = direction;
            rayState.cell = tet;
            rayState.forbiddenFace = initialForbiddenFace;
            rayState.weight = sourceWeight;
            rayState.wavelength = wavelength;
            rayState.rseBatch = rseBatch;
            auto const walkResult = ray::walk(
                acc,
                mesh,
                rayState,
                ray::RayWalkBehaviour{
                    ForwardAseCellPolicy{accumulation, tracePolicies},
                    StoreReflectionBoundary<ALPAKA_TYPEOF(reservoir), ALPAKA_TYPEOF(rng)>{
                        reservoir,
                        rng,
                        candidateIndex}});
            if(walkResult == ray::WalkResult::failed)
                CountDroppedForwardRay<ALPAKA_TYPEOF(accumulation)>{accumulation}(acc, mesh, rayState);
        }

        ALPAKA_FN_HOST_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto const tracePolicies,
            hase::data::TraceView const mesh,
            std::uint32_t const forwardRayCount,
            std::uint32_t const batch,
            double const sourceStrengthTotal,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation,
            alpaka::concepts::SpecializationOf<SurfaceReservoirSpans> auto reservoir,
            std::uint32_t const rngSeed) const
        {
            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                std::uint32_t const batchRayIndex = rayNumber;
                std::uint32_t const batchSeed = rseBatchSeed(rngSeed, batch);
                auto rng = hase::random::makeRandomEngine(batchSeed, rayHistoryId(0u, batchRayIndex));
                std::uint32_t const tet = sampleStratifiedVolumeBySourceStrength(
                    mesh,
                    sourceStrengthTotal,
                    batchRayIndex,
                    forwardRayCount,
                    rseBatchSourceStratificationOffset(rngSeed, batch),
                    rng);
                double const sourceWeight = sourceStrengthTotal > 0.0 ? 1.0 : 0.0;
                core::Point const origin = samplePointInVolume(mesh, tet, rng);
                core::Point const direction = sampleIsotropicDirection(rng);
                std::uint32_t const material = mesh.getMaterialId(tet);
                std::uint32_t const spectrumSize = mesh.crossSectionCount(material);
                std::uint32_t const spectrumIndex = stratifiedSpectrumIndex(
                    spectrumSize,
                    batchRayIndex,
                    forwardRayCount,
                    rseBatchSpectrumStratificationPhase(rngSeed, batch, spectrumSize));
                walkForwardRay(
                    acc,
                    tracePolicies,
                    mesh,
                    tet,
                    origin,
                    direction,
                    -1,
                    sourceWeight,
                    spectrumSize == 0u ? 0.0 : mesh.emissionWavelength(material, spectrumIndex),
                    batch,
                    accumulation,
                    reservoir,
                    batchRayIndex,
                    rng);
            }
        }
    };

    /** @brief Relaunch bounded SRM states from either face centroids or retained exact hits. */
    struct AccumulateSurfaceReservoirForwardPhiAse
    {
        ALPAKA_FN_HOST_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto const tracePolicies,
            hase::data::TraceView const mesh,
            std::uint32_t const forwardRayCount,
            std::uint32_t const batch,
            double const sourceWeight,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation,
            alpaka::concepts::SpecializationOf<SurfaceReservoirSpans> auto input,
            alpaka::concepts::SpecializationOf<SurfaceReservoirSamplingCdfSpans> auto sampling,
            alpaka::concepts::SpecializationOf<SurfaceReservoirSpans> auto output,
            std::uint32_t const rngSeed,
            std::uint32_t const reflectionPass) const
        {
            AccumulateForwardPhiAseSurfaceReservoir walker;
            for(auto [rayNumber] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{forwardRayCount}))
            {
                std::uint32_t const batchRayIndex = rayNumber;
                std::uint32_t const batchSeed = rseBatchSeed(rngSeed, batch);
                auto rng = hase::random::makeRandomEngine(batchSeed, rayHistoryId(reflectionPass, batchRayIndex));
                SurfaceReservoirSample const sample = sampleSurfaceReservoir(
                    input,
                    sampling,
                    mesh.numberOfCells * mesh.numberOfFacesPerCell,
                    batchRayIndex,
                    rng);
                if(!sample.valid)
                    continue;

                std::uint32_t const tet = sample.faceId / mesh.numberOfFacesPerCell;
                std::uint32_t const localFace = sample.faceId % mesh.numberOfFacesPerCell;
                core::Point const origin
                    = restoreSurfaceReservoirPosition(input.positionSpans, mesh, tet, localFace, sample.slotIndex);
                walker.walkForwardRay(
                    acc,
                    tracePolicies,
                    mesh,
                    tet,
                    origin,
                    normalize(input.directions.at(sample.slotIndex)),
                    static_cast<std::int32_t>(localFace),
                    sourceWeight,
                    input.wavelengths[sample.slotIndex],
                    batch,
                    accumulation,
                    output,
                    batchRayIndex,
                    rng);
            }
        }
    };
} // namespace hase::kernels::forward
