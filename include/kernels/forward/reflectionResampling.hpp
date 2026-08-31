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
#include <kernels/forward/policyRay.hpp>

#include <cstdint>

namespace hase::kernels::forward
{
    /** @brief Ray-indexed reflected boundary candidates for one transport pass. */
    template<
        alpaka::concepts::IView<double> T_CartesianView,
        alpaka::concepts::IView<double> T_Weights,
        alpaka::concepts::IView<double> T_Wavelengths,
        alpaka::concepts::IView<std::uint32_t> T_FaceIds>
    struct ReflectionCandidateSpans
    {
        hase::core::PositionViewSoA<T_CartesianView> positions;
        hase::core::DirectionViewSoA<T_CartesianView> directions;
        T_Weights weights;
        T_Wavelengths wavelengths;
        T_FaceIds faceIds;
    };

    /** @brief Ordered cumulative weights and their total for reflected candidates. */
    template<alpaka::concepts::IView<double> T_Cdf, alpaka::concepts::IView<double> T_TotalWeight>
    struct ReflectionSamplingSpans
    {
        T_Cdf cdf;
        T_TotalWeight totalWeight;
    };

    /** @brief Reflected boundary state offered to the next transport pass. */
    struct ReflectionBoundarySample
    {
        hase::core::Direction direction;
        double weight = 0.0;
        double wavelength = 0.0;
        bool valid = true;
    };

    /** @brief Selected reflected boundary candidate. */
    struct ReflectionCandidateSample
    {
        std::uint32_t candidateIndex = 0u;
        bool valid = false;
    };

    /**
     * @param rngSeed Seed assigned to one independent statistical batch.
     * @param pass Reflection pass index.
     * @return Deterministic systematic-resampling offset in `(0, 1)`.
     */
    [[nodiscard]] ALPAKA_FN_HOST_ACC constexpr double reflectionResamplingOffset(
        std::uint32_t const rngSeed,
        std::uint32_t const pass)
    {
        std::uint64_t value = (static_cast<std::uint64_t>(rngSeed) << 32u) | pass;
        value ^= value >> 30u;
        value *= 0xbf58'476d'1ce4'e5b9ull;
        value ^= value >> 27u;
        value *= 0x94d0'49bb'1331'11ebull;
        value ^= value >> 31u;
        std::uint32_t const sample = static_cast<std::uint32_t>(value >> 32u);
        return (static_cast<double>(sample) + 0.5) * (1.0 / 4294967296.0);
    }

    /**
     * @brief Store one ray-owned reflected boundary candidate and terminate the walk.
     * @param mesh Device trace defining the flat face layout.
     * @param rayState Current ray state supplying the boundary position.
     * @param cell Current cell index.
     * @param localFace Hit face within the current cell.
     * @param sample Direction, weight, wavelength, and validity to store.
     * @param candidates Output candidate views.
     * @param candidateIndex Unique ray-owned candidate slot for this pass.
     * @return `BoundaryResult::stop` after handling the boundary.
     */
    ALPAKA_FN_ACC ray::BoundaryResult storeReflectionCandidate(
        hase::data::TraceView const& mesh,
        ray::State auto& rayState,
        std::uint32_t const cell,
        std::uint32_t const localFace,
        ReflectionBoundarySample const sample,
        alpaka::concepts::SpecializationOf<ReflectionCandidateSpans> auto candidates,
        std::uint32_t const candidateIndex)
    {
        if(!sample.valid || sample.weight <= 0.0 || !alpaka::math::isfinite(sample.weight))
            return ray::BoundaryResult::stop;

        candidates.positions.x[candidateIndex] = rayState.position.x;
        candidates.positions.y[candidateIndex] = rayState.position.y;
        candidates.positions.z[candidateIndex] = rayState.position.z;
        candidates.directions.x[candidateIndex] = sample.direction.x;
        candidates.directions.y[candidateIndex] = sample.direction.y;
        candidates.directions.z[candidateIndex] = sample.direction.z;
        candidates.weights[candidateIndex] = sample.weight;
        candidates.wavelengths[candidateIndex] = sample.wavelength;
        candidates.faceIds[candidateIndex] = cell * mesh.numberOfFacesPerCell + localFace;
        return ray::BoundaryResult::stop;
    }

    /** @brief Map invalid reflected-candidate weights to zero before building the sampling CDF. */
    struct FilterReflectionSamplingWeight
    {
        [[nodiscard]] ALPAKA_FN_HOST_ACC constexpr double operator()(double const weight) const
        {
            return weight > 0.0 && alpaka::math::isfinite(weight) ? weight : 0.0;
        }
    };

    /** @brief Capture the final unnormalized reflected-candidate prefix. */
    struct CaptureReflectionSamplingTotalWeight
    {
        ALPAKA_FN_HOST_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const&,
            std::uint32_t const candidateCount,
            alpaka::concepts::SpecializationOf<ReflectionSamplingSpans> auto sampling) const
        {
            sampling.totalWeight[0u] = candidateCount == 0u ? 0.0 : sampling.cdf[candidateCount - 1u];
        }
    };

    /**
     * @brief Locate one deterministic systematic-resampling target in an ordered CDF.
     * @param cdf Unnormalized non-decreasing cumulative candidate weights.
     * @param candidateCount Number of candidate entries in `cdf`.
     * @param totalWeight Final cumulative weight.
     * @param rayIndex Relaunched ray index within the pass.
     * @param rayCount Total number of relaunched rays.
     * @param offset Seed-derived systematic offset in `(0, 1)`.
     * @return Candidate index whose cumulative interval contains the target.
     */
    template<typename T_Cdf>
    [[nodiscard]] ALPAKA_FN_HOST_ACC std::uint32_t reflectionCandidateIndex(
        T_Cdf const& cdf,
        std::uint32_t const candidateCount,
        double const totalWeight,
        std::uint32_t const rayIndex,
        std::uint32_t const rayCount,
        double const offset)
    {
        double const target = (static_cast<double>(rayIndex) + offset) * (totalWeight / static_cast<double>(rayCount));
        std::uint32_t lower = 0u;
        std::uint32_t upper = candidateCount;
        while(lower < upper)
        {
            std::uint32_t const middle = lower + (upper - lower) / 2u;
            if(cdf[middle] <= target)
                lower = middle + 1u;
            else
                upper = middle;
        }
        return lower < candidateCount ? lower : candidateCount - 1u;
    }

    /**
     * @brief Select one candidate by deterministic systematic weighted resampling.
     * @param candidates Candidate weights defining valid entries.
     * @param sampling Unnormalized ordered cumulative weights and their total.
     * @param candidateCount Number of candidate slots populated for the pass.
     * @param rayIndex Relaunched ray index within the pass.
     * @param rayCount Total number of relaunched rays.
     * @param offset Seed-derived systematic offset in `(0, 1)`.
     * @return Selected candidate, or an invalid sample if the total weight is zero.
     */
    [[nodiscard]] ALPAKA_FN_HOST_ACC ReflectionCandidateSample sampleReflectionCandidate(
        alpaka::concepts::SpecializationOf<ReflectionCandidateSpans> auto candidates,
        alpaka::concepts::SpecializationOf<ReflectionSamplingSpans> auto sampling,
        std::uint32_t const candidateCount,
        std::uint32_t const rayIndex,
        std::uint32_t const rayCount,
        double const offset)
    {
        double const totalWeight = sampling.totalWeight[0u];
        if(candidateCount == 0u || rayCount == 0u || totalWeight <= 0.0)
            return {};

        std::uint32_t const candidateIndex
            = reflectionCandidateIndex(sampling.cdf, candidateCount, totalWeight, rayIndex, rayCount, offset);
        return {candidateIndex, candidates.weights[candidateIndex] > 0.0};
    }
} // namespace hase::kernels::forward
