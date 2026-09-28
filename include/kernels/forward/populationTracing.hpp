#pragma once

#include <core/forwardPopulationPlan.hpp>
#include <kernels/forward/directBoundary.hpp>

#include <cstdint>

namespace hase::core
{
    /** Transport record only: no per-history field or moment storage. */
    struct ForwardPopulationRay
    {
        Position position{};
        Direction direction{};
        Point faceBarycentric{};
        double weight{};
        double wavelength{};
        data::DomainId domainId{};
        std::uint32_t cell{};
        std::int32_t face{-1};
        std::uint32_t rayPopulationId{};
        std::uint32_t depth{};
        std::uint64_t historyId{};
        std::uint32_t ordinal{};
    };
} // namespace hase::core

namespace hase::kernels::forward
{
    struct PrepareForwardPopulation
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            data::TraceView const mesh,
            data::AseDomainSourceView const sources,
            core::ForwardPopulationBatch const work,
            std::uint32_t const seed,
            alpaka::concepts::IView<core::ForwardPopulationRay> auto output) const
        {
            auto const begin = sources.offsets[work.domainId];
            auto const end = sources.offsets[work.domainId + 1u];
            auto const total = sources.sourceStrengthTotals[work.domainId];
            for(auto [index] :
                alpaka::onAcc::makeIdxMap(acc, alpaka::onAcc::worker::threadsInGrid, alpaka::IdxRange{work.rayCount}))
            {
                auto const ray = work.rayOffset + index;
                auto const id = rayHistoryId(work.domainId, ray);
                // A zero-source allocation never dereferences an empty source CDF.
                // The invalid-cell sentinel is skipped by transport because its weight is zero.
                core::ForwardPopulationRay prepared{
                    .weight = work.sourceWeight,
                    .domainId = work.domainId,
                    .cell = mesh.numberOfCells,
                    .rayPopulationId = work.rayPopulationId,
                    .historyId = id,
                    .ordinal = work.populationOffset + ray};
                if(total > 0.0 && begin < end && work.sourceWeight > 0.0)
                {
                    auto rng = alpaka::rand::engine::Philox4x32x10{rayPopulationSeed(seed, work.rayPopulationId), id};
                    auto const target = total
                                        * stratifiedUnitInterval(
                                            ray,
                                            work.domainRayCount,
                                            rayPopulationSourceStratificationOffset(seed, work.rayPopulationId));
                    auto lower = begin;
                    auto upper = end;
                    while(lower < upper)
                    {
                        auto const middle = lower + (upper - lower) / 2u;
                        if(sources.sourceStrengthPrefix[middle] <= target)
                            lower = middle + 1u;
                        else
                            upper = middle;
                    }
                    prepared.cell = sources.globalCells[lower < end ? lower : end - 1u];
                    auto const material = mesh.getMaterialId(prepared.cell);
                    auto const spectrumSize = mesh.crossSectionCount(material);
                    auto const spectrum = stratifiedSpectrumIndex(
                        spectrumSize,
                        ray,
                        work.domainRayCount,
                        rayPopulationSpectrumStratificationPhase(seed, work.rayPopulationId, spectrumSize),
                        rayPopulationSpectrumPermutationSeed(seed, work.rayPopulationId));
                    prepared.position = samplePointInVolume(mesh, prepared.cell, rng);
                    prepared.direction = sampleIsotropicDirection(rng);
                    prepared.wavelength = spectrumSize == 0u ? 0.0 : mesh.emissionWavelength(material, spectrum);
                }
                output[index] = prepared;
            }
        }
    };

    struct TraceForwardPopulation
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto policies,
            data::TraceView const mesh,
            alpaka::concepts::IView<core::ForwardPopulationRay> auto input,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation,
            alpaka::concepts::SpecializationOf<core::BoundaryRaySpans> auto candidates,
            data::AseDomainInterfaceView const interfaces,
            bool const reflections,
            std::uint32_t const candidateOffset) const
        {
            for(auto [index] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{input.getExtents().x()}))
            {
                auto const ray = input[index];
                if(ray.weight <= 0.0)
                    continue;
                AccumulateForwardPhiAseDirect{}.walk(
                    acc,
                    policies,
                    mesh,
                    ray.cell,
                    ray.position,
                    ray.direction,
                    ray.face,
                    ray.weight,
                    ray.wavelength,
                    ray.rayPopulationId,
                    accumulation,
                    candidates,
                    candidateOffset + static_cast<std::uint32_t>(index),
                    interfaces,
                    ray.domainId,
                    ray.depth + 1u,
                    ray.historyId,
                    reflections);
            }
        }
    };

    /** Trace already stratified primary rays into one logical batch's bounded reservoir. */
    struct TracePreparedForwardSrm
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            concepts::TracePolicyList auto policies,
            data::TraceView const mesh,
            alpaka::concepts::IView<core::ForwardPopulationRay> auto input,
            alpaka::concepts::SpecializationOf<ForwardAccumulationSpans> auto accumulation,
            alpaka::concepts::SpecializationOf<SurfaceReservoirSpans> auto reservoir,
            data::AseDomainInterfaceView const interfaces,
            std::uint32_t const seed,
            bool const reflections) const
        {
            for(auto [index] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{input.getExtents().x()}))
            {
                auto const ray = input[index];
                if(ray.weight <= 0.0)
                    continue;
                auto rng = random::makeRandomEngine(seed, ray.historyId);
                AccumulateForwardPhiAseReservoir{}.walkForwardRay(
                    acc,
                    policies,
                    mesh,
                    ray.cell,
                    ray.position,
                    ray.direction,
                    ray.face,
                    ray.weight,
                    ray.wavelength,
                    ray.rayPopulationId,
                    accumulation,
                    reservoir,
                    static_cast<std::uint32_t>(index),
                    rng,
                    interfaces,
                    reflections);
            }
        }
    };

    /** Pack only initialized candidate payloads, preserving the global candidate ordering. */
    struct PackForwardPopulationCandidates
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            alpaka::concepts::IView<core::ForwardPopulationRay> auto input,
            alpaka::concepts::SpecializationOf<core::BoundaryRaySpans> auto candidates,
            alpaka::concepts::IView<core::ForwardPopulationRay> auto output) const
        {
            for(auto [index] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{output.getExtents().x()}))
            {
                core::ForwardPopulationRay ray{
                    .ordinal = 2u * input[index / 2u].ordinal + static_cast<std::uint32_t>(index % 2u)};
                ray.weight = candidates.weights[index];
                if(ray.weight > 0.0)
                {
                    ray.position = candidates.positions.at(index);
                    ray.direction = candidates.directions.at(index);
                    ray.faceBarycentric = candidates.faceBarycentric.at(index);
                    ray.wavelength = candidates.wavelengths[index];
                    ray.domainId = candidates.targetDomains[index];
                    ray.cell = candidates.targetCells[index];
                    ray.face = static_cast<std::int32_t>(candidates.targetFaces[index]);
                    ray.rayPopulationId = candidates.rayPopulationIds[index];
                    ray.depth = candidates.reflectionDepths[index];
                    ray.historyId = candidates.historyIds[index];
                }
                output[index] = ray;
            }
        }
    };

    struct UnpackForwardPopulationCandidates
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            alpaka::concepts::IView<core::ForwardPopulationRay> auto input,
            alpaka::concepts::SpecializationOf<core::BoundaryRaySpans> auto output) const
        {
            for(auto [i] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{input.getExtents().x()}))
            {
                auto const ray = input[i];
                output.positions.x[i] = ray.position.x;
                output.positions.y[i] = ray.position.y;
                output.positions.z[i] = ray.position.z;
                output.directions.x[i] = ray.direction.x;
                output.directions.y[i] = ray.direction.y;
                output.directions.z[i] = ray.direction.z;
                output.faceBarycentric.x[i] = ray.faceBarycentric.x;
                output.faceBarycentric.y[i] = ray.faceBarycentric.y;
                output.faceBarycentric.z[i] = ray.faceBarycentric.z;
                output.weights[i] = ray.weight;
                output.wavelengths[i] = ray.wavelength;
                output.targetDomains[i] = ray.domainId;
                output.targetCells[i] = ray.cell;
                output.targetFaces[i] = static_cast<std::uint32_t>(ray.face);
                output.rayPopulationIds[i] = ray.rayPopulationId;
                output.reflectionDepths[i] = ray.depth;
                output.historyIds[i] = ray.historyId;
            }
        }
    };

    struct SelectForwardPopulationRays
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            alpaka::concepts::IView<core::ForwardPopulationRay> auto candidates,
            alpaka::concepts::IView<std::uint32_t> auto selected,
            alpaka::concepts::IView<double> auto weights,
            alpaka::concepts::IView<core::ForwardPopulationRay> auto output,
            std::uint32_t const pass) const
        {
            for(auto [i] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{output.getExtents().x()}))
            {
                auto ray = candidates[selected[i]];
                ray.weight = weights[i];
                ray.direction = normalize(ray.direction);
                ray.ordinal = static_cast<std::uint32_t>(i);
                ray.historyId ^= rayHistoryId(pass, ray.ordinal);
                output[i] = ray;
            }
        }
    };

    /** Select from resident SoA candidates without packing through a host-visible record. */
    struct SelectResidentForwardPopulationRays
    {
        ALPAKA_FN_ACC void operator()(
            alpaka::onAcc::concepts::Acc auto const& acc,
            alpaka::concepts::SpecializationOf<core::BoundaryRaySpans> auto candidates,
            alpaka::concepts::IView<std::uint32_t> auto selected,
            alpaka::concepts::IView<double> auto weights,
            alpaka::concepts::IView<core::ForwardPopulationRay> auto output,
            std::uint32_t const pass) const
        {
            for(auto [i] : alpaka::onAcc::makeIdxMap(
                    acc,
                    alpaka::onAcc::worker::threadsInGrid,
                    alpaka::IdxRange{output.getExtents().x()}))
            {
                auto const candidate = selected[i];
                core::ForwardPopulationRay ray{
                    .position = candidates.positions.at(candidate),
                    .direction = normalize(candidates.directions.at(candidate)),
                    .faceBarycentric = candidates.faceBarycentric.at(candidate),
                    .weight = weights[i],
                    .wavelength = candidates.wavelengths[candidate],
                    .domainId = candidates.targetDomains[candidate],
                    .cell = candidates.targetCells[candidate],
                    .face = static_cast<std::int32_t>(candidates.targetFaces[candidate]),
                    .rayPopulationId = candidates.rayPopulationIds[candidate],
                    .depth = candidates.reflectionDepths[candidate],
                    .historyId = candidates.historyIds[candidate],
                    .ordinal = static_cast<std::uint32_t>(i)};
                ray.historyId ^= rayHistoryId(pass, ray.ordinal);
                output[i] = ray;
            }
        }
    };
} // namespace hase::kernels::forward
