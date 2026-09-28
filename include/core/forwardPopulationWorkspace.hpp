#pragma once

#include <core/forwardDirect.hpp>
#include <kernels/forward/populationTracing.hpp>

#include <array>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace hase::core
{
    /** Persistent transport scratch. Statistical scores remain in the device context. */
    template<alpaka::onHost::concepts::Device T_Device>
    class ForwardPopulationWorkspace
    {
        using T_Rays
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<ForwardPopulationRay>(std::declval<T_Device&>(), std::size_t{1u}));

        struct Range
        {
            ForwardPopulationBatch work;
            std::uint32_t offset;
        };

    public:
        ForwardPopulationWorkspace(T_Device device) : m_device(std::move(device))
        {
        }

        void prepare(
            concepts::Queue auto const& queue,
            data::TraceView const mesh,
            data::AseDomainSourceView const sources,
            std::span<ForwardPopulationBatch const> const work,
            std::uint32_t const seed)
        {
            m_ranges.clear();
            std::uint32_t count = 0u;
            for(auto const& item : work)
            {
                m_ranges.push_back({item, count});
                count += item.rayCount;
            }
            if(!m_primaries || m_primaryCapacity < count)
            {
                m_primaries = std::make_unique<T_Rays>(
                    alpaka::onHost::alloc<ForwardPopulationRay>(m_device, std::size_t{count}));
                m_primaryCapacity = count;
            }
            for(auto const& range : m_ranges)
                queue.enqueue(
                    getRayFrameSpec(range.work.rayCount, queue),
                    alpaka::KernelBundle{
                        kernels::forward::PrepareForwardPopulation{},
                        mesh,
                        sources,
                        range.work,
                        seed,
                        m_primaries->getView().getSubView(
                            alpaka::Vec{std::size_t{range.offset}},
                            alpaka::Vec{std::size_t{range.work.rayCount}})});
            alpaka::onHost::wait(queue);
        }

        /** Scheduling may change ownership, but never combine two logical reservoirs. */
        void forEachPreparedBatch(auto const& trace)
        {
            for(auto const& range : m_ranges)
                trace(
                    range.work,
                    m_primaries->getView().getSubView(
                        alpaka::Vec{std::size_t{range.offset}},
                        alpaka::Vec{std::size_t{range.work.rayCount}}));
        }

        std::vector<ForwardPopulationRay> trace(
            concepts::Queue auto const& queue,
            data::TraceView const mesh,
            data::AseDomainInterfaceView const interfaces,
            AseTraceControls const& controls,
            alpaka::concepts::SpecializationOf<kernels::forward::ForwardAccumulationSpans> auto accumulation,
            std::uint32_t const population,
            std::span<ForwardPopulationRay const> const incoming,
            bool const primary)
        {
            std::vector<ForwardPopulationRay> result;
            auto traceRange = [&](alpaka::concepts::IView<ForwardPopulationRay> auto input, std::uint32_t const count)
            {
                if(count == 0u)
                    return;
                reserve(count, forwardPopulationDomainCount(controls.domainCount));
                auto const candidateCount = boundaryCandidateCount(count);
                alpaka::onHost::fill(queue, m_scratch->first.weights, 0.0, alpaka::Vec{std::size_t{candidateCount}});
                auto launch = [&](auto diagnostics)
                {
                    queue.enqueue(
                        getRayFrameSpec(count, queue),
                        alpaka::KernelBundle{
                            kernels::forward::TraceForwardPopulation{},
                            kernels::forward::TracePolicyList{
                                kernels::forward::tracePolicy::source::volume,
                                kernels::forward::tracePolicy::cell::forwardAse,
                                kernels::forward::tracePolicy::boundary::boundaryCandidates,
                                kernels::forward::tracePolicy::position::exact,
                                diagnostics},
                            mesh,
                            input,
                            accumulation,
                            m_scratch->first.view(),
                            interfaces,
                            controls.useReflections,
                            std::uint32_t{0u}});
                };
                if(controls.enableDiagnostics)
                    launch(kernels::forward::tracePolicy::diagnostics::enabled);
                else
                    launch(kernels::forward::tracePolicy::diagnostics::none);
                if(!controls.useReflections && controls.domainCount <= 1u)
                    return;
                auto output = m_records->getView().getSubView(alpaka::Vec{std::size_t{candidateCount}});
                queue.enqueue(
                    getRayFrameSpec(candidateCount, queue),
                    alpaka::KernelBundle{
                        kernels::forward::PackForwardPopulationCandidates{},
                        input,
                        m_scratch->first.view(),
                        output});
                auto const previous = result.size();
                result.resize(previous + candidateCount);
                auto host = alpaka::makeView(
                    alpaka::api::host,
                    result.data() + previous,
                    alpaka::Vec{std::size_t{candidateCount}});
                alpaka::onHost::memcpy(queue, host, output);
                // The result vector may grow on the next chunk; complete its outstanding copy first.
                alpaka::onHost::wait(queue);
            };
            if(primary)
            {
                for(auto const& range : m_ranges)
                    if(range.work.rayPopulationId == population)
                        traceRange(
                            m_primaries->getView().getSubView(
                                alpaka::Vec{std::size_t{range.offset}},
                                alpaka::Vec{std::size_t{range.work.rayCount}}),
                            range.work.rayCount);
            }
            else if(!incoming.empty())
            {
                auto const count = static_cast<std::uint32_t>(incoming.size());
                reserve(count, forwardPopulationDomainCount(controls.domainCount));
                auto input = m_incoming->getView().getSubView(alpaka::Vec{std::size_t{count}});
                auto host = alpaka::makeView(alpaka::api::host, incoming.data(), alpaka::Vec{std::size_t{count}});
                alpaka::onHost::memcpy(queue, input, host);
                traceRange(input, count);
            }
            alpaka::onHost::wait(queue);
            return result;
        }

        /** Trace one complete population and return only its scalar surviving weight. */
        double traceResident(
            concepts::Queue auto const& queue,
            data::TraceView const mesh,
            data::AseDomainInterfaceView const interfaces,
            AseTraceControls const& controls,
            alpaka::concepts::SpecializationOf<kernels::forward::ForwardAccumulationSpans> auto accumulation,
            std::uint32_t const population,
            bool const primary)
        {
            if(primary)
            {
                m_activeCount = 0u;
                for(auto const& range : m_ranges)
                    if(range.work.rayPopulationId == population)
                        m_activeCount += range.work.rayCount;
            }
            if(m_activeCount == 0u)
                return 0.0;
            reserve(m_activeCount, forwardPopulationDomainCount(controls.domainCount));
            auto const candidateCount = boundaryCandidateCount(m_activeCount);
            alpaka::onHost::fill(queue, m_scratch->first.weights, 0.0, alpaka::Vec{std::size_t{candidateCount}});
            auto launch = [&](alpaka::concepts::IView<ForwardPopulationRay> auto input,
                              std::uint32_t const offset,
                              auto diagnostics)
            {
                queue.enqueue(
                    getRayFrameSpec(static_cast<std::uint32_t>(input.getExtents().x()), queue),
                    alpaka::KernelBundle{
                        kernels::forward::TraceForwardPopulation{},
                        kernels::forward::TracePolicyList{
                            kernels::forward::tracePolicy::source::volume,
                            kernels::forward::tracePolicy::cell::forwardAse,
                            kernels::forward::tracePolicy::boundary::boundaryCandidates,
                            kernels::forward::tracePolicy::position::exact,
                            diagnostics},
                        mesh,
                        input,
                        accumulation,
                        m_scratch->first.view(),
                        interfaces,
                        controls.useReflections,
                        offset});
            };
            auto trace = [&](auto diagnostics)
            {
                if(primary)
                {
                    for(auto const& range : m_ranges)
                        if(range.work.rayPopulationId == population)
                            launch(
                                m_primaries->getView().getSubView(
                                    alpaka::Vec{std::size_t{range.offset}},
                                    alpaka::Vec{std::size_t{range.work.rayCount}}),
                                range.work.populationOffset + range.work.rayOffset,
                                diagnostics);
                }
                else
                    launch(m_incoming->getView().getSubView(alpaka::Vec{std::size_t{m_activeCount}}), 0u, diagnostics);
            };
            if(controls.enableDiagnostics)
                trace(kernels::forward::tracePolicy::diagnostics::enabled);
            else
                trace(kernels::forward::tracePolicy::diagnostics::none);
            auto const candidateWeights
                = m_scratch->first.weights.getView().getSubView(alpaka::Vec{std::size_t{candidateCount}});
            auto& totalWeight = m_scratch->comb.totalWeight;
            detail::reduce(queue, 0.0, totalWeight, candidateWeights);
            alpaka::onHost::memcpy(queue, m_scratch->totalWeight.toDeviceView(), m_scratch->comb.totalWeight);
            m_scratch->totalWeight.toHost(queue);
            return m_scratch->totalWeightHost[0u];
        }

        /** Canonical candidate ordering precedes device-side population-wide combing. */
        std::vector<ForwardPopulationRay> select(
            alpaka::concepts::SpecializationOf<alpakaUtils::DevBundle> auto& bundle,
            concepts::Queue auto const& queue,
            std::vector<ForwardPopulationRay> candidates,
            std::uint32_t const domains,
            std::uint32_t const seed,
            std::uint32_t const pass)
        {
            std::ranges::sort(candidates, {}, &ForwardPopulationRay::ordinal);
            auto const size = static_cast<std::uint32_t>(candidates.size());
            if(size == 0u)
                return {};
            if(size % 2u != 0u)
                throw std::runtime_error("incomplete forward boundary candidate pairs");
            for(std::uint32_t i = 0u; i < size; ++i)
                if(candidates[i].ordinal != i)
                    throw std::runtime_error("missing or duplicated forward boundary candidate");
            reserve(size / 2u, domains);
            auto records = m_records->getView().getSubView(alpaka::Vec{std::size_t{size}});
            auto host = alpaka::makeView(alpaka::api::host, candidates.data(), alpaka::Vec{std::size_t{size}});
            alpaka::onHost::memcpy(queue, records, host);
            queue.enqueue(
                getRayFrameSpec(size, queue),
                alpaka::KernelBundle{
                    kernels::forward::UnpackForwardPopulationCandidates{},
                    records,
                    m_scratch->first.view()});
            auto const count = enqueueSelection(bundle, queue, size / 2u, domains, seed, pass);
            if(count == 0u)
                return {};
            auto output = m_incoming->getView().getSubView(alpaka::Vec{std::size_t{count}});
            queue.enqueue(
                getRayFrameSpec(count, queue),
                alpaka::KernelBundle{
                    kernels::forward::SelectForwardPopulationRays{},
                    records,
                    m_scratch->comb.selectedView(count),
                    m_scratch->comb.selectedWeightsView(count),
                    output,
                    pass});
            std::vector<ForwardPopulationRay> result(count);
            auto resultHost = alpaka::makeView(alpaka::api::host, result.data(), alpaka::Vec{std::size_t{count}});
            alpaka::onHost::memcpy(queue, resultHost, output);
            alpaka::onHost::wait(queue);
            return result;
        }

        /** Comb resident candidates into the next resident input without moving ray records to the host. */
        void selectResident(
            alpaka::concepts::SpecializationOf<alpakaUtils::DevBundle> auto& bundle,
            concepts::Queue auto const& queue,
            std::uint32_t const domains,
            std::uint32_t const seed,
            std::uint32_t const pass)
        {
            auto const count = enqueueSelection(bundle, queue, m_activeCount, domains, seed, pass);
            m_activeCount = count;
            if(count == 0u)
                return;
            queue.enqueue(
                getRayFrameSpec(count, queue),
                alpaka::KernelBundle{
                    kernels::forward::SelectResidentForwardPopulationRays{},
                    m_scratch->first.view(),
                    m_scratch->comb.selectedView(count),
                    m_scratch->comb.selectedWeightsView(count),
                    m_incoming->getView().getSubView(alpaka::Vec{std::size_t{count}}),
                    pass});
        }

    private:
        std::uint32_t enqueueSelection(
            alpaka::concepts::SpecializationOf<alpakaUtils::DevBundle> auto& bundle,
            concepts::Queue auto const& queue,
            std::uint32_t const parentCount,
            std::uint32_t const domains,
            std::uint32_t const seed,
            std::uint32_t const pass)
        {
            auto& scratch = *m_scratch;
            auto const size = boundaryCandidateCount(parentCount);
            scratch.comb.enqueueRouteMeasurements(
                bundle,
                queue,
                scratch.first.weights.getView(),
                scratch.first.targetDomains.getView(),
                parentCount,
                domains);
            alpaka::onHost::memcpy(queue, scratch.liveCount.toDeviceView(), scratch.comb.liveParentCount);
            alpaka::onHost::memcpy(queue, scratch.routeWeights.toDeviceView(), scratch.comb.routeWeights);
            alpaka::onHost::memcpy(
                queue,
                scratch.routeCandidateCounts.toDeviceView(),
                scratch.comb.routeCandidateCounts);
            scratch.liveCount.toHost(queue);
            scratch.routeWeights.toHost(queue);
            scratch.routeCandidateCounts.toHost(queue);
            auto const count = scratch.liveCountHost[0u];
            if(count == 0u)
                return 0u;
            auto const nonEmpty
                = std::ranges::count_if(scratch.routeCandidateCountsHost, [](auto n) { return n > 0u; });
            if(static_cast<std::uint32_t>(nonEmpty) > count)
                scratch.comb.enqueue(bundle, queue, scratch.first.weights.getView(), size, count, seed, pass);
            else
            {
                auto const counts = allocateBoundaryRoutePopulations(
                    std::span<double const>{scratch.routeWeightsHost},
                    count,
                    std::span<std::uint32_t const>{scratch.routeCandidateCountsHost});
                std::uint32_t offset = 0u;
                for(std::uint32_t domain = 0u; domain < domains; ++domain)
                    if(counts[domain] != 0u)
                    {
                        scratch.comb.enqueueSpatialDomain(
                            bundle,
                            queue,
                            scratch.first.view(),
                            size,
                            scratch.routeCandidateCountsHost[domain],
                            domain,
                            offset,
                            counts[domain],
                            seed,
                            (static_cast<std::uint64_t>(pass) << 32u) | domain);
                        offset += counts[domain];
                    }
            }
            return count;
        }

        void reserve(std::uint32_t const count, std::uint32_t const domains)
        {
            if(m_scratch && m_capacity >= count && m_domains == domains)
                return;
            auto const capacity = boundaryCandidateCount(count);
            m_scratch = std::make_unique<DirectBoundaryScratch<T_Device>>(m_device, count, domains);
            m_records = std::make_unique<T_Rays>(
                alpaka::onHost::alloc<ForwardPopulationRay>(m_device, std::size_t{capacity}));
            m_incoming
                = std::make_unique<T_Rays>(alpaka::onHost::alloc<ForwardPopulationRay>(m_device, std::size_t{count}));
            m_capacity = count;
            m_domains = domains;
        }

        T_Device m_device;
        std::unique_ptr<T_Rays> m_primaries, m_records, m_incoming;
        std::unique_ptr<DirectBoundaryScratch<T_Device>> m_scratch;
        std::vector<Range> m_ranges;
        std::uint32_t m_primaryCapacity{}, m_capacity{}, m_domains{};
        std::uint32_t m_activeCount{};
    };
} // namespace hase::core
