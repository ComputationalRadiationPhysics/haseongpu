#pragma once

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/HybridBuffer.hpp>
#include <concepts/concepts.hpp>
#include <core/Runtime.hpp>
#include <core/boundaryRayBuffer.hpp>
#include <core/domainSchedule.hpp>
#include <core/forwardSrm.hpp>
#include <core/particleComb.hpp>
#include <kernels/forward/directBoundary.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace hase::core
{
    /** @brief Persistent device buffers for population-controlled direct boundary tracing. */
    template<alpaka::onHost::concepts::Device T_Device>
    struct DirectBoundaryScratch
    {
        using T_TotalWeightBuffer = alpakaUtils::GetHybridBuffer_t<T_Device, std::array<double, 1u>>;
        using T_LiveCountBuffer = alpakaUtils::GetHybridBuffer_t<T_Device, std::array<std::uint32_t, 1u>>;
        using T_RouteWeightBuffer = alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<double>>;
        using T_RouteCountBuffer = alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<std::uint32_t>>;

        DirectBoundaryScratch(T_Device& device, std::uint32_t const capacity, std::uint32_t const domainCount)
            : first(device, boundaryCandidateCount(capacity))
            , second(device, boundaryCandidateCount(capacity))
            , comb(device, boundaryCandidateCount(capacity), capacity, domainCount)
            , totalWeightHost{}
            , totalWeight(alpakaUtils::getHybridBuffer(device, totalWeightHost))
            , liveCountHost{}
            , liveCount(alpakaUtils::getHybridBuffer(device, liveCountHost))
            , routeWeightsHost(domainCount, 0.0)
            , routeWeights(alpakaUtils::getHybridBuffer(device, routeWeightsHost))
            , routeCandidateCountsHost(domainCount, 0u)
            , routeCandidateCounts(alpakaUtils::getHybridBuffer(device, routeCandidateCountsHost))
            , m_capacity(capacity)
        {
        }

        void validate(std::uint32_t const rayCount) const
        {
            if(rayCount > m_capacity)
                throw std::out_of_range("direct boundary scratch capacity exceeded");
        }

        BoundaryRayBuffer<T_Device> first;
        BoundaryRayBuffer<T_Device> second;
        ParticleCombWorkspace<T_Device> comb;
        std::array<double, 1u> totalWeightHost;
        T_TotalWeightBuffer totalWeight;
        std::array<std::uint32_t, 1u> liveCountHost;
        T_LiveCountBuffer liveCount;
        std::vector<double> routeWeightsHost;
        T_RouteWeightBuffer routeWeights;
        std::vector<std::uint32_t> routeCandidateCountsHost;
        T_RouteCountBuffer routeCandidateCounts;

    private:
        std::uint32_t m_capacity;
    };

    /**
     * @brief Trace direct boundary passes while combing every outgoing population on the device.
     *
     * Large candidate, prefix, selection, and relaunch buffers remain device-resident. One scalar
     * total is downloaded at pass boundaries for convergence and pass termination.
     */
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Exec>
    void runForwardDirect(
        alpakaUtils::DevBundle<T_Device, T_Exec>& devBundle,
        concepts::Queue auto const& queue,
        hase::data::TraceView const mesh,
        AseTraceControls const& experiment,
        ForwardPhiAseRawResult& result,
        unsigned const rayCount,
        unsigned const rayPopulationId,
        double const sourceStrengthTotal,
        alpaka::concepts::IBuffer<double> auto& vertexPopulationScoreSum,
        alpaka::concepts::IBuffer<std::uint32_t> auto& volumeRayVisits,
        alpaka::concepts::IBuffer<std::uint32_t> auto& droppedRays,
        unsigned const rngSeed,
        DirectBoundaryScratch<T_Device>& scratch,
        hase::kernels::forward::concepts::TracePolicy auto const diagnostics,
        hase::data::AseDomainInterfaceView const interfaceMap = {},
        hase::data::AseDomainSourceView const domainSources = {},
        std::span<std::uint32_t const> const domainRayCounts = {},
        std::span<double const> const domainSourceWeights = {},
        std::span<std::uint32_t const> const domainPopulationCounts = {})
    {
        scratch.validate(rayCount);
        std::uint32_t populationRayCount = rayCount;
        std::uint32_t const domainCount
            = domainPopulationCounts.empty() ? 0u : static_cast<std::uint32_t>(domainPopulationCounts.size());
        auto accumulation = hase::kernels::forward::ForwardAccumulationSpans{
            vertexPopulationScoreSum.getMdSpan(),
            volumeRayVisits.getMdSpan(),
            droppedRays.getMdSpan()};
        auto const frameSpec = getRayFrameSpec(rayCount, queue);
        auto readTotalWeight = [&]()
        {
            alpaka::onHost::memcpy(queue, scratch.totalWeight.toDeviceView(), scratch.comb.totalWeight);
            scratch.totalWeight.toHost(queue);
            return scratch.totalWeight.getHostView()[0u];
        };

        bool useGlobalComb = domainCount == 0u;
        auto measureRoutes = [&](auto& candidates, std::uint32_t const parentCount)
        {
            if(domainCount == 0u)
                return std::vector<std::uint32_t>{parentCount};
            scratch.comb.enqueueRouteMeasurements(
                devBundle,
                queue,
                candidates.weights.getView(),
                candidates.targetDomains.getView(),
                parentCount,
                domainCount);
            alpaka::onHost::memcpy(queue, scratch.liveCount.toDeviceView(), scratch.comb.liveParentCount);
            alpaka::onHost::memcpy(queue, scratch.routeWeights.toDeviceView(), scratch.comb.routeWeights);
            alpaka::onHost::memcpy(
                queue,
                scratch.routeCandidateCounts.toDeviceView(),
                scratch.comb.routeCandidateCounts);
            scratch.liveCount.toHost(queue);
            scratch.routeWeights.toHost(queue);
            scratch.routeCandidateCounts.toHost(queue);
            auto const nonEmptyRoutes
                = std::ranges::count_if(scratch.routeCandidateCountsHost, [](auto const count) { return count > 0u; });
            useGlobalComb = static_cast<std::size_t>(nonEmptyRoutes) > scratch.liveCountHost[0u];
            if(useGlobalComb)
                return std::vector<std::uint32_t>{scratch.liveCountHost[0u]};
            return allocateBoundaryRoutePopulations(
                std::span<double const>{scratch.routeWeightsHost},
                scratch.liveCountHost[0u],
                std::span<std::uint32_t const>{scratch.routeCandidateCountsHost});
        };

        auto enqueueComb = [&](auto& candidates,
                               std::uint32_t const candidateCount,
                               std::span<std::uint32_t const> const routeCounts,
                               std::uint32_t const pass)
        {
            if(useGlobalComb)
            {
                scratch.comb.enqueue(
                    devBundle,
                    queue,
                    candidates.weights.getView(),
                    candidateCount,
                    std::accumulate(routeCounts.begin(), routeCounts.end(), std::uint32_t{0u}),
                    rngSeed,
                    static_cast<std::uint64_t>(pass));
                return;
            }
            std::uint32_t outputOffset = 0u;
            for(std::uint32_t domain = 0u; domain < routeCounts.size(); ++domain)
            {
                auto const count = routeCounts[domain];
                if(count == 0u)
                    continue;
                scratch.comb.enqueueSpatialDomain(
                    devBundle,
                    queue,
                    candidates.view(),
                    candidateCount,
                    scratch.routeCandidateCountsHost[domain],
                    domain,
                    outputOffset,
                    count,
                    rngSeed,
                    (static_cast<std::uint64_t>(pass) << 32u) | domain);
                outputOffset += count;
            }
            auto const expectedOutputCount
                = std::accumulate(routeCounts.begin(), routeCounts.end(), std::uint32_t{0u});
            if(outputOffset != expectedOutputCount)
                throw std::runtime_error("domain population counts do not match the relaunch population");
            scratch.comb.enqueueTotal(devBundle, queue, candidates.weights.getView(), candidateCount);
        };

        auto const initialCandidateCount = boundaryCandidateCount(rayCount);
        alpaka::onHost::fill(
            queue,
            scratch.first.weights,
            0.0,
            alpaka::Vec{static_cast<std::size_t>(initialCandidateCount)});
        if(domainRayCounts.empty())
            queue.enqueue(
                frameSpec,
                alpaka::KernelBundle{
                    hase::kernels::forward::AccumulateForwardPhiAseDirect{},
                    hase::kernels::forward::TracePolicyList{
                        hase::kernels::forward::tracePolicy::source::volume,
                        hase::kernels::forward::tracePolicy::cell::forwardAse,
                        hase::kernels::forward::tracePolicy::boundary::boundaryCandidates,
                        hase::kernels::forward::tracePolicy::position::exact,
                        diagnostics},
                    mesh,
                    rayCount,
                    rayPopulationId,
                    sourceStrengthTotal,
                    accumulation,
                    scratch.first.view(),
                    rngSeed,
                    interfaceMap,
                    experiment.useReflections});
        else
        {
            unsigned candidateOffset = 0u;
            for(std::uint32_t domain = 0u; domain < domainRayCounts.size(); ++domain)
            {
                auto const count = domainRayCounts[domain];
                if(count == 0u)
                    continue;
                auto const domainFrame = getRayFrameSpec(count, queue);
                double const sourceWeight = domainSourceWeights[domain];
                queue.enqueue(
                    domainFrame,
                    alpaka::KernelBundle{
                        hase::kernels::forward::AccumulateDomainForwardPhiAseDirect{},
                        hase::kernels::forward::TracePolicyList{
                            hase::kernels::forward::tracePolicy::source::volume,
                            hase::kernels::forward::tracePolicy::cell::forwardAse,
                            hase::kernels::forward::tracePolicy::boundary::boundaryCandidates,
                            hase::kernels::forward::tracePolicy::position::exact,
                            diagnostics},
                        mesh,
                        domainSources,
                        domain,
                        count,
                        candidateOffset,
                        rayPopulationId,
                        sourceWeight,
                        accumulation,
                        scratch.first.view(),
                        rngSeed,
                        interfaceMap,
                        experiment.useReflections});
                candidateOffset += count;
            }
            if(candidateOffset != rayCount)
                throw std::runtime_error("domain batch ray counts do not match the launch total");
        }
        auto routeCounts = measureRoutes(scratch.first, rayCount);
        populationRayCount = std::accumulate(routeCounts.begin(), routeCounts.end(), std::uint32_t{0u});
        scratch.validate(populationRayCount);
        enqueueComb(scratch.first, initialCandidateCount, routeCounts, 0u);

        double const initialWeight = readTotalWeight();
        result.boundaryMaxPasses = experiment.resolvedBoundaryMaxPasses();
        if(initialWeight <= 0.0)
        {
            result.boundaryStatus = data::BoundaryStatus::converged;
            return;
        }

        bool inputFirst = true;
        double previousWeight = initialWeight;
        std::vector<double> residualWeightFractions{1.0};
        unsigned growCount = 0u;
        constexpr unsigned divergenceStreak = 3u;
        result.boundaryDivergenceStreak = divergenceStreak;
        result.boundaryStatus = data::BoundaryStatus::maxPasses;
        result.boundaryRemainingFraction = 1.0;
        for(unsigned pass = 1u; pass <= result.boundaryMaxPasses; ++pass)
        {
            if(populationRayCount == 0u)
            {
                result.boundaryStatus = data::BoundaryStatus::converged;
                break;
            }
            auto const relaunchFrame = getRayFrameSpec(populationRayCount, queue);
            auto const relaunchCandidateCount = boundaryCandidateCount(populationRayCount);
            auto& output = inputFirst ? scratch.second : scratch.first;
            auto& input = inputFirst ? scratch.first : scratch.second;
            alpaka::onHost::fill(
                queue,
                output.weights,
                0.0,
                alpaka::Vec{static_cast<std::size_t>(relaunchCandidateCount)});
            queue.enqueue(
                relaunchFrame,
                alpaka::KernelBundle{
                    hase::kernels::forward::AccumulateParticleCombedForwardPhiAse{},
                    hase::kernels::forward::TracePolicyList{
                        hase::kernels::forward::tracePolicy::source::boundaryCandidates,
                        hase::kernels::forward::tracePolicy::cell::forwardAse,
                        hase::kernels::forward::tracePolicy::boundary::boundaryCandidates,
                        hase::kernels::forward::tracePolicy::position::exact,
                        diagnostics},
                    mesh,
                    populationRayCount,
                    accumulation,
                    input.view(),
                    scratch.comb.selectedView(populationRayCount),
                    scratch.comb.selectedWeightsView(populationRayCount),
                    output.view(),
                    interfaceMap,
                    pass,
                    experiment.useReflections});
            routeCounts = measureRoutes(output, populationRayCount);
            auto const nextPopulationRayCount
                = std::accumulate(routeCounts.begin(), routeCounts.end(), std::uint32_t{0u});
            enqueueComb(output, relaunchCandidateCount, routeCounts, pass);
            double const currentWeight = readTotalWeight();
            result.boundaryPasses = pass;
            result.boundaryRemainingFraction = currentWeight / initialWeight;
            residualWeightFractions.push_back(result.boundaryRemainingFraction);
            if(currentWeight <= 0.0 || result.boundaryRemainingFraction < experiment.reflectionTolerance)
            {
                result.boundaryStatus = data::BoundaryStatus::converged;
                break;
            }
            if(currentWeight > previousWeight)
            {
                ++growCount;
                if(growCount >= divergenceStreak)
                {
                    auto const tail = estimateBoundaryTail(residualWeightFractions);
                    if(tail.divergent)
                    {
                        result.boundaryStatus = data::BoundaryStatus::diverged;
                        break;
                    }
                }
            }
            else
            {
                growCount = 0u;
            }
            previousWeight = currentWeight;
            populationRayCount = nextPopulationRayCount;
            inputFirst = !inputFirst;
        }

        auto const tail = estimateBoundaryTail(residualWeightFractions);
        result.boundaryGamma = tail.gamma;
        result.boundaryGammaStandardError = tail.gammaStandardError;
        result.boundaryTailFactor = tail.tailFactor;
        result.boundaryTailClosure = tail.tailClosure;
        if(result.boundaryStatus == data::BoundaryStatus::diverged || tail.divergent)
        {
            result.boundaryStatus = data::BoundaryStatus::diverged;
        }
        // A scalar residual fit is diagnostic only: it cannot establish a spatial field tail.
    }
} // namespace hase::core
