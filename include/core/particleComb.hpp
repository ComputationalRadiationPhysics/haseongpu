#pragma once

#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <concepts/concepts.hpp>
#include <core/boundaryReduce.hpp>
#include <kernels/forward/particleCombing.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <stdexcept>
#include <utility>

namespace hase::core
{
    /** @brief Persistent device-only CDF and selection workspace for systematic particle combing. */
    template<alpaka::onHost::concepts::Device T_Device>
    class ParticleCombWorkspace
    {
        using T_DoubleBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<double>(std::declval<T_Device&>(), std::size_t{1u}));
        using T_UnsignedBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<std::uint32_t>(std::declval<T_Device&>(), std::size_t{1u}));
        using T_HistoryBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<std::uint64_t>(std::declval<T_Device&>(), std::size_t{1u}));
        using T_ByteBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<char>(std::declval<T_Device&>(), std::size_t{1u}));

    public:
        ParticleCombWorkspace(
            T_Device& device,
            std::uint32_t const maxCandidateCount,
            std::uint32_t const maxOutputCount,
            std::uint32_t const maxDomainCount = 1u)
            : cdf(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(maxCandidateCount)))
            , totalWeight(alpaka::onHost::alloc<double>(device, std::size_t{1u}))
            , systematicOffset(alpaka::onHost::alloc<double>(device, std::size_t{1u}))
            , filteredWeights(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(maxCandidateCount)))
            , selectedCandidates(
                  alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(maxOutputCount)))
            , selectedWeights(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(maxOutputCount)))
            , liveParents(alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(maxOutputCount)))
            , liveParentCount(alpaka::onHost::alloc<std::uint32_t>(device, std::size_t{1u}))
            , routeWeights(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(maxDomainCount)))
            , routeCandidateCounts(
                  alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(maxDomainCount)))
            , spatialKeys(alpaka::onHost::alloc<std::uint64_t>(device, static_cast<std::size_t>(maxCandidateCount)))
            , spatialKeysScratch(
                  alpaka::onHost::alloc<std::uint64_t>(device, static_cast<std::size_t>(maxCandidateCount)))
            , spatialIndices(alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(maxCandidateCount)))
            , spatialIndicesScratch(
                  alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(maxCandidateCount)))
            , radixFlags(alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(maxCandidateCount)))
            , radixPrefix(alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(maxCandidateCount)))
            , radixZeroCount(alpaka::onHost::alloc<std::uint32_t>(device, std::size_t{1u}))
            , scanBuffer(
                  alpaka::onHost::alloc<char>(
                      device,
                      alpaka::onHost::getScanBufferSize<double>(
                          alpaka::Vec{static_cast<std::size_t>(maxCandidateCount)})))
            , unsignedScanBuffer(
                  alpaka::onHost::alloc<char>(
                      device,
                      alpaka::onHost::getScanBufferSize<std::uint32_t>(
                          alpaka::Vec{static_cast<std::size_t>(maxCandidateCount)})))
            , m_maxCandidateCount(maxCandidateCount)
            , m_maxOutputCount(maxOutputCount)
            , m_maxDomainCount(maxDomainCount)
        {
        }

        /** @brief Reduce surviving parent count and target-domain weights without leaving the device. */
        template<alpaka::concepts::Executor T_Executor>
        void enqueueRouteMeasurements(
            alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
            concepts::Queue auto const& queue,
            alpaka::concepts::IView<double> auto candidateWeights,
            alpaka::concepts::IView<std::uint32_t> auto candidateDomains,
            std::uint32_t const parentCount,
            std::uint32_t const domainCount)
        {
            auto const candidateCount = 2u * parentCount;
            if(candidateCount > m_maxCandidateCount || parentCount > m_maxOutputCount
               || domainCount > m_maxDomainCount)
                throw std::out_of_range("boundary route measurement exceeds workspace capacity");
            auto const parentFrame = hase::alpakaUtils::getFrameSpec<std::uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{parentCount});
            queue.enqueue(
                parentFrame,
                alpaka::KernelBundle{
                    kernels::forward::MarkLiveBoundaryParents{},
                    candidateWeights.getSubView(alpaka::Vec{static_cast<std::size_t>(candidateCount)}),
                    liveParents.getView().getSubView(alpaka::Vec{static_cast<std::size_t>(parentCount)}),
                    std::uint32_t{parentCount}});
            detail::reduce(
                queue,
                std::uint32_t{0u},
                liveParentCount,
                liveParents.getView().getSubView(alpaka::Vec{static_cast<std::size_t>(parentCount)}));
            auto const candidateExtent = alpaka::Vec{static_cast<std::size_t>(candidateCount)};
            for(std::uint32_t domain = 0u; domain < domainCount; ++domain)
            {
                auto output = routeWeights.getView().getSubView(
                    alpaka::Vec{static_cast<std::size_t>(domain)},
                    alpaka::Vec{std::size_t{1u}});
                alpaka::onHost::transformReduce(
                    queue,
                    devBundle.executor,
                    0.0,
                    output,
                    std::plus{},
                    alpaka::ScalarFunc{kernels::forward::SelectBoundaryRouteWeight{domain}},
                    candidateWeights.getSubView(candidateExtent),
                    candidateDomains.getSubView(candidateExtent));
                auto countOutput = routeCandidateCounts.getView().getSubView(
                    alpaka::Vec{static_cast<std::size_t>(domain)},
                    alpaka::Vec{std::size_t{1u}});
                alpaka::onHost::transformReduce(
                    queue,
                    devBundle.executor,
                    std::uint32_t{0u},
                    countOutput,
                    std::plus{},
                    alpaka::ScalarFunc{kernels::forward::SelectBoundaryRouteCandidate{domain}},
                    candidateWeights.getSubView(candidateExtent),
                    candidateDomains.getSubView(candidateExtent));
            }
        }

        /** @brief Spatially stratify a destination route with unbiased stratum weights. */
        template<alpaka::concepts::Executor T_Executor, typename T_Candidates>
        void enqueueSpatialDomain(
            alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
            concepts::Queue auto const& queue,
            T_Candidates const candidates,
            std::uint32_t const candidateCount,
            std::uint32_t const routeCandidateCount,
            std::uint32_t const domain,
            std::uint32_t const outputOffset,
            std::uint32_t const outputCount,
            std::uint32_t const seed,
            std::uint64_t const historyId)
        {
            if(candidateCount > m_maxCandidateCount || routeCandidateCount > candidateCount
               || outputCount > routeCandidateCount || outputOffset > m_maxOutputCount
               || outputCount > m_maxOutputCount - outputOffset)
                throw std::out_of_range("spatial particle-comb route exceeds workspace capacity");
            if(outputCount == 0u)
                return;
            auto const candidateExtent = alpaka::Vec{static_cast<std::size_t>(candidateCount)};
            auto const routeExtent = alpaka::Vec{static_cast<std::size_t>(routeCandidateCount)};
            auto const candidateFrame = hase::alpakaUtils::getFrameSpec<std::uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{candidateCount});
            auto const routeFrame = hase::alpakaUtils::getFrameSpec<std::uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{routeCandidateCount});
            auto const scalarFrame = alpaka::onHost::FrameSpec{alpaka::Vec{1u}, alpaka::Vec{1u}, devBundle.executor};
            queue.enqueue(
                candidateFrame,
                alpaka::KernelBundle{
                    kernels::forward::MakeBoundarySpatialKeys{},
                    candidates.weights.getSubView(candidateExtent),
                    candidates.targetDomains.getSubView(candidateExtent),
                    candidates.targetCells.getSubView(candidateExtent),
                    candidates.targetFaces.getSubView(candidateExtent),
                    candidates.faceBarycentric.x.getSubView(candidateExtent),
                    candidates.faceBarycentric.y.getSubView(candidateExtent),
                    spatialKeys.getView().getSubView(candidateExtent),
                    spatialIndices.getView().getSubView(candidateExtent),
                    std::uint32_t{candidateCount},
                    std::uint32_t{domain}});
            for(std::uint32_t bit = 0u; bit < 64u; ++bit)
            {
                auto inputKeys
                    = (bit % 2u == 0u ? spatialKeys : spatialKeysScratch).getView().getSubView(candidateExtent);
                auto inputIndices
                    = (bit % 2u == 0u ? spatialIndices : spatialIndicesScratch).getView().getSubView(candidateExtent);
                auto outputKeys
                    = (bit % 2u == 0u ? spatialKeysScratch : spatialKeys).getView().getSubView(candidateExtent);
                auto outputIndices
                    = (bit % 2u == 0u ? spatialIndicesScratch : spatialIndices).getView().getSubView(candidateExtent);
                auto flags = radixFlags.getView().getSubView(candidateExtent);
                auto prefix = radixPrefix.getView().getSubView(candidateExtent);
                alpaka::onHost::transform(
                    queue,
                    devBundle.executor,
                    flags,
                    alpaka::ScalarFunc{kernels::forward::MarkBoundaryRadixZeros{bit}},
                    inputKeys);
                alpaka::onHost::exclusiveScan(queue, devBundle.executor, unsignedScanBuffer, prefix, flags);
                queue.enqueue(
                    scalarFrame,
                    alpaka::KernelBundle{
                        kernels::forward::CaptureBoundaryRadixZeroCount{},
                        ALPAKA_TYPEOF(flags){flags},
                        ALPAKA_TYPEOF(prefix){prefix},
                        radixZeroCount.getView(),
                        std::uint32_t{candidateCount}});
                queue.enqueue(
                    candidateFrame,
                    alpaka::KernelBundle{
                        kernels::forward::ScatterBoundaryRadixPass{},
                        ALPAKA_TYPEOF(inputKeys){inputKeys},
                        ALPAKA_TYPEOF(inputIndices){inputIndices},
                        ALPAKA_TYPEOF(flags){flags},
                        ALPAKA_TYPEOF(prefix){prefix},
                        radixZeroCount.getView(),
                        ALPAKA_TYPEOF(outputKeys){outputKeys},
                        ALPAKA_TYPEOF(outputIndices){outputIndices},
                        std::uint32_t{candidateCount}});
            }
            auto sortedIndices = spatialIndices.getView().getSubView(candidateExtent);
            auto sortedWeights = filteredWeights.getView().getSubView(routeExtent);
            queue.enqueue(
                routeFrame,
                alpaka::KernelBundle{
                    kernels::forward::GatherBoundaryRouteWeights{},
                    candidates.weights.getSubView(candidateExtent),
                    ALPAKA_TYPEOF(sortedIndices){sortedIndices},
                    ALPAKA_TYPEOF(sortedWeights){sortedWeights},
                    std::uint32_t{routeCandidateCount}});
            auto cdfView = cdf.getView().getSubView(routeExtent);
            alpaka::onHost::inclusiveScan(queue, devBundle.executor, scanBuffer, cdfView, sortedWeights);
            auto const outputFrame = hase::alpakaUtils::getFrameSpec<std::uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{outputCount});
            queue.enqueue(
                outputFrame,
                alpaka::KernelBundle{
                    kernels::forward::SelectSpatialBoundaryCandidates{},
                    ALPAKA_TYPEOF(cdfView){cdfView},
                    ALPAKA_TYPEOF(sortedIndices){sortedIndices},
                    selectedCandidates.getView(),
                    selectedWeights.getView(),
                    std::uint32_t{routeCandidateCount},
                    std::uint32_t{outputOffset},
                    std::uint32_t{outputCount},
                    std::uint32_t{seed},
                    std::uint64_t{historyId}});
        }

        /** @brief Enqueue scan, total capture, random offset, and candidate selection without waiting. */
        template<alpaka::concepts::Executor T_Executor>
        void enqueue(
            alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
            concepts::Queue auto const& queue,
            alpaka::concepts::IView<double> auto candidateWeights,
            std::uint32_t const candidateCount,
            std::uint32_t const outputCount,
            std::uint32_t const seed,
            std::uint64_t const historyId)
        {
            if(candidateCount > m_maxCandidateCount || outputCount > m_maxOutputCount)
                throw std::out_of_range("particle-comb workspace capacity exceeded");
            auto const extent = alpaka::Vec{static_cast<std::size_t>(candidateCount)};
            auto weightView = candidateWeights.getSubView(extent);
            enqueueSelection(devBundle, queue, weightView, candidateCount, 0u, outputCount, seed, historyId);
        }

        /** @brief Enqueue systematic combing for one target domain into a selected-output segment. */
        template<alpaka::concepts::Executor T_Executor>
        void enqueueDomain(
            alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
            concepts::Queue auto const& queue,
            alpaka::concepts::IView<double> auto candidateWeights,
            alpaka::concepts::IView<std::uint32_t> auto candidateDomains,
            std::uint32_t const candidateCount,
            std::uint32_t const domain,
            std::uint32_t const outputOffset,
            std::uint32_t const outputCount,
            std::uint32_t const seed,
            std::uint64_t const historyId)
        {
            if(candidateCount > m_maxCandidateCount || outputOffset > m_maxOutputCount
               || outputCount > m_maxOutputCount - outputOffset)
                throw std::out_of_range("particle-comb domain segment exceeds workspace capacity");
            auto const candidateExtent = alpaka::Vec{static_cast<std::size_t>(candidateCount)};
            auto const candidateFrame = hase::alpakaUtils::getFrameSpec<std::uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{candidateCount});
            auto filteredView = filteredWeights.getView().getSubView(candidateExtent);
            queue.enqueue(
                candidateFrame,
                alpaka::KernelBundle{
                    kernels::forward::FilterParticleCombDomainWeights{},
                    candidateWeights.getSubView(candidateExtent),
                    candidateDomains.getSubView(candidateExtent),
                    ALPAKA_TYPEOF(filteredView){filteredView},
                    std::uint32_t{candidateCount},
                    std::uint32_t{domain}});
            enqueueSelection(
                devBundle,
                queue,
                filteredView,
                candidateCount,
                outputOffset,
                outputCount,
                seed,
                historyId);
        }

        /** @brief Enqueue one domain segment when candidates are mesh faces. */
        template<alpaka::concepts::Executor T_Executor>
        void enqueueFaceDomain(
            alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
            concepts::Queue auto const& queue,
            alpaka::concepts::IView<double> auto faceWeights,
            std::span<std::uint32_t const> const cellDomains,
            std::uint32_t const faceCount,
            std::uint32_t const facesPerCell,
            std::uint32_t const domain,
            std::uint32_t const outputOffset,
            std::uint32_t const outputCount,
            std::uint32_t const seed,
            std::uint64_t const historyId)
        {
            if(faceCount > m_maxCandidateCount || outputOffset > m_maxOutputCount
               || outputCount > m_maxOutputCount - outputOffset)
                throw std::out_of_range("particle-comb face-domain segment exceeds workspace capacity");
            auto const candidateExtent = alpaka::Vec{static_cast<std::size_t>(faceCount)};
            auto const candidateFrame = hase::alpakaUtils::getFrameSpec<std::uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{faceCount});
            auto filteredView = filteredWeights.getView().getSubView(candidateExtent);
            queue.enqueue(
                candidateFrame,
                alpaka::KernelBundle{
                    kernels::forward::FilterParticleCombFaceDomainWeights{},
                    faceWeights.getSubView(candidateExtent),
                    cellDomains,
                    ALPAKA_TYPEOF(filteredView){filteredView},
                    std::uint32_t{faceCount},
                    std::uint32_t{facesPerCell},
                    std::uint32_t{domain}});
            enqueueSelection(devBundle, queue, filteredView, faceCount, outputOffset, outputCount, seed, historyId);
        }

        /** @brief Enqueue only the candidate scan and total capture used for convergence. */
        template<alpaka::concepts::Executor T_Executor>
        void enqueueTotal(
            alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
            concepts::Queue auto const& queue,
            alpaka::concepts::IView<double> auto candidateWeights,
            std::uint32_t const candidateCount)
        {
            if(candidateCount > m_maxCandidateCount)
                throw std::out_of_range("particle-comb candidate count exceeds workspace capacity");
            auto const extent = alpaka::Vec{static_cast<std::size_t>(candidateCount)};
            auto cdfView = cdf.getView().getSubView(extent);
            auto weightView = candidateWeights.getSubView(extent);
            alpaka::onHost::inclusiveScan(queue, devBundle.executor, scanBuffer, cdfView, weightView);
            auto const scalarFrame = alpaka::onHost::FrameSpec{alpaka::Vec{1u}, alpaka::Vec{1u}, devBundle.executor};
            queue.enqueue(
                scalarFrame,
                alpaka::KernelBundle{
                    kernels::forward::CaptureParticleCombTotalWeight{},
                    ALPAKA_TYPEOF(cdfView){cdfView},
                    totalWeight.getView(),
                    std::uint32_t{candidateCount}});
        }

        [[nodiscard]] auto selectedView(std::uint32_t const count)
        {
            return selectedCandidates.getView().getSubView(alpaka::Vec{static_cast<std::size_t>(count)});
        }

        [[nodiscard]] auto totalWeightView()
        {
            return totalWeight.getView();
        }

        [[nodiscard]] auto selectedWeightsView(std::uint32_t const count)
        {
            return selectedWeights.getView().getSubView(alpaka::Vec{static_cast<std::size_t>(count)});
        }

        T_DoubleBuffer cdf;
        T_DoubleBuffer totalWeight;
        T_DoubleBuffer systematicOffset;
        T_DoubleBuffer filteredWeights;
        T_UnsignedBuffer selectedCandidates;
        T_DoubleBuffer selectedWeights;
        T_UnsignedBuffer liveParents;
        T_UnsignedBuffer liveParentCount;
        T_DoubleBuffer routeWeights;
        T_UnsignedBuffer routeCandidateCounts;
        T_HistoryBuffer spatialKeys;
        T_HistoryBuffer spatialKeysScratch;
        T_UnsignedBuffer spatialIndices;
        T_UnsignedBuffer spatialIndicesScratch;
        T_UnsignedBuffer radixFlags;
        T_UnsignedBuffer radixPrefix;
        T_UnsignedBuffer radixZeroCount;
        T_ByteBuffer scanBuffer;
        T_ByteBuffer unsignedScanBuffer;

    private:
        template<alpaka::concepts::Executor T_Executor>
        void enqueueSelection(
            alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
            concepts::Queue auto const& queue,
            alpaka::concepts::IView<double> auto candidateWeights,
            std::uint32_t const candidateCount,
            std::uint32_t const outputOffset,
            std::uint32_t const outputCount,
            std::uint32_t const seed,
            std::uint64_t const historyId)
        {
            auto const extent = alpaka::Vec{static_cast<std::size_t>(candidateCount)};
            auto cdfView = cdf.getView().getSubView(extent);
            alpaka::onHost::inclusiveScan(queue, devBundle.executor, scanBuffer, cdfView, candidateWeights);
            auto const scalarFrame = alpaka::onHost::FrameSpec{alpaka::Vec{1u}, alpaka::Vec{1u}, devBundle.executor};
            queue.enqueue(
                scalarFrame,
                alpaka::KernelBundle{
                    kernels::forward::CaptureParticleCombTotalWeight{},
                    ALPAKA_TYPEOF(cdfView){cdfView},
                    totalWeight.getView(),
                    std::uint32_t{candidateCount}});
            queue.enqueue(
                scalarFrame,
                alpaka::KernelBundle{
                    kernels::forward::GenerateParticleCombOffset{},
                    totalWeight.getView(),
                    systematicOffset.getView(),
                    std::uint32_t{outputCount},
                    std::uint32_t{seed},
                    std::uint64_t{historyId}});
            auto const outputFrame = alpaka::onHost::FrameSpec{
                alpaka::Vec{(outputCount + 127u) / 128u},
                alpaka::Vec{128u},
                devBundle.executor};
            queue.enqueue(
                outputFrame,
                alpaka::KernelBundle{
                    kernels::forward::SelectParticleCombCandidates{},
                    ALPAKA_TYPEOF(cdfView){cdfView},
                    totalWeight.getView(),
                    systematicOffset.getView(),
                    selectedCandidates.getView(),
                    std::uint32_t{candidateCount},
                    std::uint32_t{outputCount},
                    std::uint32_t{outputOffset}});
            queue.enqueue(
                outputFrame,
                alpaka::KernelBundle{
                    kernels::forward::AssignParticleCombSelectedWeights{},
                    totalWeight.getView(),
                    selectedWeights.getView(),
                    std::uint32_t{outputCount},
                    std::uint32_t{outputOffset}});
        }

        std::uint32_t m_maxCandidateCount;
        std::uint32_t m_maxOutputCount;
        std::uint32_t m_maxDomainCount;
    };
} // namespace hase::core
