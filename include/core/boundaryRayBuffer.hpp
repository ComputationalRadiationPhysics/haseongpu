#pragma once

#include <alpaka/alpaka.hpp>

#include <alpakaUtils/DevBundle.hpp>
#include <alpakaUtils/HybridBuffer.hpp>
#include <alpakaUtils/memory.hpp>
#include <concepts/concepts.hpp>
#include <core/geometry.hpp>
#include <data/AseDomainGraph.hpp>
#include <kernels/forwardPhiAseMapping.hpp>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace hase::core
{
    /** @return Candidate capacity required for deterministic reflection/transmission splitting. */
    [[nodiscard]] inline std::uint32_t boundaryCandidateCount(std::uint32_t const rayCount)
    {
        if(rayCount > std::numeric_limits<std::uint32_t>::max() / 2u)
            throw std::overflow_error("boundary candidate count exceeds the 32-bit index range");
        return 2u * rayCount;
    }

    /** @brief Persistent device copy of one domain's compact face-routing table. */
    template<alpaka::onHost::concepts::Device T_Device>
    class ResidentAseDomainInterfaces
    {
    public:
        ResidentAseDomainInterfaces(
            T_Device device,
            std::vector<hase::data::DomainId>& targetDomainValues,
            std::vector<std::uint32_t>& targetCellValues,
            std::vector<std::uint32_t>& targetFaceValues,
            std::vector<hase::data::DomainId>& cellDomainValues,
            std::vector<float>& reflectivityValues,
            std::vector<float>& sourceRefractiveIndexValues,
            std::vector<float>& targetRefractiveIndexValues)
            : m_device(std::move(device))
            , targetDomains(hase::alpakaUtils::getHybridBuffer(m_device, targetDomainValues))
            , targetCells(hase::alpakaUtils::getHybridBuffer(m_device, targetCellValues))
            , targetFaces(hase::alpakaUtils::getHybridBuffer(m_device, targetFaceValues))
            , cellDomains(hase::alpakaUtils::getHybridBuffer(m_device, cellDomainValues))
            , reflectivities(hase::alpakaUtils::getHybridBuffer(m_device, reflectivityValues))
            , sourceRefractiveIndices(hase::alpakaUtils::getHybridBuffer(m_device, sourceRefractiveIndexValues))
            , targetRefractiveIndices(hase::alpakaUtils::getHybridBuffer(m_device, targetRefractiveIndexValues))
        {
        }

        ResidentAseDomainInterfaces(T_Device device, hase::data::AseDomain& domain)
            : ResidentAseDomainInterfaces(
                  std::move(device),
                  domain.boundaryTargetDomains,
                  domain.boundaryTargetCells,
                  domain.boundaryTargetFaces,
                  domain.cellDomains,
                  domain.boundaryReflectivities,
                  domain.boundarySourceRefractiveIndices,
                  domain.boundaryTargetRefractiveIndices)
        {
        }

        ResidentAseDomainInterfaces(T_Device device, hase::data::AseDomainGraph& graph)
            : ResidentAseDomainInterfaces(
                  std::move(device),
                  graph.globalBoundaryTargetDomains,
                  graph.globalBoundaryTargetCells,
                  graph.globalBoundaryTargetFaces,
                  graph.globalCellDomains,
                  graph.globalBoundaryReflectivities,
                  graph.globalBoundarySourceRefractiveIndices,
                  graph.globalBoundaryTargetRefractiveIndices)
        {
        }

        /** @brief Enqueue all routing-table uploads without synchronizing the queue. */
        void toDevice(concepts::Queue auto const& queue)
        {
            targetDomains.toDevice(queue);
            targetCells.toDevice(queue);
            targetFaces.toDevice(queue);
            cellDomains.toDevice(queue);
            reflectivities.toDevice(queue);
            sourceRefractiveIndices.toDevice(queue);
            targetRefractiveIndices.toDevice(queue);
        }

        [[nodiscard]] hase::data::AseDomainInterfaceView view() const
        {
            return {
                std::span<hase::data::DomainId const>(
                    targetDomains.toDeviceView().data(),
                    targetDomains.getExtents().x()),
                std::span<std::uint32_t const>(targetCells.toDeviceView().data(), targetCells.getExtents().x()),
                std::span<std::uint32_t const>(targetFaces.toDeviceView().data(), targetFaces.getExtents().x()),
                std::span<hase::data::DomainId const>(cellDomains.toDeviceView().data(), cellDomains.getExtents().x()),
                std::span<float const>(reflectivities.toDeviceView().data(), reflectivities.getExtents().x()),
                std::span<float const>(
                    sourceRefractiveIndices.toDeviceView().data(),
                    sourceRefractiveIndices.getExtents().x()),
                std::span<float const>(
                    targetRefractiveIndices.toDeviceView().data(),
                    targetRefractiveIndices.getExtents().x())};
        }

    private:
        T_Device m_device;

    public:
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<hase::data::DomainId>> targetDomains;
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<std::uint32_t>> targetCells;
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<std::uint32_t>> targetFaces;
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<hase::data::DomainId>> cellDomains;
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<float>> reflectivities;
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<float>> sourceRefractiveIndices;
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<float>> targetRefractiveIndices;
    };

    /** @brief Persistent device copy of precomputed per-domain source sampling tables. */
    template<alpaka::onHost::concepts::Device T_Device>
    class ResidentAseDomainSources
    {
        using T_ByteBuffer = ALPAKA_TYPEOF(alpaka::onHost::alloc<char>(std::declval<T_Device&>(), std::size_t{1u}));

    public:
        ResidentAseDomainSources(T_Device device, hase::data::AseDomainGraph& graph)
            : m_device(std::move(device))
            , offsets(hase::alpakaUtils::getHybridBuffer(m_device, graph.domainCellOffsets))
            , globalCells(hase::alpakaUtils::getHybridBuffer(m_device, graph.domainGlobalCells))
            , sourceStrengthPrefix(hase::alpakaUtils::getHybridBuffer(m_device, graph.domainSourceStrengthPrefix))
            , sourceStrengthTotals(hase::alpakaUtils::getHybridBuffer(m_device, graph.domainSourceStrengthTotals))
            , scanBuffer(
                  alpaka::onHost::alloc<char>(
                      m_device,
                      alpaka::onHost::getScanBufferSize<double>(alpaka::Vec{graph.domainSourceStrengthPrefix.size()})))
        {
        }

        void toDevice(concepts::Queue auto const& queue)
        {
            offsets.toDevice(queue);
            globalCells.toDevice(queue);
            sourceStrengthPrefix.toDevice(queue);
            sourceStrengthTotals.toDevice(queue);
        }

        /** @brief Rebuild independent domain-local CDFs in one flattened device buffer. */
        template<alpaka::concepts::Executor T_Executor>
        void rebuild(
            alpakaUtils::DevBundle<T_Device, T_Executor>& devBundle,
            concepts::Queue auto const& queue,
            hase::data::TraceView const mesh)
        {
            auto const cellCount = static_cast<std::uint32_t>(globalCells.getExtents().x());
            auto const domainCount = static_cast<std::uint32_t>(sourceStrengthTotals.getExtents().x());
            if(cellCount == 0u)
            {
                auto totals = sourceStrengthTotals.toDeviceView();
                alpaka::onHost::fill(queue, totals, 0.0);
                return;
            }
            auto sourceStrengthPrefixView = sourceStrengthPrefix.toDeviceView();
            auto globalCellsView = globalCells.toDeviceView();
            alpaka::onHost::transform(
                queue,
                devBundle.executor,
                sourceStrengthPrefixView,
                hase::kernels::DomainSourceStrength{mesh},
                globalCellsView);
            // Topology offsets are already available on the host. Only scans execute
            // here; source values and prefixes never leave the device. Starting each
            // CDF at zero preserves weak sources next to much stronger domains.
            auto const hostOffsets = offsets.getHostView();
            for(std::uint32_t domain = 0u; domain < domainCount; ++domain)
            {
                auto const begin = static_cast<std::size_t>(hostOffsets[domain]);
                auto const end = static_cast<std::size_t>(hostOffsets[domain + 1u]);
                if(begin != end)
                {
                    auto localPrefix
                        = sourceStrengthPrefixView.getSubView(alpaka::Vec{begin}, alpaka::Vec{end - begin});
                    alpaka::onHost::inclusiveScanInPlace(queue, devBundle.executor, scanBuffer, localPrefix);
                }
            }
            auto const domainFrame = hase::alpakaUtils::getFrameSpec<std::uint32_t>(
                devBundle.device,
                devBundle.executor,
                alpaka::Vec{domainCount});
            queue.enqueue(
                domainFrame,
                alpaka::KernelBundle{
                    hase::kernels::ComputeDomainSourceStrengthTotals{},
                    offsets.toDeviceView(),
                    sourceStrengthPrefix.toDeviceView(),
                    sourceStrengthTotals.toDeviceView(),
                    domainCount});
        }

        /** @return Current per-domain source totals after one small synchronized download. */
        [[nodiscard]] std::vector<double> downloadSourceStrengthTotals(concepts::Queue auto const& queue)
        {
            sourceStrengthTotals.toHost(queue);
            auto const host = sourceStrengthTotals.getHostView();
            std::vector<double> result(host.getExtents().product());
            for(std::size_t index = 0u; index < result.size(); ++index)
                result[index] = host[index];
            return result;
        }

        [[nodiscard]] hase::data::AseDomainSourceView view() const
        {
            return {
                std::span<std::uint32_t const>(offsets.toDeviceView().data(), offsets.getExtents().x()),
                std::span<std::uint32_t const>(globalCells.toDeviceView().data(), globalCells.getExtents().x()),
                std::span<double const>(
                    sourceStrengthPrefix.toDeviceView().data(),
                    sourceStrengthPrefix.getExtents().x()),
                std::span<double const>(
                    sourceStrengthTotals.toDeviceView().data(),
                    sourceStrengthTotals.getExtents().x())};
        }

    private:
        T_Device m_device;

    public:
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<std::uint32_t>> offsets;
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<std::uint32_t>> globalCells;
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<double>> sourceStrengthPrefix;
        hase::alpakaUtils::GetHybridBuffer_t<T_Device, std::vector<double>> sourceStrengthTotals;
        T_ByteBuffer scanBuffer;
    };

    /** @brief Trivially-copyable device views of compact domain-boundary histories. */
    template<
        alpaka::concepts::IView<double> T_DoubleView,
        alpaka::concepts::IView<std::uint32_t> T_UnsignedView,
        alpaka::concepts::IView<std::uint64_t> T_HistoryView>
    struct BoundaryRaySpans
    {
        PositionViewSoA<T_DoubleView> positions;
        DirectionViewSoA<T_DoubleView> directions;
        PositionViewSoA<T_DoubleView> faceBarycentric;
        T_DoubleView weights;
        T_DoubleView wavelengths;
        T_UnsignedView targetDomains;
        T_UnsignedView targetCells;
        T_UnsignedView targetFaces;
        T_UnsignedView rayPopulationIds;
        T_UnsignedView reflectionDepths;
        T_HistoryView historyIds;
    };

    /** @brief Owning SoA device buffer reused as a domain input, candidate, or route queue. */
    template<alpaka::onHost::concepts::Device T_Device>
    class BoundaryRayBuffer
    {
        using T_DoubleBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<double>(std::declval<T_Device&>(), std::size_t{1u}));
        using T_UnsignedBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<std::uint32_t>(std::declval<T_Device&>(), std::size_t{1u}));
        using T_HistoryBuffer
            = ALPAKA_TYPEOF(alpaka::onHost::alloc<std::uint64_t>(std::declval<T_Device&>(), std::size_t{1u}));

    public:
        BoundaryRayBuffer(T_Device& device, std::uint32_t const capacity)
            : positions(device, static_cast<std::size_t>(capacity))
            , directions(device, static_cast<std::size_t>(capacity))
            , faceBarycentric(device, static_cast<std::size_t>(capacity))
            , weights(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(capacity)))
            , wavelengths(alpaka::onHost::alloc<double>(device, static_cast<std::size_t>(capacity)))
            , targetDomains(alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(capacity)))
            , targetCells(alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(capacity)))
            , targetFaces(alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(capacity)))
            , rayPopulationIds(alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(capacity)))
            , reflectionDepths(alpaka::onHost::alloc<std::uint32_t>(device, static_cast<std::size_t>(capacity)))
            , historyIds(alpaka::onHost::alloc<std::uint64_t>(device, static_cast<std::size_t>(capacity)))
            , m_capacity(capacity)
        {
        }

        [[nodiscard]] auto view()
        {
            auto result = BoundaryRaySpans{
                positions.view(),
                directions.view(),
                faceBarycentric.view(),
                weights.getView(),
                wavelengths.getView(),
                targetDomains.getView(),
                targetCells.getView(),
                targetFaces.getView(),
                rayPopulationIds.getView(),
                reflectionDepths.getView(),
                historyIds.getView()};
            static_assert(std::is_trivially_copyable_v<ALPAKA_TYPEOF(result)>);
            static_assert(alpaka::concepts::KernelArg<ALPAKA_TYPEOF(result)>);
            return result;
        }

        [[nodiscard]] std::uint32_t capacity() const
        {
            return m_capacity;
        }

        PositionBufferSoA<T_Device> positions;
        DirectionBufferSoA<T_Device> directions;
        PositionBufferSoA<T_Device> faceBarycentric;
        T_DoubleBuffer weights;
        T_DoubleBuffer wavelengths;
        T_UnsignedBuffer targetDomains;
        T_UnsignedBuffer targetCells;
        T_UnsignedBuffer targetFaces;
        T_UnsignedBuffer rayPopulationIds;
        T_UnsignedBuffer reflectionDepths;
        T_HistoryBuffer historyIds;

    private:
        std::uint32_t m_capacity;
    };
} // namespace hase::core
