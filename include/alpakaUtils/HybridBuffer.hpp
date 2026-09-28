

/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <concepts/concepts.hpp>

#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace hase::internal::hybridBuffer
{
    /** @brief Host storage eligible for ownership transfer into a HybridBuffer. */
    template<typename T>
    concept RvalueStorage = !std::is_lvalue_reference_v<T>;

    /** @brief Storage that is not an Alpaka host-side device object. */
    template<typename T>
    concept HostStorage = !alpaka::onHost::concepts::Device<std::remove_cvref_t<T>>;

    /**
     * @param storage Host container or view whose lifetime remains external.
     * @return Non-owning Alpaka host view with matching data and extents.
     */
    template<HostStorage T_Storage>
    [[nodiscard]] auto makeHostView(T_Storage& storage)
    {
        if constexpr(alpaka::concepts::IView<T_Storage>)
            return alpaka::makeView(storage);
        else
            return alpaka::makeView(
                alpaka::api::host,
                alpaka::onHost::data(storage),
                alpaka::onHost::getExtents(storage));
    }

    /** @param storage Fixed-size host array. @return Non-owning one-dimensional host view. */
    template<typename T, std::size_t T_Size>
    [[nodiscard]] auto makeHostView(T (&storage)[T_Size])
    {
        return alpaka::makeView(alpaka::api::host, storage, alpaka::Vec{T_Size});
    }

    /**
     * @param hostBuffer Host representation to compare.
     * @param deviceBuffer Device representation to compare.
     * @throws std::invalid_argument If their runtime extents differ.
     */
    template<HostStorage T_HostBuffer, alpaka::concepts::IView T_DeviceBuffer>
    void validateMatchingBuffers(T_HostBuffer const& hostBuffer, T_DeviceBuffer const& deviceBuffer)
    {
        static_assert(
            alpaka::trait::getDim_v<T_HostBuffer> == alpaka::trait::getDim_v<T_DeviceBuffer>,
            "hybrid buffer host and device dimensionality must match");
        static_assert(
            std::same_as<
                std::remove_const_t<alpaka::GetValueType_t<T_HostBuffer>>,
                std::remove_const_t<alpaka::GetValueType_t<T_DeviceBuffer>>>,
            "hybrid buffer host and device value types must match");

        auto const hostExtents = alpaka::onHost::getExtents(hostBuffer);
        auto const deviceExtents = alpaka::onHost::getExtents(deviceBuffer);
        for(std::uint32_t dimension = 0u; dimension < alpaka::trait::getDim_v<T_HostBuffer>; ++dimension)
            if(hostExtents[dimension] != deviceExtents[dimension])
                throw std::invalid_argument("hybrid buffer host and device extents must match");
    }
} // namespace hase::internal::hybridBuffer

namespace hase::alpakaUtils
{
    /**
     * @brief One logical buffer with explicit host and device representations.
     *
     * Construction allocates or associates storage but never transfers data.
     * Lvalue storage is represented by a non-owning Alpaka view; rvalue host
     * storage is owned by the HybridBuffer. Non-owning storage must retain its
     * allocation and outlive the HybridBuffer and all queued operations using
     * its views.
     */
    template<hase::internal::hybridBuffer::HostStorage T_HostBuffer, alpaka::concepts::IView T_DeviceBuffer>
    class HybridBuffer
    {
    public:
        using HostBuffer = T_HostBuffer;
        using DeviceBuffer = T_DeviceBuffer;

        /**
         * @brief Associate host storage with a newly allocated device buffer.
         *
         * @param device Device on which matching storage is allocated.
         * @param hostBuffer Host representation to own or reference according
         * to `T_HostBuffer`.
         *
         * Construction does not transfer data.
         */
        template<alpaka::onHost::concepts::Device T_Device>
        HybridBuffer(T_Device const& device, T_HostBuffer hostBuffer)
            : m_hostBuffer(std::move(hostBuffer))
            , m_deviceBuffer(alpaka::onHost::allocLike(device, m_hostBuffer))
        {
        }

        /**
         * @brief Associate existing host and device representations without copying.
         *
         * @param hostBuffer Host representation to own or reference.
         * @param deviceBuffer Device view to associate.
         * @throws std::invalid_argument If the runtime extents differ.
         *
         * Element type and dimensionality mismatches are rejected at compile time.
         */
        HybridBuffer(T_HostBuffer hostBuffer, T_DeviceBuffer deviceBuffer)
            : m_hostBuffer(std::move(hostBuffer))
            , m_deviceBuffer(std::move(deviceBuffer))
        {
            hase::internal::hybridBuffer::validateMatchingBuffers(m_hostBuffer, m_deviceBuffer);
        }

        HybridBuffer(HybridBuffer const&) = delete;
        HybridBuffer& operator=(HybridBuffer const&) = delete;
        HybridBuffer(HybridBuffer&&) = default;
        HybridBuffer& operator=(HybridBuffer&&) = default;

        /**
         * @brief Enqueue a host-to-device copy.
         *
         * @param queue Queue targeting the associated device allocation.
         * @note This method does not wait. The buffer and referenced host
         * storage must remain valid until the queued copy completes.
         */
        void toDevice(concepts::Queue auto const& queue)
        {
            alpaka::onHost::memcpy(queue, m_deviceBuffer, m_hostBuffer);
        }

        /**
         * @brief Copy the device representation to the host representation.
         *
         * @param queue Queue targeting the associated device allocation.
         * @note This method waits for completion before returning.
         */
        void toHost(concepts::Queue auto const& queue)
        {
            alpaka::onHost::memcpy(queue, m_hostBuffer, m_deviceBuffer);
            alpaka::onHost::wait(queue);
        }

        /** @return Writable non-owning view of the host representation. */
        [[nodiscard]] auto getHostView()
        {
            return hase::internal::hybridBuffer::makeHostView(m_hostBuffer);
        }

        /** @return Non-owning host view obtained from a const logical buffer. */
        [[nodiscard]] auto getHostView() const
        {
            return hase::internal::hybridBuffer::makeHostView(m_hostBuffer);
        }

        /** @return Writable non-owning view of the device representation. */
        [[nodiscard]] auto toDeviceView()
        {
            return alpaka::makeView(m_deviceBuffer);
        }

        /** @return Non-owning device view obtained from a const logical buffer. */
        [[nodiscard]] auto toDeviceView() const
        {
            return alpaka::makeView(m_deviceBuffer);
        }

        /** @return Common host and device extents. */
        [[nodiscard]] auto getExtents() const
        {
            return alpaka::onHost::getExtents(m_hostBuffer);
        }

    private:
        T_HostBuffer m_hostBuffer;
        T_DeviceBuffer m_deviceBuffer;
    };

} // namespace hase::alpakaUtils

namespace hase::internal::hybridBuffer
{
    /**
     * @param device Device on which matching storage is allocated.
     * @param hostBuffer External host storage referenced by the result.
     * @return HybridBuffer with a non-owning host view and owning device allocation.
     */
    template<alpaka::onHost::concepts::Device T_Device, HostStorage T_HostBuffer>
    [[nodiscard]] auto makeHybridBuffer(T_Device const& device, T_HostBuffer& hostBuffer)
    {
        using HostBuffer = ALPAKA_TYPEOF(makeHostView(hostBuffer));
        using DeviceBuffer = ALPAKA_TYPEOF(
            alpaka::onHost::allocLike(std::declval<T_Device const&>(), std::declval<HostBuffer const&>()));
        return hase::alpakaUtils::HybridBuffer<HostBuffer, DeviceBuffer>(device, makeHostView(hostBuffer));
    }

    /**
     * @param device Device on which matching storage is allocated.
     * @param hostBuffer Rvalue host storage moved into the result.
     * @return HybridBuffer owning both host storage and its device allocation.
     */
    template<alpaka::onHost::concepts::Device T_Device, RvalueStorage T_HostBuffer>
    requires HostStorage<T_HostBuffer>
    [[nodiscard]] auto makeHybridBuffer(T_Device const& device, T_HostBuffer&& hostBuffer)
    {
        using HostBuffer = std::remove_cvref_t<T_HostBuffer>;
        using DeviceBuffer = ALPAKA_TYPEOF(
            alpaka::onHost::allocLike(std::declval<T_Device const&>(), std::declval<HostBuffer const&>()));
        return hase::alpakaUtils::HybridBuffer<HostBuffer, DeviceBuffer>(
            device,
            std::forward<T_HostBuffer>(hostBuffer));
    }

    /**
     * @param hostBuffer External host storage referenced by the result.
     * @param deviceBuffer Mutable external device allocation referenced by the result.
     * @return Non-owning HybridBuffer over both representations.
     */
    [[nodiscard]] auto makeHybridBuffer(HostStorage auto& hostBuffer, alpaka::concepts::IView auto& deviceBuffer)
    {
        using HostBuffer = ALPAKA_TYPEOF(makeHostView(hostBuffer));
        using DeviceBuffer = ALPAKA_TYPEOF(alpaka::makeView(deviceBuffer));
        return hase::alpakaUtils::HybridBuffer<HostBuffer, DeviceBuffer>(
            makeHostView(hostBuffer),
            alpaka::makeView(deviceBuffer));
    }

    /**
     * @param hostBuffer External host storage referenced by the result.
     * @param deviceBuffer Const external device view referenced by the result.
     * @return Non-owning HybridBuffer over both representations.
     */
    [[nodiscard]] auto makeHybridBuffer(HostStorage auto& hostBuffer, alpaka::concepts::IView auto const& deviceBuffer)
    {
        using HostBuffer = ALPAKA_TYPEOF(makeHostView(hostBuffer));
        using DeviceBuffer = ALPAKA_TYPEOF(alpaka::makeView(deviceBuffer));
        return hase::alpakaUtils::HybridBuffer<HostBuffer, DeviceBuffer>(
            makeHostView(hostBuffer),
            alpaka::makeView(deviceBuffer));
    }
} // namespace hase::internal::hybridBuffer

namespace hase::alpakaUtils
{

    /**
     * @brief Specializable operation used by `getHybridBuffer`.
     * @tparam T_First Decayed first dispatch argument.
     * @tparam T_Second Decayed second dispatch argument.
     */
    template<typename T_First, typename T_Second>
    struct GetHybridBuffer
    {
        using type = ALPAKA_TYPEOF(
            hase::internal::hybridBuffer::makeHybridBuffer(std::declval<T_First&>(), std::declval<T_Second&>()));

        template<typename T_FirstInput, typename T_SecondInput>
        requires(
            std::same_as<std::remove_cvref_t<T_FirstInput>, T_First>
            && std::same_as<std::remove_cvref_t<T_SecondInput>, T_Second>)
        /**
         * @param first Device or host storage accepted by the selected overload.
         * @param second Host storage or device view accepted by the selected overload.
         * @return HybridBuffer preserving the value category of both arguments.
         */
        [[nodiscard]] auto operator()(T_FirstInput&& first, T_SecondInput&& second) const
        {
            return hase::internal::hybridBuffer::makeHybridBuffer(
                std::forward<T_FirstInput>(first),
                std::forward<T_SecondInput>(second));
        }
    };

    template<typename T_First, typename T_Second>
    using GetHybridBuffer_t =
        typename GetHybridBuffer<std::remove_cvref_t<T_First>, std::remove_cvref_t<T_Second>>::type;

    /**
     * @brief Construct a HybridBuffer through the specializable GetHybridBuffer operation.
     *
     * @param first Device or host storage selected by the matching operation.
     * @param second Host storage or device view selected by the matching operation.
     * @return Logical buffer that owns moved host storage and otherwise retains
     * non-owning views of caller-provided storage.
     */
    template<typename T_First, typename T_Second>
    [[nodiscard]] auto getHybridBuffer(T_First&& first, T_Second&& second)
    {
        return GetHybridBuffer<std::remove_cvref_t<T_First>, std::remove_cvref_t<T_Second>>{}(
            std::forward<T_First>(first),
            std::forward<T_Second>(second));
    }
} // namespace hase::alpakaUtils
