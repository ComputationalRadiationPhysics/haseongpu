/**
 * Copyright 2026 Tim Hanel
 *
 */
#pragma once
#include <alpaka/alpaka.hpp>
using Vec1D = alpaka::Vec<uint32_t, 1>;
using Vec2D = alpaka::Vec<uint32_t, 2>;
using Vec3D = alpaka::Vec<uint32_t, 3>;

namespace hase::alpakaUtils
{
    ALPAKA_FN_ACC auto getLinGlobalIdx(alpaka::onAcc::concepts::Acc auto const& acc)
    {
        auto idxMd = acc.getIdxWithin(alpaka::onAcc::origin::grid, alpaka::onAcc::unit::threads);
        auto extentMd = acc.getExtentsOf(alpaka::onAcc::origin::grid, alpaka::onAcc::unit::threads);
        return alpaka::linearize(extentMd, idxMd);
    }

    namespace internal
    {
        template<typename T>
        struct GetApiFromDevice;

        template<alpaka::concepts::Api T_Api, alpaka::concepts::DeviceKind T_DeviceKind>
        struct GetApiFromDevice<alpaka::onHost::Device<T_Api, T_DeviceKind>>
        {
            using type = T_Api;
        };
    } // namespace internal

    template<typename T_Device>
    using GetApiFromDevice = typename internal::GetApiFromDevice<T_Device>::type;

    template<typename T_Device>
    constexpr auto getApiFromDevice()
    {
        return std::declval<GetApiFromDevice<T_Device>>();
    }
} // namespace hase::alpakaUtils
