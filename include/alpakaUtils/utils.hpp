/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once
#include <alpaka/alpaka.hpp>

#include <utility>

namespace hase::alpakaUtils
{
    template<typename T>
    struct GetApiFromDevice;

    template<alpaka::concepts::Api T_Api, alpaka::concepts::DeviceKind T_DeviceKind>
    struct GetApiFromDevice<alpaka::onHost::Device<T_Api, T_DeviceKind>>
    {
        using type = T_Api;
    };

    using Vec1D = alpaka::Vec<uint32_t, 1>;
    using Vec2D = alpaka::Vec<uint32_t, 2>;
    using Vec3D = alpaka::Vec<uint32_t, 3>;

    ALPAKA_FN_ACC auto getLinGlobalIdx(alpaka::onAcc::concepts::Acc auto const& acc)
    {
        auto idxMd = acc.getIdxWithin(alpaka::onAcc::origin::grid, alpaka::onAcc::unit::threads);
        auto extentMd = acc.getExtentsOf(alpaka::onAcc::origin::grid, alpaka::onAcc::unit::threads);
        return alpaka::linearize(extentMd, idxMd);
    }

    template<typename T_Device>
    using ApiFromDevice = typename hase::alpakaUtils::GetApiFromDevice<T_Device>::type;

    template<typename T_Device>
    constexpr auto getApiFromDevice()
    {
        return std::declval<ApiFromDevice<T_Device>>();
    }
} // namespace hase::alpakaUtils
