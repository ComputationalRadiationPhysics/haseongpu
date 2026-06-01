/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once
#include <alpaka/alpaka.hpp>

namespace hase::alpakaUtils
{
    template<alpaka::onHost::concepts::Device T_Device, alpaka::concepts::Executor T_Executor>
    struct DevBundle
    {
        T_Device device;
        T_Executor executor;

        DevBundle(T_Device const& device, T_Executor const& executor) : device(device), executor(executor)
        {
        }
    };
} // namespace hase::alpakaUtils
