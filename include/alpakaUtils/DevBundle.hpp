/**
 * Copyright 2026 Tim Hanel
 *
 */
#pragma once
#include <alpaka/alpaka.hpp>

template<alpaka::onHost::concepts::Device T_Device, typename T_Executor>
struct DevBundle
{
    T_Device device;
    T_Executor executor;

    DevBundle(T_Device const& device, T_Executor const& executor) : device(device), executor(executor)
    {
    }
};
