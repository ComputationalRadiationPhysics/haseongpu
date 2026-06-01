/**
 * Copyright 2026 Tim Hanel
 *
 */
#pragma once
#include <alpaka/alpaka.hpp>

#include <type_traits>

namespace hase::concepts
{
    template<typename T>
    struct IsQueue : std::false_type
    {
    };

    template<
        alpaka::concepts::Api T_Api,
        alpaka::concepts::DeviceKind T_DeviceKind,
        alpaka::concepts::QueueKind T_QueueKind>
    struct IsQueue<alpaka::onHost::Queue<alpaka::onHost::Device<T_Api, T_DeviceKind>, T_QueueKind>> : std::true_type
    {
    };

    template<typename T>
    concept Queue = IsQueue<std::remove_cvref_t<T>>::value;
} // namespace hase::concepts
