/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <atomic>
#include <stdexcept>

namespace hase::core
{
    class OperationCancelled : public std::runtime_error
    {
    public:
        OperationCancelled() : std::runtime_error("HASEonGPU operation cancelled")
        {
        }
    };

    inline std::atomic_bool& cancellationRequested()
    {
        static std::atomic_bool requested{false};
        return requested;
    }

    inline void requestCancellation()
    {
        cancellationRequested().store(true, std::memory_order_relaxed);
    }

    inline void clearCancellation()
    {
        cancellationRequested().store(false, std::memory_order_relaxed);
    }

    inline void throwIfCancellationRequested()
    {
        if(cancellationRequested().load(std::memory_order_relaxed))
        {
            throw OperationCancelled{};
        }
    }
} // namespace hase::core
