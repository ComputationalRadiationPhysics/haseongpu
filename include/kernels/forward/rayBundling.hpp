/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <alpaka/alpaka.hpp>

#include <core/geometry.hpp>

#include <cstdint>

namespace hase::kernels::forward
{
    // Two equal-area z bands and four equal-width azimuth sectors.  Keeping the
    // number of strata small is intentional: a (cell, direction) bucket still
    // needs enough histories to fill a useful block-sized bundle.
    inline constexpr std::uint32_t forwardDirectionZStrata = 2u;
    inline constexpr std::uint32_t forwardDirectionPhiStrata = 4u;
    inline constexpr std::uint32_t forwardDirectionStrata = forwardDirectionZStrata * forwardDirectionPhiStrata;

    [[nodiscard]] inline ALPAKA_FN_HOST_ACC unsigned directionStratum(hase::core::Point const direction)
    {
        unsigned const zStratum = direction.z >= 0.0 ? 1u : 0u;
        unsigned phiStratum = 0u;
        if(direction.y >= 0.0)
        {
            phiStratum = direction.x >= 0.0 ? 0u : 1u;
        }
        else
        {
            phiStratum = direction.x < 0.0 ? 2u : 3u;
        }
        return zStratum * forwardDirectionPhiStrata + phiStratum;
    }

} // namespace hase::kernels::forward
