/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/types.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace hase::core
{
    struct SimulationSnapshot
    {
        unsigned step = 0u;
        double time = 0.0;
        std::vector<double> betaVolume;
        Result aseResult;
        std::vector<double> dndtPump;
        std::vector<double> dndtAse;
        std::vector<std::string> fields;

        [[nodiscard]] bool contains(std::string const& field) const
        {
            return std::ranges::find(fields, field) != fields.end();
        }
    };
} // namespace hase::core
