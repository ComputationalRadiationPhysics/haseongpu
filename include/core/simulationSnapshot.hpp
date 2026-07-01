/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/mesh.hpp>
#include <core/types.hpp>

#include <vector>

namespace hase::core
{
    struct SimulationSnapshot
    {
        unsigned step = 0u;
        double time = 0.0;
        HostMesh mesh;
        Result aseResult;
        std::vector<double> dndtPump;
        std::vector<double> dndtAse;
    };
} // namespace hase::core
