/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/mesh.hpp>
#include <core/simulationRunControl.hpp>
#include <core/types.hpp>

namespace hase::core
{
    struct SimulationContext
    {
        core::ExperimentParameters experiment;
        core::ComputeParameters compute;
        core::HostMesh mesh;
        core::Result result;
        core::SimulationRunControl run;
    };
} // namespace hase::core
