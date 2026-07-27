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
#include <random/random.hpp>

namespace hase::core
{
    inline double calcVolumeDndtAse(
        HostMesh const& mesh,
        double const sigmaA,
        double const sigmaE,
        float const phiAse,
        unsigned const volume)
    {
        double const gainPerDensity = mesh.betaVolume[volume] * (sigmaE + sigmaA) - sigmaA;
        return gainPerDensity * phiAse;
    }

    inline unsigned baseRngSeed(ComputeParameters const& compute)
    {
        if(compute.rngSeed == ComputeParameters::unspecifiedRngSeed)
            return random::SeedGenerator::get().getSeed();
        return compute.rngSeed;
    }
} // namespace hase::core
