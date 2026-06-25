/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/mesh.hpp>

namespace hase::kernels::forward
{
    [[nodiscard]] ALPAKA_FN_ACC unsigned sampleVolumeByVolume(
        hase::core::DeviceMeshView const& mesh,
        double const totalVolume,
        alpaka::rand::engine::Philox4x32x10& rndEngine)
    {
        if(mesh.numberOfCells == 0u)
        {
            return 0u;
        }
        if(totalVolume <= 0.0)
        {
            auto const value = alpaka::rand::distribution::UniformReal{
                0.0f,
                static_cast<float>(mesh.numberOfCells),
                alpaka::rand::interval::oc}(rndEngine);
            unsigned const index = static_cast<unsigned>(value);
            return index < mesh.numberOfCells ? index : mesh.numberOfCells - 1u;
        }

        double const target = alpaka::rand::distribution::UniformReal<double>{}(rndEngine) * totalVolume;
        unsigned lower = 0u;
        unsigned upper = mesh.numberOfCells;
        while(lower < upper)
        {
            unsigned const middle = lower + (upper - lower) / 2u;
            if(target <= mesh.cellVolumePrefix[middle])
            {
                upper = middle;
            }
            else
            {
                lower = middle + 1u;
            }
        }
        return lower < mesh.numberOfCells ? lower : mesh.numberOfCells - 1u;
    }

    [[nodiscard]] ALPAKA_FN_ACC hase::core::Point samplePointInVolume(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        alpaka::rand::engine::Philox4x32x10& rndEngine)
    {
        return mesh.genRndPointInTetra(
            mesh.getCellPoint(tet, 0u),
            mesh.getCellPoint(tet, 1u),
            mesh.getCellPoint(tet, 2u),
            mesh.getCellPoint(tet, 3u),
            rndEngine);
    }

    [[nodiscard]] ALPAKA_FN_ACC hase::core::Point sampleIsotropicDirection(
        alpaka::rand::engine::Philox4x32x10& rndEngine)
    {
        constexpr double pi = 3.14159265358979323846;
        double const z = 2.0 * alpaka::rand::distribution::UniformReal<double>{}(rndEngine) - 1.0;
        double const phi = 2.0 * pi * alpaka::rand::distribution::UniformReal<double>{}(rndEngine);
        double const radius = alpaka::math::sqrt(alpaka::math::max(0.0, 1.0 - z * z));
        return hase::core::Point{radius * alpaka::math::cos(phi), radius * alpaka::math::sin(phi), z};
    }
} // namespace hase::kernels::forward
