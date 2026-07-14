/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/mesh.hpp>

#include <cstdint>

namespace hase::kernels::forward
{
    // Direct histories sample a spectral bin, a beta-volume-weighted source,
    // and an isotropic direction.  Stratification preserves those densities.
    // Cover the spectrum as evenly as possible; phase randomizes surplus bins.
    [[nodiscard]] inline ALPAKA_FN_ACC unsigned stratifiedSpectrumIndex(
        unsigned const spectrumSize,
        unsigned const globalRayIndex,
        unsigned const globalRayCount,
        unsigned const phase)
    {
        if(spectrumSize == 0u || globalRayCount == 0u)
        {
            return 0u;
        }
        unsigned const evenlySpacedIndex
            = static_cast<unsigned>(static_cast<std::uint64_t>(globalRayIndex) * spectrumSize / globalRayCount);
        return (evenlySpacedIndex + phase % spectrumSize) % spectrumSize;
    }

    [[nodiscard]] inline ALPAKA_FN_ACC double stratifiedUnitInterval(
        unsigned const globalRayIndex,
        unsigned const globalRayCount,
        double const shift)
    {
        return globalRayCount == 0u
                   ? 0.0
                   : (static_cast<double>(globalRayIndex) + shift) / static_cast<double>(globalRayCount);
    }

    [[nodiscard]] inline ALPAKA_FN_ACC unsigned sampleVolumeByVolume(
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

        double const target = alpaka::rand::distribution::UniformReal<double>{}(rndEngine) *totalVolume;
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

    [[nodiscard]] inline ALPAKA_FN_ACC unsigned sampleVolumeByBetaVolumeTarget(
        hase::core::DeviceMeshView const& mesh,
        double const betaVolumeTotal,
        double const target)
    {
        unsigned lower = 0u;
        unsigned upper = mesh.numberOfCells;
        while(lower < upper)
        {
            unsigned const middle = lower + (upper - lower) / 2u;
            if(target < mesh.betaVolumePrefix[middle])
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

    // Source probability is betaVolume[cell] * volume[cell] / total.
    [[nodiscard]] inline ALPAKA_FN_ACC unsigned sampleVolumeByBetaVolume(
        hase::core::DeviceMeshView const& mesh,
        double const betaVolumeTotal,
        alpaka::rand::engine::Philox4x32x10& rndEngine)
    {
        if(mesh.numberOfCells == 0u || betaVolumeTotal <= 0.0 || mesh.betaVolumePrefix.size() != mesh.numberOfCells)
        {
            return sampleVolumeByVolume(
                mesh,
                mesh.cellVolumePrefix.empty() ? 0.0 : mesh.cellVolumePrefix.back(),
                rndEngine);
        }

        double const target = alpaka::rand::distribution::UniformReal<double>{}(rndEngine) *betaVolumeTotal;
        return sampleVolumeByBetaVolumeTarget(mesh, betaVolumeTotal, target);
    }

    // One randomized systematic CDF point per globally assigned source ray.
    [[nodiscard]] inline ALPAKA_FN_ACC unsigned sampleStratifiedVolumeByBetaVolume(
        hase::core::DeviceMeshView const& mesh,
        double const betaVolumeTotal,
        unsigned const globalRayIndex,
        unsigned const globalRayCount,
        double const shift,
        alpaka::rand::engine::Philox4x32x10& rndEngine)
    {
        if(mesh.numberOfCells == 0u || betaVolumeTotal <= 0.0 || mesh.betaVolumePrefix.size() != mesh.numberOfCells
           || globalRayCount == 0u)
        {
            return sampleVolumeByBetaVolume(mesh, betaVolumeTotal, rndEngine);
        }
        return sampleVolumeByBetaVolumeTarget(
            mesh,
            betaVolumeTotal,
            stratifiedUnitInterval(globalRayIndex, globalRayCount, shift) * betaVolumeTotal);
    }

    [[nodiscard]] inline ALPAKA_FN_ACC hase::core::Point samplePointInVolume(
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

    // Geometry-independent physical source distribution; no inferred axis.
    [[nodiscard]] inline ALPAKA_FN_ACC hase::core::Point sampleIsotropicDirection(
        alpaka::rand::engine::Philox4x32x10& rndEngine)
    {
        constexpr double pi = 3.14159265358979323846;
        double const z = 2.0 * alpaka::rand::distribution::UniformReal<double>{}(rndEngine) -1.0;
        double const phi = 2.0 * pi * alpaka::rand::distribution::UniformReal<double>{}(rndEngine);
        double const radius = alpaka::math::sqrt(alpaka::math::max(0.0, 1.0 - z * z));
        return hase::core::Point{radius * alpaka::math::cos(phi), radius * alpaka::math::sin(phi), z};
    }
} // namespace hase::kernels::forward
