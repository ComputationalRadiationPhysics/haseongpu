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

#include <array>
#include <cmath>
#include <limits>

namespace hase::kernels::forward
{
    using BarycentricTet4 = std::array<double, 4u>;

    [[nodiscard]] inline ALPAKA_FN_HOST_ACC double signedTetVolume6(
        hase::core::Point const a,
        hase::core::Point const b,
        hase::core::Point const c,
        hase::core::Point const d)
    {
        return hase::core::dot(hase::core::cross(b - a, c - a), d - a);
    }

    [[nodiscard]] inline ALPAKA_FN_HOST_ACC BarycentricTet4 barycentricCoordinates(
        hase::core::Point const point,
        hase::core::Point const a,
        hase::core::Point const b,
        hase::core::Point const c,
        hase::core::Point const d)
    {
        double const denominator = signedTetVolume6(a, b, c, d);
        if(alpaka::math::abs(denominator) <= std::numeric_limits<double>::epsilon())
        {
            return {0.25, 0.25, 0.25, 0.25};
        }

        double const invDenominator = 1.0 / denominator;
        return {
            signedTetVolume6(point, b, c, d) * invDenominator,
            signedTetVolume6(a, point, c, d) * invDenominator,
            signedTetVolume6(a, b, point, d) * invDenominator,
            signedTetVolume6(a, b, c, point) * invDenominator};
    }

    [[nodiscard]] inline ALPAKA_FN_HOST_ACC BarycentricTet4 barycentricCoordinates(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        hase::core::Point const point)
    {
        return barycentricCoordinates(
            point,
            mesh.getCellPoint(tet, 0u),
            mesh.getCellPoint(tet, 1u),
            mesh.getCellPoint(tet, 2u),
            mesh.getCellPoint(tet, 3u));
    }

    [[nodiscard]] inline ALPAKA_FN_HOST_ACC double centerProximityWeight(BarycentricTet4 const& barycentric)
    {
        double distanceSquared = 0.0;
        for(double const coordinate : barycentric)
        {
            double const delta = coordinate - 0.25;
            distanceSquared += delta * delta;
        }
        constexpr double maxCenterDistance = 0.86602540378443864676;
        double const normalizedDistance = alpaka::math::sqrt(distanceSquared) / maxCenterDistance;
        return alpaka::math::max(0.0, 1.0 - normalizedDistance);
    }
} // namespace hase::kernels::forward
