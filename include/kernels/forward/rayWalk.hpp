/**
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <core/mesh.hpp>
#include <kernels/forward/barycentric.hpp>

#include <limits>

namespace hase::kernels::forward
{
    [[nodiscard]] ALPAKA_FN_ACC hase::core::Point advance(
        hase::core::Point const point,
        hase::core::Point const direction,
        double const length)
    {
        return point + direction * length;
    }

    [[nodiscard]] ALPAKA_FN_ACC bool pointInTriangle(
        hase::core::Point const point,
        hase::core::Point const a,
        hase::core::Point const b,
        hase::core::Point const c)
    {
        hase::core::Point const v0 = c - a;
        hase::core::Point const v1 = b - a;
        hase::core::Point const v2 = point - a;
        double const dot00 = hase::core::dot(v0, v0);
        double const dot01 = hase::core::dot(v0, v1);
        double const dot02 = hase::core::dot(v0, v2);
        double const dot11 = hase::core::dot(v1, v1);
        double const dot12 = hase::core::dot(v1, v2);
        double const denominator = dot00 * dot11 - dot01 * dot01;
        if(alpaka::math::abs(denominator) <= std::numeric_limits<double>::epsilon())
        {
            return false;
        }
        double const invDenominator = 1.0 / denominator;
        double const u = (dot11 * dot02 - dot01 * dot12) * invDenominator;
        double const v = (dot00 * dot12 - dot01 * dot02) * invDenominator;
        constexpr double tolerance = 1.0e-10;
        return u >= -tolerance && v >= -tolerance && (u + v) <= 1.0 + tolerance;
    }

    [[nodiscard]] ALPAKA_FN_ACC bool pointInFace(
        hase::core::Point const point,
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        unsigned const localFace)
    {
        int const p0 = mesh.getCellFacePoint(tet, localFace, 0u);
        int const p1 = mesh.getCellFacePoint(tet, localFace, 1u);
        int const p2 = mesh.getCellFacePoint(tet, localFace, 2u);
        if(p0 < 0 || p1 < 0 || p2 < 0)
        {
            return false;
        }
        return pointInTriangle(
            point,
            mesh.getPoint(static_cast<unsigned>(p0)),
            mesh.getPoint(static_cast<unsigned>(p1)),
            mesh.getPoint(static_cast<unsigned>(p2)));
    }

    [[nodiscard]] ALPAKA_FN_ACC double faceIntersectionLength(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        unsigned const localFace,
        hase::core::Point const origin,
        hase::core::Point const direction,
        double const maxLength)
    {
        int const p0 = mesh.getCellFacePoint(tet, localFace, 0u);
        int const p1 = mesh.getCellFacePoint(tet, localFace, 1u);
        int const p2 = mesh.getCellFacePoint(tet, localFace, 2u);
        if(p0 < 0 || p1 < 0 || p2 < 0)
        {
            return 0.0;
        }

        hase::core::Point const a = mesh.getPoint(static_cast<unsigned>(p0));
        hase::core::Point const b = mesh.getPoint(static_cast<unsigned>(p1));
        hase::core::Point const c = mesh.getPoint(static_cast<unsigned>(p2));
        hase::core::Point const normal = hase::core::cross(b - a, c - a);
        double const denominator = hase::core::dot(normal, direction);
        if(alpaka::math::abs(denominator) <= std::numeric_limits<double>::epsilon())
        {
            return 0.0;
        }
        double const length = hase::core::dot(normal, a - origin) / denominator;
        constexpr double tolerance = 1.0e-10;
        if(length <= tolerance || length > maxLength)
        {
            return 0.0;
        }
        return pointInFace(advance(origin, direction, length), mesh, tet, localFace) ? length : 0.0;
    }

    [[nodiscard]] ALPAKA_FN_ACC int nextFaceIntersection(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        hase::core::Point const origin,
        hase::core::Point const direction,
        int const forbiddenFace,
        double& length)
    {
        int nextFace = -1;
        for(unsigned localFace = 0u; localFace < mesh.numberOfFacesPerCell; ++localFace)
        {
            if(static_cast<int>(localFace) == forbiddenFace)
            {
                continue;
            }
            double const candidate = faceIntersectionLength(mesh, tet, localFace, origin, direction, length);
            if(candidate > 0.0)
            {
                length = candidate;
                nextFace = static_cast<int>(localFace);
            }
        }
        return nextFace;
    }

    [[nodiscard]] ALPAKA_FN_ACC double localGainCoefficient(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        double const sigmaA,
        double const sigmaE)
    {
        if(mesh.getCellType(tet) == mesh.claddingNumber)
        {
            return -mesh.claddingAbsorption;
        }
        double const gainPerDensity = mesh.getBetaVolume(tet) * (sigmaE + sigmaA) - sigmaA;
        return static_cast<double>(mesh.nTot) * gainPerDensity;
    }

    [[nodiscard]] ALPAKA_FN_ACC double localSegmentGain(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        double const length,
        double const sigmaA,
        double const sigmaE)
    {
        return alpaka::math::exp(localGainCoefficient(mesh, tet, sigmaA, sigmaE) * length);
    }

    [[nodiscard]] ALPAKA_FN_ACC double localSegmentTrackLengthIntegral(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        double const length,
        double const sigmaA,
        double const sigmaE)
    {
        double const gainCoefficient = localGainCoefficient(mesh, tet, sigmaA, sigmaE);
        double const gainLength = gainCoefficient * length;
        if(alpaka::math::abs(gainLength) < 1.0e-8)
        {
            return length;
        }
        return (alpaka::math::exp(gainLength) - 1.0) / gainCoefficient;
    }

    [[nodiscard]] ALPAKA_FN_ACC double segmentCenterWeight(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        hase::core::Point const midpoint)
    {
        return centerProximityWeight(barycentricCoordinates(mesh, tet, midpoint));
    }
} // namespace hase::kernels::forward
