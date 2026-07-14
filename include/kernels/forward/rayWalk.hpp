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
    [[nodiscard]] inline ALPAKA_FN_ACC hase::core::Point advance(
        hase::core::Point const point,
        hase::core::Point const direction,
        double const length)
    {
        return point + direction * length;
    }

    [[nodiscard]] inline ALPAKA_FN_ACC hase::core::Point normalize(hase::core::Point const value)
    {
        double const length = value.euclidLength();
        if(length <= std::numeric_limits<double>::epsilon())
        {
            return hase::core::Point{0.0, 0.0, 0.0};
        }
        return value * (1.0 / length);
    }

    [[nodiscard]] inline ALPAKA_FN_ACC hase::core::Point faceCentroid(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        unsigned const localFace)
    {
        hase::core::Point sum{0.0, 0.0, 0.0};
        for(unsigned localVertex = 0u; localVertex < hase::core::tet4FaceWidth; ++localVertex)
        {
            int const point = mesh.getCellFacePoint(tet, localFace, localVertex);
            if(point < 0)
            {
                return sum;
            }
            sum = sum + mesh.getPoint(static_cast<unsigned>(point));
        }
        return sum * (1.0 / static_cast<double>(hase::core::tet4FaceWidth));
    }

    [[nodiscard]] inline ALPAKA_FN_ACC hase::core::Point outwardFaceNormal(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        unsigned const localFace)
    {
        int const p0 = mesh.getCellFacePoint(tet, localFace, 0u);
        int const p1 = mesh.getCellFacePoint(tet, localFace, 1u);
        int const p2 = mesh.getCellFacePoint(tet, localFace, 2u);
        if(p0 < 0 || p1 < 0 || p2 < 0)
        {
            return hase::core::Point{0.0, 0.0, 0.0};
        }
        hase::core::Point const a = mesh.getPoint(static_cast<unsigned>(p0));
        hase::core::Point const b = mesh.getPoint(static_cast<unsigned>(p1));
        hase::core::Point const c = mesh.getPoint(static_cast<unsigned>(p2));
        hase::core::Point normal = normalize(hase::core::cross(b - a, c - a));
        hase::core::Point const centroid = (a + b + c) * (1.0 / 3.0);
        if(hase::core::dot(normal, mesh.getCellCenterPoint(tet) - centroid) > 0.0)
        {
            normal = normal * -1.0;
        }
        return normal;
    }

    [[nodiscard]] inline ALPAKA_FN_ACC hase::core::Point reflectedDirection(
        hase::core::Point const direction,
        hase::core::Point const outwardNormal)
    {
        return normalize(direction - outwardNormal * (2.0 * hase::core::dot(direction, outwardNormal)));
    }

    [[nodiscard]] inline ALPAKA_FN_HOST_ACC double barycentricFaceIntersectionLength(
        double const coordinate,
        double const directionalChange,
        double const maxLength)
    {
        if(directionalChange >= -std::numeric_limits<double>::epsilon())
        {
            return 0.0;
        }
        double const length = -coordinate / directionalChange;
        constexpr double tolerance = 1.0e-10;
        return length > tolerance && length <= maxLength ? length : 0.0;
    }

    [[nodiscard]] inline ALPAKA_FN_ACC int nextFaceIntersection(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        hase::core::Point const origin,
        hase::core::Point const direction,
        int const forbiddenFace,
        double& length)
    {
        int nextFace = -1;
        // The first decreasing face coordinate to reach zero is the Tet4 exit face.
        for(unsigned localFace = 0u; localFace < mesh.numberOfFacesPerCell; ++localFace)
        {
            if(static_cast<int>(localFace) == forbiddenFace)
            {
                continue;
            }
            double const candidate = barycentricFaceIntersectionLength(
                mesh.getFaceBarycentricCoordinate(tet, localFace, origin),
                mesh.getFaceBarycentricDirection(tet, localFace, direction),
                length);
            if(candidate > 0.0)
            {
                length = candidate;
                nextFace = static_cast<int>(localFace);
            }
        }
        return nextFace;
    }

    [[nodiscard]] inline ALPAKA_FN_ACC double localGainCoefficient(
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

    [[nodiscard]] inline ALPAKA_FN_ACC double localSegmentGain(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        double const length,
        double const sigmaA,
        double const sigmaE)
    {
        return alpaka::math::exp(localGainCoefficient(mesh, tet, sigmaA, sigmaE) * length);
    }

    [[nodiscard]] inline ALPAKA_FN_ACC double localSegmentTrackLengthIntegral(
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

    [[nodiscard]] inline ALPAKA_FN_ACC double segmentCenterWeight(
        hase::core::DeviceMeshView const& mesh,
        unsigned const tet,
        hase::core::Point const midpoint)
    {
        return centerProximityWeight(barycentricCoordinates(mesh, tet, midpoint));
    }
} // namespace hase::kernels::forward
