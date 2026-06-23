/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 * Copyright 2026 Tim Hanel
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */
#include <kernels/propagateRay.hpp>
#include <cassert> /* assert */
#include <cstdio>

namespace hase::kernels
{
    ALPAKA_FN_ACC hase::core::Ray calcNextRay(hase::core::Ray ray, double const length)
    {
        ray.p.x = ray.p.x + length * ray.dir.x;
        ray.p.y = ray.p.y + length * ray.dir.y;
        ray.p.z = ray.p.z + length * ray.dir.z;
        return ray;
    }

    ALPAKA_FN_ACC bool pointInTriangle3D(
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

    ALPAKA_FN_ACC bool pointInFace(
        hase::core::Point const point,
        hase::core::DeviceMeshView const& mesh,
        unsigned const cell,
        unsigned const localFace)
    {
        int const p0 = mesh.getCellFacePoint(cell, localFace, 0u);
        int const p1 = mesh.getCellFacePoint(cell, localFace, 1u);
        int const p2 = mesh.getCellFacePoint(cell, localFace, 2u);
        int const p3 = mesh.getCellFacePoint(cell, localFace, 3u);
        if(p0 < 0 || p1 < 0 || p2 < 0)
        {
            return false;
        }
        hase::core::Point const a = mesh.getPoint(static_cast<unsigned>(p0));
        hase::core::Point const b = mesh.getPoint(static_cast<unsigned>(p1));
        hase::core::Point const c = mesh.getPoint(static_cast<unsigned>(p2));
        if(pointInTriangle3D(point, a, b, c))
        {
            return true;
        }
        if(p3 < 0)
        {
            return false;
        }
        hase::core::Point const d = mesh.getPoint(static_cast<unsigned>(p3));
        return pointInTriangle3D(point, a, c, d);
    }

    ALPAKA_FN_ACC double checkCellFace(
        unsigned const cell,
        unsigned const localFace,
        hase::core::Ray const ray,
        hase::core::DeviceMeshView const& mesh,
        double const length)
    {
        int const p0 = mesh.getCellFacePoint(cell, localFace, 0u);
        int const p1 = mesh.getCellFacePoint(cell, localFace, 1u);
        int const p2 = mesh.getCellFacePoint(cell, localFace, 2u);
        if(p0 < 0 || p1 < 0 || p2 < 0)
        {
            return 0.0;
        }
        hase::core::Point const a = mesh.getPoint(static_cast<unsigned>(p0));
        hase::core::Point const b = mesh.getPoint(static_cast<unsigned>(p1));
        hase::core::Point const c = mesh.getPoint(static_cast<unsigned>(p2));
        hase::core::Point const normal = hase::core::cross(b - a, c - a);
        double const denominator = hase::core::dot(normal, ray.dir);
        if(alpaka::math::abs(denominator) <= std::numeric_limits<double>::epsilon())
        {
            return 0.0;
        }
        double const lengthTmp = hase::core::dot(normal, a - ray.p) / denominator;
        constexpr double tolerance = 1.0e-10;
        if(lengthTmp <= tolerance || lengthTmp > length)
        {
            return 0.0;
        }
        hase::core::Point const intersection = calcNextRay(ray, lengthTmp).p;
        return pointInFace(intersection, mesh, cell, localFace) ? lengthTmp : 0.0;
    }

    ALPAKA_FN_ACC int calcCellRayIntersection(
        double* length,
        unsigned const cell,
        hase::core::Ray const ray,
        int const forbiddenFace,
        hase::core::DeviceMeshView const& mesh)
    {
        int nextFace = -1;
        for(unsigned localFace = 0u; localFace < mesh.numberOfFacesPerCell; ++localFace)
        {
            if(static_cast<int>(localFace) == forbiddenFace)
            {
                continue;
            }
            double const lengthTmp = checkCellFace(cell, localFace, ray, mesh, *length);
            if(lengthTmp > 0.0)
            {
                *length = lengthTmp;
                nextFace = static_cast<int>(localFace);
            }
        }
        return nextFace;
    }

    ALPAKA_FN_ACC double calcCellGain(
        unsigned const cell,
        double const length,
        hase::core::DeviceMeshView const& mesh,
        double const sigmaA,
        double const sigmaE)
    {
        if(mesh.getCellType(cell) == mesh.claddingNumber)
        {
            return alpaka::math::exp(-(mesh.claddingAbsorption) * length);
        }
        return static_cast<double>(
            alpaka::math::exp(mesh.nTot * (mesh.getBetaVolume(cell) * (sigmaE + sigmaA) - sigmaA) * length));
    }

    ALPAKA_FN_ACC double propagateRay(
        hase::core::Ray nextRay,
        unsigned* nextCell,
        hase::core::DeviceMeshView const& mesh,
        double sigmaA,
        double sigmaE)
    {
        double const distanceTotal = nextRay.length;
        double distanceRemaining = nextRay.length;
        double length = 0.0;
        double gain = 1.0;
        int nextForbiddenFace = -1;
        constexpr double numberOfEpsShifts = 64.0;
        constexpr double boundaryNudgeFactor = numberOfEpsShifts * std::numeric_limits<double>::epsilon();

        if(distanceTotal < SMALL)
        {
            return 1.0;
        }

        nextRay = hase::core::normalizeRay(nextRay);
        while(alpaka::math::abs(distanceRemaining) > SMALL)
        {
            assert(*nextCell < mesh.numberOfCells);
            length = distanceRemaining;
            int const nextFace = calcCellRayIntersection(&length, *nextCell, nextRay, nextForbiddenFace, mesh);
            nextRay = calcNextRay(nextRay, length);
            gain *= calcCellGain(*nextCell, length, mesh, sigmaA, sigmaE);
            assert(length >= 0.0);
            distanceRemaining -= length;

            if(nextFace < 0)
            {
                break;
            }

            int const neighbor = mesh.getCellNeighbor(*nextCell, static_cast<unsigned>(nextFace));
            if(neighbor < 0)
            {
                break;
            }
            nextForbiddenFace = mesh.getCellNeighborLocalFace(*nextCell, static_cast<unsigned>(nextFace));
            *nextCell = static_cast<unsigned>(neighbor);

            double const boundaryNudge = boundaryNudgeFactor * distanceTotal;
            static_assert(SMALL > boundaryNudgeFactor);
            if(length < boundaryNudge && distanceRemaining > boundaryNudge)
            {
                nextRay = calcNextRay(nextRay, boundaryNudge);
                distanceRemaining -= boundaryNudge;
            }
        }
        return gain;
    }

    ALPAKA_FN_ACC double propagateRayWithReflection(
        hase::core::Point startPoint,
        hase::core::Point const endPoint,
        unsigned const reflections,
        hase::core::ReflectionPlane,
        unsigned startCell,
        hase::core::DeviceMeshView const& mesh,
        double const sigmaA,
        double const sigmaE)
    {
        (void) reflections;
        assert(reflections == 0u);
        hase::core::Ray ray = hase::core::generateRay(startPoint, endPoint);
        double gain = propagateRay(ray, &startCell, mesh, sigmaA, sigmaE);
        return gain / (ray.length * ray.length);
    }

} // namespace hase::kernels
