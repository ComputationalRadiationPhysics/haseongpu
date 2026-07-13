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


#include <core/geometry.hpp>
#include <core/mesh.hpp>
#include <kernels/reflection.hpp>

#include <cassert>
#include <cmath> /* M_PI */

namespace hase::kernels
{

    ALPAKA_FN_ACC double calcIntersectionAngle(hase::core::Ray const ray, double* reflectionAngle)
    {
        // Calc intesection angle with z-plane
        double nominator = alpaka::math::abs(ray.dir.z);
        double denominator
            = alpaka::math::sqrt((ray.dir.x * ray.dir.x) + (ray.dir.y * ray.dir.y) + (ray.dir.z * ray.dir.z));
        if(denominator != 0.0)
        {
            double radian = alpaka::math::acos(nominator / denominator);
            *reflectionAngle = ((180. / M_PI) * radian);
            return 0;
        }
        return 1;
    }

    ALPAKA_FN_ACC int calcPlaneIntersectionPoint(
        hase::core::Ray const reflectionRay,
        hase::core::ReflectionPlane const reflectionPlane,
        hase::core::DeviceMeshView const& mesh,
        hase::core::Point* intersectionPoint)
    {
        // Assume that mesh is on x/y axis and parallel to x/y axis
        double planeZ = 0.0;
        if(reflectionPlane == hase::core::TOP_REFLECTION)
        {
            // Reflection on TOP plane
            planeZ = mesh.thickness * (mesh.numberOfLevels - 1u);
        }
        double denominator = reflectionRay.dir.z;
        if(denominator != 0.0)
        {
            double nominator = planeZ - reflectionRay.p.z;
            double length = nominator / denominator;
            if(length > 0)
            {
                intersectionPoint->x = reflectionRay.p.x + length * reflectionRay.dir.x;
                intersectionPoint->y = reflectionRay.p.y + length * reflectionRay.dir.y;
                intersectionPoint->z = planeZ;
                return 0;
            }
        }
        return 1;
    }

    /*
     * TOP_REFLECTION = 1
     * BOTTOM_REFLECTION = -1
     * defined in reflection.hpp
     */
    ALPAKA_FN_ACC hase::core::Ray generateReflectionRay(
        hase::core::Point const startPoint,
        hase::core::Point endPoint,
        int const reflectionsLeft,
        hase::core::ReflectionPlane const reflectionPlane,
        hase::core::DeviceMeshView const& mesh)
    {
        float mirrorPlaneZ = 0;
        double const reflectionSign = reflectionPlane;
        if(reflectionsLeft % 2 == 0)
        {
            // Even reflectionCount is postponement

            endPoint.z = endPoint.z + reflectionSign * (reflectionsLeft * mesh.thickness * (mesh.numberOfLevels - 1u));
        }
        else
        {
            // Odd reflectionsCount is reflection

            if(reflectionPlane == hase::core::TOP_REFLECTION)
            {
                mirrorPlaneZ
                    = alpaka::math::ceil(reflectionsLeft / (double) 2) * mesh.thickness * (mesh.numberOfLevels - 1u);
            }
            else
            {
                mirrorPlaneZ = alpaka::math::floor(reflectionsLeft / (double) 2) * mesh.thickness
                               * (mesh.numberOfLevels - 1u) * reflectionSign;
            }

            endPoint.z = reflectionSign * alpaka::math::abs((mirrorPlaneZ + mirrorPlaneZ - endPoint.z));
        }
        return hase::core::generateRay(startPoint, endPoint);
    }

    ALPAKA_FN_ACC int calcNextReflection(
        hase::core::Point const startPoint,
        hase::core::Point const endPoint,
        unsigned const reflectionsLeft,
        hase::core::ReflectionPlane const reflectionPlane,
        hase::core::Point* reflectionPoint,
        double* reflectionAngle,
        hase::core::DeviceMeshView const& mesh)
    {
        hase::core::Ray reflectionRay
            = generateReflectionRay(startPoint, endPoint, reflectionsLeft, reflectionPlane, mesh);
        if(calcPlaneIntersectionPoint(reflectionRay, reflectionPlane, mesh, reflectionPoint))
            return 1;
        if(calcIntersectionAngle(reflectionRay, reflectionAngle))
            return 1;

        return 0;
    }

} // namespace hase::kernels
