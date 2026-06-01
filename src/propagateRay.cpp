/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
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
#include <propagateRay.hpp>
#include <reflection.hpp> /* calcNextReflection */

#include <cassert> /* assert */

/**
 * @brief Checks a level-plane(currentLevel * thickness) for intersection with an ray (zPos, zVec).
 *        If the intersection-length is greater then length.
 *        Than the intersection-length will be returned.
 *        Otherwise 0 will be returned.
 *
 * @return intersection-length if intersection-length <= length
 * @return 0 if intersection-length > length
 *
 **/
ALPAKA_FN_ACC double checkSurface(
    int const currentLevel,
    double const zPos,
    double const zVec,
    double const length,
    double const thickness)
{
    double denominator = zVec;
    if(denominator != 0.0)
    {
        double nominator = currentLevel * thickness - zPos;
        double lengthTmp = nominator / denominator;
        if(lengthTmp <= length && lengthTmp > 0.0)
        {
            return lengthTmp;
        }
    }
    return 0;
}

/**
 * @brief Checks an edges of the given triangle/prism for an intersection
 *        with ray and calculates the intersection-length. If the intersection-length
 *        is greater then length. Than the intersection-length will be
 *        returned. Otherwise 0 will be returned.
 *
 * @return intersection-length if intersection-length <= length
 * @return 0 if intersection-length > length
 **/
ALPAKA_FN_ACC double checkEdge(
    unsigned const triangle,
    int const edge,
    Ray const ray,
    DeviceMeshView const& mesh,
    double const length)
{
    NormalRay normal = mesh.getNormal(triangle, edge);
    double denominator = normal.dir.x * ray.dir.x + normal.dir.y * ray.dir.y;

    if(denominator != 0.0)
    {
        double nominator
            = normal.dir.x * normal.p.x + normal.dir.y * normal.p.y - normal.dir.x * ray.p.x - normal.dir.y * ray.p.y;

        double lengthTmp = nominator / denominator;
        if(lengthTmp <= length && lengthTmp > 0.0)
        {
            return lengthTmp;
        }
    }

    return 0;
}

/**
 * @brief Calculates the intersection-length for the propagated ray and
 *        the current triangle.
 *
 * @return edge number of the intesected edge (-1 for no intersection)
 *
 **/
ALPAKA_FN_ACC int calcTriangleRayIntersection(
    double* length,
    unsigned const triangle,
    Ray const ray,
    unsigned const level,
    int const forbiddenEdge,
    DeviceMeshView const& mesh)
{
    int edge = -1;
// Check 3 edges of triangle
#pragma unroll 3
    for(int edge_i = 0; edge_i < 3; ++edge_i)
    {
        if(edge_i != forbiddenEdge)
        {
            double lengthTmp = checkEdge(triangle, edge_i, ray, mesh, *length);
            if(lengthTmp)
            {
                *length = lengthTmp;
                edge = edge_i;
            }
        }
    }

    // check the upper surface
    if(forbiddenEdge != 3)
    {
        double lengthTmp = checkSurface(level + 1, ray.p.z, ray.dir.z, *length, mesh.thickness);
        if(lengthTmp)
        {
            *length = lengthTmp;
            edge = 3;
        }
    }

    // check the lower surface
    if(forbiddenEdge != 4)
    {
        double lengthTmp = checkSurface(level, ray.p.z, ray.dir.z, *length, mesh.thickness);
        if(lengthTmp)
        {
            *length = lengthTmp;
            edge = 4;
        }
    }
    return edge;
}

/**
 * @brief This is simple vector calculation. The startpoint
 *        of ray will be moved by length.
 *
 * @return ray is the ray with moved startpoint
 *
 **/
ALPAKA_FN_ACC Ray calcNextRay(Ray ray, double const length)
{
    ray.p.x = ray.p.x + length * ray.dir.x;
    ray.p.y = ray.p.y + length * ray.dir.y;
    ray.p.z = ray.p.z + length * ray.dir.z;

    return ray;
}

/**
 * @brief Calculates the gain for the given prism(triangle and level) and
 *        the intersection-length of the ray.
 *
 * @return gain
 *
 **/
ALPAKA_FN_ACC double calcPrismGain(
    unsigned const triangle,
    unsigned const level,
    double const length,
    DeviceMeshView const& mesh,
    double const sigmaA,
    double const sigmaE)
{
    if(mesh.getCellType(triangle) == mesh.claddingNumber)
    {
        return exp(-(mesh.claddingAbsorption) * length);
    }
    else
    {
        return (double) exp(mesh.nTot * (mesh.getBetaVolume(triangle, level) * (sigmaE + sigmaA) - sigmaA) * length);
    }
}

/**
 * @brief Sets the next triangle, next forbiddenEdge
 *        and next level depending on the cutted edge of
 *        the current triangle and the propagated ray.
 *
 **/
ALPAKA_FN_ACC void updateFromEdge(
    unsigned* triangle,
    int* forbiddenEdge,
    unsigned* level,
    DeviceMeshView const& mesh,
    int const edge)
{
    switch(edge)
    {
    case 0:
    case 1:
    case 2:
        // One of three edges
        *forbiddenEdge = mesh.getForbiddenEdge(*triangle, edge);
        *triangle = mesh.getNeighbor(*triangle, edge);
        break;

    case 3:
        // Upper surface
        *forbiddenEdge = 4;
        if(*level != (mesh.numberOfLevels - 2))
            (*level)++;
        break;

    case 4:
        // Lower surface
        *forbiddenEdge = 3;
        if(*level != 0)
            (*level)--;
        break;
    }
}

ALPAKA_FN_ACC double propagateRay(
    Ray nextRay,
    unsigned* nextLevel,
    unsigned* nextTriangle,
    DeviceMeshView const& mesh,
    double sigmaA,
    double sigmaE)
{
    double distanceTotal = nextRay.length;
    double distanceRemaining = nextRay.length;
    double length = 0;
    double gain = 1;
    int nextForbiddenEdge = -1;
    int nextEdge = -1;

    // Length to small, could be same points
    if(distanceTotal < SMALL)
        return 1;

    nextRay = normalizeRay(nextRay);
    while(fabs(distanceRemaining) > SMALL)
    {
        assert(*nextLevel <= mesh.numberOfLevels);
        // Calc gain for triangle intersection
        length = distanceRemaining;
        nextEdge = calcTriangleRayIntersection(&length, *nextTriangle, nextRay, *nextLevel, nextForbiddenEdge, mesh);
        nextRay = calcNextRay(nextRay, length);
        double gainTmp = calcPrismGain(*nextTriangle, *nextLevel, length, mesh, sigmaA, sigmaE);
        gain *= gainTmp;
        assert(length >= 0);

        distanceRemaining -= length;

        // Calc nextTriangle, nextForbiddenEdge and nextLevel
        if(nextEdge != -1)
        {
            updateFromEdge(nextTriangle, &nextForbiddenEdge, nextLevel, mesh, nextEdge);
        }
    }

    return gain;
}

ALPAKA_FN_ACC double propagateRayWithReflection(
    Point startPoint,
    Point const endPoint,
    unsigned const reflections,
    ReflectionPlane reflectionPlane,
    unsigned startLevel,
    unsigned startTriangle,
    DeviceMeshView const& mesh,
    double const sigmaA,
    double const sigmaE)
{
    double distanceTotal = 0;
    double gain = 1.0;

    for(unsigned reflection = 0; reflection < reflections; ++reflection)
    {
        float reflectivity = getReflectivity(reflectionPlane, startTriangle, mesh);
        float totalReflectionAngle = getReflectionAngle(reflectionPlane, mesh);
        Point reflectionPoint = {0, 0, 0};
        double reflectionAngle = 0;

        // Calc reflectionPoint and reflectionAngle
        calcNextReflection(
            startPoint,
            endPoint,
            (reflections - reflection),
            reflectionPlane,
            &reflectionPoint,
            &reflectionAngle,
            mesh);
        Ray reflectionRay = generateRay(startPoint, reflectionPoint);
        distanceTotal += reflectionRay.length;
        gain *= propagateRay(reflectionRay, &startLevel, &startTriangle, mesh, sigmaA, sigmaE);

        assert(reflectionAngle <= 90);
        assert(reflectionAngle >= 0);

        if(reflectionAngle <= totalReflectionAngle)
        {
            gain *= reflectivity;
            if(gain == 0)
            {
                return 0;
            }
        }

        startPoint = reflectionPoint;
        reflectionPlane = reflectionPlane == TOP_REFLECTION ? BOTTOM_REFLECTION : TOP_REFLECTION;
    }

    Ray ray = generateRay(startPoint, endPoint);
    gain *= propagateRay(ray, &startLevel, &startTriangle, mesh, sigmaA, sigmaE);

    distanceTotal += ray.length;

    return gain / (distanceTotal * distanceTotal);
}
