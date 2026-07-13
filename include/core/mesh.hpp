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

/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once

#include <alpaka/alpaka.hpp>
#include <alpaka/core/common.hpp>

#include <alpakaUtils/memory.hpp>
#include <alpakaUtils/utils.hpp>
#include <core/geometry.hpp>

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>
#define REFLECTION_SMALL 1E-3
#define SMALL 1E-5
#define VERY_SMALL 0.0

namespace hase::core
{

    template<class T, class B, class E>
    inline void assertRange(
        [[maybe_unused]] std::vector<T> const& v,
        [[maybe_unused]] B const minElement,
        [[maybe_unused]] E const maxElement,
        [[maybe_unused]] bool const equals)
    {
        if(equals)
        {
            assert(*std::min_element(v.begin(), v.end()) == minElement);
            assert(*std::max_element(v.begin(), v.end()) == maxElement);
        }
        else
        {
            assert(*std::min_element(v.begin(), v.end()) >= minElement);
            assert(*std::max_element(v.begin(), v.end()) <= maxElement);
        }
    }

    template<class T, class B>
    inline void assertMin(
        [[maybe_unused]] std::vector<T> const& v,
        [[maybe_unused]] B const minElement,
        [[maybe_unused]] bool const equals)
    {
        if(equals)
        {
            assert(*std::min_element(v.begin(), v.end()) == minElement);
        }
        else
        {
            assert(*std::min_element(v.begin(), v.end()) >= minElement);
        }
    }

    inline double distance2D(TwoDimPoint const p1, TwoDimPoint const p2)
    {
        return std::abs(std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y)));
    }

    inline double getMaxDistance(std::vector<TwoDimPoint> const& points)
    {
        double maxDistance = -1.0;

        for(unsigned p1 = 0; p1 < points.size(); ++p1)
        {
            for(unsigned p2 = p1; p2 < points.size(); ++p2)
            {
                maxDistance = std::max(maxDistance, distance2D(points[p1], points[p2]));
            }
        }

        return maxDistance;
    }

    template<typename Mesh>
    [[nodiscard]] constexpr float getReflectivity(ReflectionPlane reflectionPlane, unsigned triangle, Mesh const& mesh)
    {
        switch(reflectionPlane)
        {
        case BOTTOM_REFLECTION:
            return mesh.reflectivities[triangle];
        case TOP_REFLECTION:
            return mesh.reflectivities[triangle + mesh.numberOfTriangles];
        }
        return 0.0f;
    }

    template<typename Mesh>
    [[nodiscard]] constexpr float getReflectionAngle(ReflectionPlane reflectionPlane, Mesh const& mesh)
    {
        switch(reflectionPlane)
        {
        case BOTTOM_REFLECTION:
            return mesh.totalReflectionAngles[0];
        case TOP_REFLECTION:
            return mesh.totalReflectionAngles[1];
        }
        return 0.0f;
    }

    inline double calculateMaxDiameter(double const* points, unsigned const offset)
    {
        TwoDimPoint minX = {DBL_MAX, 0};
        TwoDimPoint minY = {0, DBL_MAX};
        TwoDimPoint maxX = {DBL_MIN, 0};
        TwoDimPoint maxY = {0, DBL_MIN};

        for(unsigned p = 0; p < offset; ++p)
        {
            TwoDimPoint np = {points[p], points[p + offset]};
            minX = (points[p] < minX.x) ? np : minX;
            maxX = (points[p] > maxX.x) ? np : maxX;
        }

        for(unsigned p = offset; p < 2 * offset; ++p)
        {
            TwoDimPoint np = {points[p - offset], points[p]};
            minY = points[p] < minY.y ? np : minY;
            maxY = points[p] > maxY.y ? np : maxY;
        }

        std::vector<TwoDimPoint> extrema;
        extrema.push_back(minX);
        extrema.push_back(minY);
        extrema.push_back(maxX);
        extrema.push_back(maxY);

        return getMaxDistance(extrema);
    }

    struct DeviceMeshView
    {
        std::span<double const> points;
        std::span<double const> betaVolume;
        std::span<double const> normalVec;
        std::span<double const> centers;
        std::span<float const> triangleSurfaces;
        std::span<int const> forbiddenEdge;
        std::span<double const> betaCells;
        std::span<unsigned const> claddingCellTypes;
        std::span<float const> refractiveIndices;
        std::span<float const> reflectivities;
        std::span<float const> totalReflectionAngles;
        std::span<unsigned const> trianglePointIndices;
        std::span<int const> triangleNeighbors;
        std::span<unsigned const> triangleNormalPoint;

        double claddingAbsorption;
        float surfaceTotal;
        float thickness;
        float nTot;
        float crystalTFluo;
        unsigned numberOfTriangles;
        unsigned numberOfLevels;
        unsigned numberOfPrisms;
        unsigned numberOfPoints;
        unsigned numberOfSamples;
        unsigned claddingNumber;

        [[nodiscard]] ALPAKA_FN_ACC int getNeighbor(unsigned triangle, int edge) const
        {
            return triangleNeighbors[triangle + edge * numberOfTriangles];
        }

        ALPAKA_FN_ACC Point genRndPoint(
            Point& origin,
            unsigned triangle,
            unsigned level,
            alpaka::rand::engine::Philox4x32x10& rndEngine) const
        {
            Point startPoint = {0, 0, 0};

            double u = alpaka::rand::distribution::UniformReal<double>{}(rndEngine);
            double v = alpaka::rand::distribution::UniformReal<double>{}(rndEngine);

            if((u + v) > 1.0)
            {
                u = 1.0 - u;
                v = 1.0 - v;
            }

            double w = 1.0 - u - v;
            unsigned t1 = trianglePointIndices[triangle];
            unsigned t2 = trianglePointIndices[triangle + numberOfTriangles];
            unsigned t3 = trianglePointIndices[triangle + 2 * numberOfTriangles];

            startPoint.z = level * thickness + alpaka::rand::distribution::UniformReal<double>{}(rndEngine) *thickness;
            startPoint.x = (points[t1] * u) + (points[t2] * v) + (points[t3] * w);
            startPoint.y = (points[t1 + numberOfPoints] * u) + (points[t2 + numberOfPoints] * v)
                           + (points[t3 + numberOfPoints] * w);
            if((origin - startPoint).euclidLength() < SMALL)
            {
                // regenerate ray point if two close to origin point
                // -> fixes numerical instability 1/inf issue yet causes small calculation bias
                return genRndPoint(origin, triangle, level, rndEngine);
            }
            return startPoint;
        }

        [[nodiscard]] ALPAKA_FN_ACC double getBetaVolume(unsigned triangle, unsigned level) const
        {
            return betaVolume[triangle + level * numberOfTriangles];
        }

        [[nodiscard]] ALPAKA_FN_ACC double getBetaVolume(unsigned prism) const
        {
            return betaVolume[prism];
        }

        ALPAKA_FN_ACC NormalRay getNormal(unsigned triangle, int edge) const
        {
            NormalRay ray = {{0, 0}, {0, 0}};
            int offset = edge * numberOfTriangles + triangle;

            ray.p.x = points[triangleNormalPoint[offset]];
            ray.p.y = points[triangleNormalPoint[offset] + numberOfPoints];
            ray.dir.x = normalVec[offset];
            ray.dir.y = normalVec[offset + 3 * numberOfTriangles];

            return ray;
        }

        ALPAKA_FN_ACC Point getSamplePoint(unsigned sample_i) const
        {
            Point p = {0, 0, 0};
            unsigned level = sample_i / numberOfPoints;
            unsigned pos = sample_i - (numberOfPoints * level);

            p.z = level * thickness;
            p.x = points[pos];
            p.y = points[pos + numberOfPoints];

            return p;
        }

        ALPAKA_FN_ACC alpaka::concepts::Simd auto getSimdCenterPoint(
            alpaka::concepts::Simd auto triangles,
            alpaka::concepts::Simd auto level) const
        {
            constexpr uint32_t width = ALPAKA_TYPEOF(triangles)::width();
            alpaka::concepts::Simd auto points = alpaka::Simd<Point, width>::fill(Point{0, 0, 0});
            for(uint32_t laneIdx = 0; laneIdx < width; ++laneIdx)
            {
                points[laneIdx].x = centers[triangles[laneIdx]];
                points[laneIdx].y = centers[triangles[laneIdx] + numberOfTriangles];
                points[laneIdx].z = (level[laneIdx] + 0.5f) * thickness;
            }
            return points;
        }

        ALPAKA_FN_ACC Point getCenterPoint(unsigned triangle, unsigned level) const
        {
            Point p = {0, 0, (level + 0.5f) * thickness};
            p.x = centers[triangle];
            p.y = centers[triangle + numberOfTriangles];
            return p;
        }

        [[nodiscard]] ALPAKA_FN_ACC int getForbiddenEdge(unsigned triangle, int edge) const
        {
            return forbiddenEdge[edge * numberOfTriangles + triangle];
        }

        [[nodiscard]] ALPAKA_FN_ACC unsigned getCellType(unsigned triangle) const
        {
            return claddingCellTypes[triangle];
        }

        ALPAKA_FN_ACC void test() const
        {
            printf("Constants:\n");
            printf("claddingAbsorption: %f\n", claddingAbsorption);
            printf("surfaceTotal: %f\n", surfaceTotal);
            printf("thickness: %f\n", thickness);
            printf("nTot: %f\n", nTot);
            printf("crystalTFluo: %f\n", crystalTFluo);
            printf("numberOfTriangles: %u\n", numberOfTriangles);
            printf("numberOfLevels: %u\n", numberOfLevels);
            printf("numberOfPrisms: %u\n", numberOfPrisms);
            printf("numberOfPoints: %u\n", numberOfPoints);
            printf("numberOfSamples: %u\n", numberOfSamples);
            printf("claddingNumber: %u\n", claddingNumber);
        }
    };
    template<alpaka::onHost::concepts::Device T_Device>
    class DeviceMeshContainer;

    /**
     * @brief fills a device mesh with the correct datastructures
     *
     * See parseMultiGPU for details on the parameters
     */
    template<alpaka::onHost::concepts::Device T_Device>
    DeviceMeshContainer<T_Device> createMesh(
        T_Device& device,
        std::vector<unsigned> const& triangleIndices,
        unsigned const numberOfTriangles,
        unsigned const numberOfLevels,
        unsigned const numberOfPoints,
        float const thicknessOfPrism,
        std::vector<double>& pointsVector,
        std::vector<double>& xOfTriangleCenter,
        std::vector<double>& yOfTriangleCenter,
        std::vector<unsigned>& positionsOfNormalVectors,
        std::vector<double>& xOfNormals,
        std::vector<double>& yOfNormals,
        std::vector<int>& forbiddenVector,
        std::vector<int>& neighborsVector,
        std::vector<float>& surfacesVector,
        std::vector<double>& betaValuesVector,
        std::vector<double>& betaCells,
        std::vector<unsigned>& cellTypes,
        std::vector<float>& refractiveIndices,
        std::vector<float>& reflectivities,
        std::vector<float>& totalReflectionAngles,
        float const nTot,
        float const crystalFluorescence,
        unsigned const cladNumber,
        double const cladAbsorption

    )
    {
        // GPU variables
        double totalSurface = 0.;
        totalSurface += std::reduce(surfacesVector.begin(), surfacesVector.begin() + numberOfTriangles, 0.0);

        // Vector Preprocessing
        std::vector<double> hostNormalVec(xOfNormals.begin(), xOfNormals.end());
        hostNormalVec.insert(hostNormalVec.end(), yOfNormals.begin(), yOfNormals.end());
        std::vector<double> hostCenters(xOfTriangleCenter.begin(), xOfTriangleCenter.end());
        hostCenters.insert(hostCenters.end(), yOfTriangleCenter.begin(), yOfTriangleCenter.end());


        return {
            device,
            cladAbsorption,
            static_cast<float>(totalSurface), // this cast might needlessly cause precision loss
            thicknessOfPrism,
            nTot,
            crystalFluorescence,
            numberOfTriangles,
            numberOfLevels,
            numberOfTriangles * (numberOfLevels - 1),
            numberOfPoints,
            numberOfPoints * numberOfLevels,
            cladNumber,
            pointsVector,
            hostNormalVec,
            betaValuesVector,
            hostCenters,
            surfacesVector,
            forbiddenVector,
            betaCells,
            cellTypes,
            refractiveIndices,
            reflectivities,
            totalReflectionAngles,
            triangleIndices,
            neighborsVector,
            positionsOfNormalVectors};
    }

    /**
     * @brief Contains the structure of the crystal
     *
     * All the fixed values of how the crystal is meshed
     *
     * points The coordinates of the triangle vertices
     *        All x coordinates followed by all of the y coordinates of the triangle vertices
     *        structure: [x_1, x_2, ... x_n, y_1, y_2, ... y_n] (n == numberOfPoints)
     *
     * betaVolume beta values for all prisms ordered accordingly to the prismIDs:
     *            prismID = triangleID + layer * numberOfTriangles;
     *            therefore, all betaVolume for a layer are grouped together
     *
     * normalVec the normal vectors for each triangle edge
     *           first half (size: 3*numberOfTriangles -> one for each side) contains
     *           the x components of each vector, second half contains the y components.
     *           the each half is ordered as follows:
     *           [ triangle1edge0, triangle2edge0, ... triangleNedge0, triangle1edge1, triangle2edge1, ... ]
     *           i.e. all first edges of each triangle, followed by all second edges of each triangle and so on.
     *
     * centers the coordinates of the center points for each triangle
     *         All x coordinates followed by all y coordinates of the triangle vertices
     *         similar to "points"
     *
     * triangleSurfaces the sizes of the surfaces of each triangle, ordered by the triangleID
     *
     * forbiddenEdge  describes the relation of edge indices of adjacent triangles
     *           -1 means, there is no adjacent triangle to that edge
     *           0,1,2 describes the index of the edge as seen from the ADJACENT triangle
     *
     *           order of data is similar to normalVec:
     *           [ triangle1edge0, triangle2edge0, ... triangleNedge0, triangle1edge1, triangle2edge1, ... ]
     *           i.e. all first edges of each triangle, followed by all second edges of each triangle and so on.
     *
     * trianglesPointIndices contains the indices to access the "points" datastructure
     *           (each triangle has 3 points as vertices). Each entry is an
     *           index from 0 to numberOfPoints, corresponding to the positions
     *           of a vertex in "points".
     *           structure is similar to "forbiddenEdge":
     *           [ triangle1A, triangle2A, ... triangleNA, triangle1B, triangle2B, ... triangleNB, triangle1C, ... ]
     *           i.e. for triangles with vertices A,B,C there are all the indices
     *           of the A-vertices, followed by all the B and C vertices.
     *
     * triangleNeighbors describes the relation of triangle indices to each other.
     *           Each entry corresponds to a triangleID (see "trianglePointIndices") which
     *           is adjacent to the current triangle.
     *           structure is similar to "forbiddenEdge":
     *           [ triangle1edge0, triangle2edge0, ... triangleNedge0, triangle1edge1, triangle2edge1, ... ]
     *
     * triangleNormalPoint contains indices to the x and y components of the positions where the
     *             normalVectors start (see normalVec). For each Triangle 3 points (3 sides)
     *             are stored in this list.
     *             Indices point to locations in "points" (i.e. normal vectors start at
     *             triangle vertices!)
     *             structure is VERY similar to trianglePointIndices:
     *             [ triangle1p0, triangle2p0, ... triangleNp0, triangle1p1, triangle2p1, ... ]
     *
     * refractiveIndices [0]->bottomInside, [1]->bottomOutside, [2]->topInside, [3]->topOutside
     *
     * reflectivities Contains the reflectivities for upper and lower surface of gain medium
     *                Structure is based on 2 layers of triangles:
     *                [refl_tri1_bott, refl_tri2_bott, ...,refl_triN_bott, refl_tri1_top, refl_tri2_top, ...,
     * refl_triN_top]
     *
     * totalReflectionAngles [0]-> bottomTotalReflectionAngle, [1]-> topTotalReflectionAngle
     */
    template<alpaka::onHost::concepts::Device T_Device>
    class DeviceMeshContainer
    {
        using T_Extent = alpaka::Vec<std::size_t, 1>;
        using T_Queue = ALPAKA_TYPEOF(std::declval<T_Device>().makeQueue(alpaka::queueKind::blocking));

    public:
        DeviceMeshContainer(
            T_Device device,
            double claddingAbsorption,
            float surfaceTotal,
            float thickness,
            float nTot,
            float crystalTFluo,
            unsigned numberOfTriangles,
            unsigned numberOfLevels,
            unsigned numberOfPrisms,
            unsigned numberOfPoints,
            unsigned numberOfSamples,
            unsigned claddingNumber,
            std::vector<double> points,
            std::vector<double> normalVec,
            std::vector<double> betaVolume,
            std::vector<double> centers,
            std::vector<float> triangleSurfaces,
            std::vector<int> forbiddenEdge,
            std::vector<double> betaCells,
            std::vector<unsigned> claddingCellTypes,
            std::vector<float> refractiveIndices,
            std::vector<float> reflectivities,
            std::vector<float> totalReflectionAngles,
            std::vector<unsigned> trianglePointIndices,
            std::vector<int> triangleNeighbors,
            std::vector<unsigned> triangleNormalPoint)
            : m_device(device)
            , m_queue(device.makeQueue(alpaka::queueKind::blocking))
            , points(hase::alpakaUtils::toDevice(m_queue, points))
            , betaVolume(hase::alpakaUtils::toDevice(m_queue, betaVolume))
            , normalVec(hase::alpakaUtils::toDevice(m_queue, normalVec))
            , centers(hase::alpakaUtils::toDevice(m_queue, centers))
            , triangleSurfaces(hase::alpakaUtils::toDevice(m_queue, triangleSurfaces))
            , forbiddenEdge(hase::alpakaUtils::toDevice(m_queue, forbiddenEdge))
            , betaCells(hase::alpakaUtils::toDevice(m_queue, betaCells))
            , claddingCellTypes(hase::alpakaUtils::toDevice(m_queue, claddingCellTypes))
            , refractiveIndices(hase::alpakaUtils::toDevice(m_queue, refractiveIndices))
            , reflectivities(hase::alpakaUtils::toDevice(m_queue, reflectivities))
            , totalReflectionAngles(hase::alpakaUtils::toDevice(m_queue, totalReflectionAngles))
            , trianglePointIndices(hase::alpakaUtils::toDevice(m_queue, trianglePointIndices))
            , triangleNeighbors(hase::alpakaUtils::toDevice(m_queue, triangleNeighbors))
            , triangleNormalPoint(hase::alpakaUtils::toDevice(m_queue, triangleNormalPoint))
            , claddingAbsorption(claddingAbsorption)
            , surfaceTotal(surfaceTotal)
            , thickness(thickness)
            , nTot(nTot)
            , crystalTFluo(crystalTFluo)
            , numberOfTriangles(numberOfTriangles)
            , numberOfLevels(numberOfLevels)
            , numberOfPrisms(numberOfPrisms)
            , numberOfPoints(numberOfPoints)
            , numberOfSamples(numberOfSamples)
            , claddingNumber(claddingNumber)
        {
        }

        ~DeviceMeshContainer() = default;

        [[nodiscard]] auto toView() const -> DeviceMeshView
        {
            return {
                std::span<double const>(points.data(), points.getMdSpan().getExtents().x()),
                std::span<double const>(betaVolume.data(), betaVolume.getMdSpan().getExtents().x()),
                std::span<double const>(normalVec.data(), normalVec.getMdSpan().getExtents().x()),
                std::span<double const>(centers.data(), centers.getMdSpan().getExtents().x()),
                std::span<float const>(triangleSurfaces.data(), triangleSurfaces.getMdSpan().getExtents().x()),
                std::span<int const>(forbiddenEdge.data(), forbiddenEdge.getMdSpan().getExtents().x()),
                std::span<double const>(betaCells.data(), betaCells.getMdSpan().getExtents().x()),
                std::span<unsigned const>(claddingCellTypes.data(), claddingCellTypes.getMdSpan().getExtents().x()),
                std::span<float const>(refractiveIndices.data(), refractiveIndices.getMdSpan().getExtents().x()),
                std::span<float const>(reflectivities.data(), reflectivities.getMdSpan().getExtents().x()),
                std::span<float const>(
                    totalReflectionAngles.data(),
                    totalReflectionAngles.getMdSpan().getExtents().x()),
                std::span<unsigned const>(
                    trianglePointIndices.data(),
                    trianglePointIndices.getMdSpan().getExtents().x()),
                std::span<int const>(triangleNeighbors.data(), triangleNeighbors.getMdSpan().getExtents().x()),
                std::span<unsigned const>(
                    triangleNormalPoint.data(),
                    triangleNormalPoint.getMdSpan().getExtents().x()),
                claddingAbsorption,
                surfaceTotal,
                thickness,
                nTot,
                crystalTFluo,
                numberOfTriangles,
                numberOfLevels,
                numberOfPrisms,
                numberOfPoints,
                numberOfSamples,
                claddingNumber};
        }

        T_Device m_device;

    private:
        T_Queue m_queue;

    public:
        template<typename T_Data>
        using T_Buffer = std::remove_cvref_t<decltype(hase::alpakaUtils::toDevice(
            std::declval<T_Queue const&>(),
            std::declval<std::vector<T_Data> const&>()))>;
        T_Buffer<double> points;
        T_Buffer<double> betaVolume;
        T_Buffer<double> normalVec;
        T_Buffer<double> centers;
        T_Buffer<float> triangleSurfaces;
        T_Buffer<int> forbiddenEdge;
        T_Buffer<double> betaCells;
        T_Buffer<unsigned> claddingCellTypes;
        T_Buffer<float> refractiveIndices;
        T_Buffer<float> reflectivities;
        T_Buffer<float> totalReflectionAngles;
        T_Buffer<unsigned> trianglePointIndices;
        T_Buffer<int> triangleNeighbors;
        T_Buffer<unsigned> triangleNormalPoint;

        double claddingAbsorption;
        float surfaceTotal;
        float thickness;
        float nTot;
        float crystalTFluo;
        unsigned numberOfTriangles;
        unsigned numberOfLevels;
        unsigned numberOfPrisms;
        unsigned numberOfPoints;
        unsigned numberOfSamples;
        unsigned claddingNumber;
    };

    class HostMesh
    {
    public:
        std::vector<unsigned> trianglePointIndices;
        unsigned numberOfTriangles;
        unsigned numberOfLevels;
        unsigned numberOfPoints;
        float thickness;
        std::vector<double> points;
        std::vector<double> triangleCenterX;
        std::vector<double> triangleCenterY;
        std::vector<unsigned> triangleNormalPoint;
        std::vector<double> triangleNormalsX;
        std::vector<double> triangleNormalsY;
        std::vector<int> forbiddenEdge;
        std::vector<int> triangleNeighbors;
        std::vector<float> triangleSurfaces;
        std::vector<double> betaVolume;
        std::vector<double> betaCells;
        std::vector<unsigned> claddingCellTypes;
        std::vector<float> refractiveIndices;
        std::vector<float> reflectivities;
        float nTot;
        float crystalTFluo;
        unsigned claddingNumber;
        double claddingAbsorption;
        std::vector<float> totalReflectionAngles;

        HostMesh() = default;

        HostMesh(
            std::vector<unsigned> trianglePointIndices,
            unsigned numberOfTriangles,
            unsigned numberOfLevels,
            unsigned numberOfPoints,
            float thickness,
            std::vector<double> points,
            std::vector<double> triangleCenterX,
            std::vector<double> triangleCenterY,
            std::vector<unsigned> triangleNormalPoint,
            std::vector<double> triangleNormalsX,
            std::vector<double> triangleNormalsY,
            std::vector<int> forbiddenEdge,
            std::vector<int> triangleNeighbors,
            std::vector<float> triangleSurfaces,
            std::vector<double> betaVolume,
            std::vector<double> betaCells,
            std::vector<unsigned> claddingCellTypes,
            std::vector<float> refractiveIndices,
            std::vector<float> reflectivities,
            float nTot,
            float crystalTFluo,
            unsigned claddingNumber,
            double claddingAbsorption)
            : trianglePointIndices(std::move(trianglePointIndices))
            , numberOfTriangles(numberOfTriangles)
            , numberOfLevels(numberOfLevels)
            , numberOfPoints(numberOfPoints)
            , thickness(thickness)
            , points(std::move(points))
            , triangleCenterX(std::move(triangleCenterX))
            , triangleCenterY(std::move(triangleCenterY))
            , triangleNormalPoint(std::move(triangleNormalPoint))
            , triangleNormalsX(std::move(triangleNormalsX))
            , triangleNormalsY(std::move(triangleNormalsY))
            , forbiddenEdge(std::move(forbiddenEdge))
            , triangleNeighbors(std::move(triangleNeighbors))
            , triangleSurfaces(std::move(triangleSurfaces))
            , betaVolume(std::move(betaVolume))
            , betaCells(std::move(betaCells))
            , claddingCellTypes(std::move(claddingCellTypes))
            , refractiveIndices(std::move(refractiveIndices))
            , reflectivities(std::move(reflectivities))
            , nTot(nTot)
            , crystalTFluo(crystalTFluo)
            , claddingNumber(claddingNumber)
            , claddingAbsorption(claddingAbsorption)
        {
        }

        void calcTotalReflectionAngles()
        {
            std::vector<float> totalReflectionAngles_(refractiveIndices.size() / 2, 0);
            for(unsigned i = 0; i < refractiveIndices.size(); i += 2)
            {
                totalReflectionAngles_.at(i / 2)
                    = (180. / M_PI * alpaka::math::asin(refractiveIndices.at(i + 1) / refractiveIndices.at(i)));
            }
            totalReflectionAngles = std::move(totalReflectionAngles_);
        }

        template<typename T_Device>
        [[nodiscard]] DeviceMeshContainer<T_Device> toDevice(T_Device& device)
        {
            return createMesh(
                device,
                trianglePointIndices,
                numberOfTriangles,
                numberOfLevels,
                numberOfPoints,
                thickness,
                points,
                triangleCenterX,
                triangleCenterY,
                triangleNormalPoint,
                triangleNormalsX,
                triangleNormalsY,
                forbiddenEdge,
                triangleNeighbors,
                triangleSurfaces,
                betaVolume,
                betaCells,
                claddingCellTypes,
                refractiveIndices,
                reflectivities,
                totalReflectionAngles,
                nTot,
                crystalTFluo,
                claddingNumber,
                claddingAbsorption);
        }

        [[nodiscard]] unsigned getMaxReflections(ReflectionPlane reflectionPlane) const
        {
            double d = calculateMaxDiameter(points.data(), numberOfPoints);
            float alpha = getReflectionAngle(reflectionPlane, *this) * static_cast<float>(M_PI) / 180.0f;
            double h = (numberOfLevels - 1u) * thickness;
            double z = d / std::tan(alpha);
            return static_cast<unsigned>(std::ceil(z / h));
        }

        [[nodiscard]] unsigned getMaxReflections() const
        {
            unsigned top = getMaxReflections(TOP_REFLECTION);
            unsigned bottom = getMaxReflections(BOTTOM_REFLECTION);
            return std::max(top, bottom);
        }
    };

} // namespace hase::core
