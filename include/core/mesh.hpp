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
#include <limits>
#include <utility>
#include <vector>
#define REFLECTION_SMALL 1E-3
#define SMALL 1E-5
#define VERY_SMALL 0.0

namespace hase::core
{
    constexpr unsigned tet4VertexCount = 4u;
    constexpr unsigned tet4FaceCount = 4u;
    constexpr unsigned tet4FaceWidth = 3u;
    constexpr unsigned vtkTetraCellType = 10u;

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
        std::span<double const> betaCells;
        std::span<unsigned const> claddingCellTypes;
        std::span<float const> refractiveIndices;
        std::span<float const> reflectivities;
        std::span<float const> totalReflectionAngles;
        std::span<unsigned const> cellPointIndices;
        std::span<unsigned const> cellTypes;
        std::span<int const> cellFaces;
        std::span<int const> cellNeighborCells;
        std::span<int const> cellNeighborLocalFaces;
        std::span<int const> cellFaceBoundaries;
        std::span<float const> cellVolumes;
        std::span<double const> cellCenters;
        std::span<double const> samplePoints;

        double claddingAbsorption;
        float nTot;
        float crystalTFluo;
        unsigned claddingNumber;
        unsigned numberOfCells;
        unsigned numberOfPrisms;
        unsigned numberOfPoints;
        unsigned numberOfSamples;
        unsigned numberOfFacesPerCell;
        unsigned numberOfCellVertices;

        [[nodiscard]] ALPAKA_FN_ACC Point getPoint(unsigned pointIndex) const
        {
            return Point{
                points[pointIndex],
                points[pointIndex + numberOfPoints],
                points[pointIndex + 2u * numberOfPoints]};
        }

        [[nodiscard]] ALPAKA_FN_ACC Point getCellPoint(unsigned cell, unsigned localVertex) const
        {
            return getPoint(cellPointIndices[cell * numberOfCellVertices + localVertex]);
        }

        [[nodiscard]] ALPAKA_FN_ACC int getCellFacePoint(unsigned cell, unsigned localFace, unsigned localVertex) const
        {
            return cellFaces[(cell * numberOfFacesPerCell + localFace) * tet4FaceWidth + localVertex];
        }

        [[nodiscard]] ALPAKA_FN_ACC int getCellNeighbor(unsigned cell, unsigned localFace) const
        {
            return cellNeighborCells[cell * numberOfFacesPerCell + localFace];
        }

        [[nodiscard]] ALPAKA_FN_ACC int getCellNeighborLocalFace(unsigned cell, unsigned localFace) const
        {
            return cellNeighborLocalFaces[cell * numberOfFacesPerCell + localFace];
        }

        [[nodiscard]] ALPAKA_FN_ACC double getBetaVolume(unsigned cell) const
        {
            return betaVolume[cell];
        }

        [[nodiscard]] ALPAKA_FN_ACC unsigned getCellType(unsigned cell) const
        {
            return claddingCellTypes[cell];
        }

        [[nodiscard]] ALPAKA_FN_ACC double tetraVolume(
            Point const a,
            Point const b,
            Point const c,
            Point const d) const
        {
            return alpaka::math::abs(dot(cross(b - a, c - a), d - a)) / 6.0;
        }

        ALPAKA_FN_ACC Point genRndPointInTetra(
            Point const a,
            Point const b,
            Point const c,
            Point const d,
            alpaka::rand::engine::Philox4x32x10& rndEngine) const
        {
            double r0 = alpaka::rand::distribution::UniformReal<double>{}(rndEngine);
            double r1 = alpaka::rand::distribution::UniformReal<double>{}(rndEngine);
            double r2 = alpaka::rand::distribution::UniformReal<double>{}(rndEngine);
            double r3 = alpaka::rand::distribution::UniformReal<double>{}(rndEngine);
            r0 = -alpaka::math::log(alpaka::math::max(r0, std::numeric_limits<double>::min()));
            r1 = -alpaka::math::log(alpaka::math::max(r1, std::numeric_limits<double>::min()));
            r2 = -alpaka::math::log(alpaka::math::max(r2, std::numeric_limits<double>::min()));
            r3 = -alpaka::math::log(alpaka::math::max(r3, std::numeric_limits<double>::min()));
            double const invSum = 1.0 / (r0 + r1 + r2 + r3);
            return a * (r0 * invSum) + b * (r1 * invSum) + c * (r2 * invSum) + d * (r3 * invSum);
        }

        ALPAKA_FN_ACC Point genRndPointInCell(
            Point& origin,
            unsigned cell,
            alpaka::rand::engine::Philox4x32x10& rndEngine) const
        {
            Point const p0 = getCellPoint(cell, 0u);
            Point const p1 = getCellPoint(cell, 1u);
            Point const p2 = getCellPoint(cell, 2u);
            Point const p3 = getCellPoint(cell, 3u);

            Point startPoint = genRndPointInTetra(p0, p1, p2, p3, rndEngine);
            if((origin - startPoint).euclidLength() < SMALL)
            {
                return genRndPointInCell(origin, cell, rndEngine);
            }
            return startPoint;
        }

        [[nodiscard]] ALPAKA_FN_ACC Point getSamplePoint(unsigned sampleIndex) const
        {
            return Point{
                samplePoints[sampleIndex],
                samplePoints[sampleIndex + numberOfSamples],
                samplePoints[sampleIndex + 2u * numberOfSamples]};
        }

        [[nodiscard]] ALPAKA_FN_ACC Point getCellCenterPoint(unsigned cell) const
        {
            return Point{
                cellCenters[cell],
                cellCenters[cell + numberOfCells],
                cellCenters[cell + 2u * numberOfCells]};
        }

        [[nodiscard]] ALPAKA_FN_ACC double getCellVolume(unsigned cell) const
        {
            return static_cast<double>(cellVolumes[cell]);
        }
    };

    template<alpaka::onHost::concepts::Device T_Device>
    class DeviceMeshContainer
    {
        using T_Queue = ALPAKA_TYPEOF(std::declval<T_Device>().makeQueue(alpaka::queueKind::blocking));

    public:
        DeviceMeshContainer(
            T_Device device,
            double claddingAbsorption,
            float nTot,
            float crystalTFluo,
            unsigned claddingNumber,
            unsigned numberOfCells,
            unsigned numberOfPoints,
            unsigned numberOfSamples,
            unsigned numberOfFacesPerCell,
            unsigned numberOfCellVertices,
            std::vector<double> points,
            std::vector<double> betaVolume,
            std::vector<double> betaCells,
            std::vector<unsigned> claddingCellTypes,
            std::vector<float> refractiveIndices,
            std::vector<float> reflectivities,
            std::vector<float> totalReflectionAngles,
            std::vector<unsigned> cellPointIndices,
            std::vector<unsigned> cellTypes,
            std::vector<int> cellFaces,
            std::vector<int> cellNeighborCells,
            std::vector<int> cellNeighborLocalFaces,
            std::vector<int> cellFaceBoundaries,
            std::vector<float> cellVolumes,
            std::vector<double> cellCenters,
            std::vector<double> samplePoints)
            : m_device(device)
            , m_queue(device.makeQueue(alpaka::queueKind::blocking))
            , points(hase::alpakaUtils::toDevice(m_queue, points))
            , betaVolume(hase::alpakaUtils::toDevice(m_queue, betaVolume))
            , betaCells(hase::alpakaUtils::toDevice(m_queue, betaCells))
            , claddingCellTypes(hase::alpakaUtils::toDevice(m_queue, claddingCellTypes))
            , refractiveIndices(hase::alpakaUtils::toDevice(m_queue, refractiveIndices))
            , reflectivities(hase::alpakaUtils::toDevice(m_queue, reflectivities))
            , totalReflectionAngles(hase::alpakaUtils::toDevice(m_queue, totalReflectionAngles))
            , cellPointIndices(hase::alpakaUtils::toDevice(m_queue, cellPointIndices))
            , cellTypes(hase::alpakaUtils::toDevice(m_queue, cellTypes))
            , cellFaces(hase::alpakaUtils::toDevice(m_queue, cellFaces))
            , cellNeighborCells(hase::alpakaUtils::toDevice(m_queue, cellNeighborCells))
            , cellNeighborLocalFaces(hase::alpakaUtils::toDevice(m_queue, cellNeighborLocalFaces))
            , cellFaceBoundaries(hase::alpakaUtils::toDevice(m_queue, cellFaceBoundaries))
            , cellVolumes(hase::alpakaUtils::toDevice(m_queue, cellVolumes))
            , cellCenters(hase::alpakaUtils::toDevice(m_queue, cellCenters))
            , samplePoints(hase::alpakaUtils::toDevice(m_queue, samplePoints))
            , claddingAbsorption(claddingAbsorption)
            , nTot(nTot)
            , crystalTFluo(crystalTFluo)
            , claddingNumber(claddingNumber)
            , numberOfCells(numberOfCells)
            , numberOfPrisms(numberOfCells)
            , numberOfPoints(numberOfPoints)
            , numberOfSamples(numberOfSamples)
            , numberOfFacesPerCell(numberOfFacesPerCell)
            , numberOfCellVertices(numberOfCellVertices)
        {
        }

        [[nodiscard]] auto toView() const -> DeviceMeshView
        {
            return {
                std::span<double const>(points.data(), points.getMdSpan().getExtents().x()),
                std::span<double const>(betaVolume.data(), betaVolume.getMdSpan().getExtents().x()),
                std::span<double const>(betaCells.data(), betaCells.getMdSpan().getExtents().x()),
                std::span<unsigned const>(claddingCellTypes.data(), claddingCellTypes.getMdSpan().getExtents().x()),
                std::span<float const>(refractiveIndices.data(), refractiveIndices.getMdSpan().getExtents().x()),
                std::span<float const>(reflectivities.data(), reflectivities.getMdSpan().getExtents().x()),
                std::span<float const>(
                    totalReflectionAngles.data(),
                    totalReflectionAngles.getMdSpan().getExtents().x()),
                std::span<unsigned const>(cellPointIndices.data(), cellPointIndices.getMdSpan().getExtents().x()),
                std::span<unsigned const>(cellTypes.data(), cellTypes.getMdSpan().getExtents().x()),
                std::span<int const>(cellFaces.data(), cellFaces.getMdSpan().getExtents().x()),
                std::span<int const>(cellNeighborCells.data(), cellNeighborCells.getMdSpan().getExtents().x()),
                std::span<int const>(
                    cellNeighborLocalFaces.data(),
                    cellNeighborLocalFaces.getMdSpan().getExtents().x()),
                std::span<int const>(cellFaceBoundaries.data(), cellFaceBoundaries.getMdSpan().getExtents().x()),
                std::span<float const>(cellVolumes.data(), cellVolumes.getMdSpan().getExtents().x()),
                std::span<double const>(cellCenters.data(), cellCenters.getMdSpan().getExtents().x()),
                std::span<double const>(samplePoints.data(), samplePoints.getMdSpan().getExtents().x()),
                claddingAbsorption,
                nTot,
                crystalTFluo,
                claddingNumber,
                numberOfCells,
                numberOfCells,
                numberOfPoints,
                numberOfSamples,
                numberOfFacesPerCell,
                numberOfCellVertices};
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
        T_Buffer<double> betaCells;
        T_Buffer<unsigned> claddingCellTypes;
        T_Buffer<float> refractiveIndices;
        T_Buffer<float> reflectivities;
        T_Buffer<float> totalReflectionAngles;
        T_Buffer<unsigned> cellPointIndices;
        T_Buffer<unsigned> cellTypes;
        T_Buffer<int> cellFaces;
        T_Buffer<int> cellNeighborCells;
        T_Buffer<int> cellNeighborLocalFaces;
        T_Buffer<int> cellFaceBoundaries;
        T_Buffer<float> cellVolumes;
        T_Buffer<double> cellCenters;
        T_Buffer<double> samplePoints;

        double claddingAbsorption;
        float nTot;
        float crystalTFluo;
        unsigned claddingNumber;
        unsigned numberOfCells;
        unsigned numberOfPrisms;
        unsigned numberOfPoints;
        unsigned numberOfSamples;
        unsigned numberOfFacesPerCell;
        unsigned numberOfCellVertices;
    };

    class HostMesh
    {
    public:
        std::vector<double> points;
        std::vector<double> betaVolume;
        std::vector<double> betaCells;
        std::vector<unsigned> claddingCellTypes;
        std::vector<float> refractiveIndices;
        std::vector<float> reflectivities;
        std::vector<float> totalReflectionAngles;
        std::vector<unsigned> cellPointIndices;
        std::vector<unsigned> cellTypes;
        std::vector<int> cellFaces;
        std::vector<int> cellNeighborCells;
        std::vector<int> cellNeighborLocalFaces;
        std::vector<int> cellFaceBoundaries;
        std::vector<float> cellVolumes;
        std::vector<double> cellCenters;
        std::vector<double> samplePoints;
        float nTot = 0.0f;
        float crystalTFluo = 0.0f;
        unsigned claddingNumber = 1u;
        double claddingAbsorption = 0.0;
        unsigned numberOfCells = 0u;
        unsigned numberOfPrisms = 0u;
        unsigned numberOfPoints = 0u;
        unsigned numberOfSamples = 0u;
        unsigned numberOfFacesPerCell = tet4FaceCount;
        unsigned numberOfCellVertices = tet4VertexCount;
        bool resultAtVolumes = false;

        HostMesh() = default;

        HostMesh(
            std::vector<unsigned> cellPointIndices,
            std::vector<unsigned> cellTypes,
            std::vector<int> cellFaces,
            std::vector<int> cellNeighborCells,
            std::vector<int> cellNeighborLocalFaces,
            std::vector<int> cellFaceBoundaries,
            std::vector<float> cellVolumes,
            std::vector<double> points,
            std::vector<double> samplePoints,
            std::vector<double> cellCenters,
            std::vector<double> betaVolume,
            std::vector<double> betaCells,
            std::vector<unsigned> claddingCellTypes,
            std::vector<float> refractiveIndices,
            std::vector<float> reflectivities,
            float nTot,
            float crystalTFluo,
            unsigned claddingNumber,
            double claddingAbsorption)
            : points(std::move(points))
            , betaVolume(std::move(betaVolume))
            , betaCells(std::move(betaCells))
            , claddingCellTypes(std::move(claddingCellTypes))
            , refractiveIndices(std::move(refractiveIndices))
            , reflectivities(std::move(reflectivities))
            , cellPointIndices(std::move(cellPointIndices))
            , cellTypes(std::move(cellTypes))
            , cellFaces(std::move(cellFaces))
            , cellNeighborCells(std::move(cellNeighborCells))
            , cellNeighborLocalFaces(std::move(cellNeighborLocalFaces))
            , cellFaceBoundaries(std::move(cellFaceBoundaries))
            , cellVolumes(std::move(cellVolumes))
            , cellCenters(std::move(cellCenters))
            , samplePoints(std::move(samplePoints))
            , nTot(nTot)
            , crystalTFluo(crystalTFluo)
            , claddingNumber(claddingNumber)
            , claddingAbsorption(claddingAbsorption)
            , numberOfCells(static_cast<unsigned>(this->cellTypes.size()))
            , numberOfPrisms(numberOfCells)
            , numberOfPoints(static_cast<unsigned>(this->points.size() / 3u))
            , numberOfSamples(static_cast<unsigned>(this->samplePoints.size() / 3u))
        {
        }

        void calcTotalReflectionAngles()
        {
            std::vector<float> angles(refractiveIndices.size() / 2u, 0.0f);
            for(unsigned i = 0; i + 1u < refractiveIndices.size(); i += 2u)
            {
                angles.at(i / 2u)
                    = 180.0f / static_cast<float>(M_PI)
                      * alpaka::math::asin(refractiveIndices.at(i + 1u) / refractiveIndices.at(i));
            }
            totalReflectionAngles = std::move(angles);
        }

        template<typename T_Device>
        [[nodiscard]] DeviceMeshContainer<T_Device> toDevice(T_Device& device)
        {
            return DeviceMeshContainer<T_Device>{
                device,
                claddingAbsorption,
                nTot,
                crystalTFluo,
                claddingNumber,
                numberOfCells,
                numberOfPoints,
                numberOfSamples,
                numberOfFacesPerCell,
                numberOfCellVertices,
                points,
                betaVolume,
                betaCells,
                claddingCellTypes,
                refractiveIndices,
                reflectivities,
                totalReflectionAngles,
                 cellPointIndices,
                 cellTypes,
                 cellFaces,
                cellNeighborCells,
                cellNeighborLocalFaces,
                cellFaceBoundaries,
                 cellVolumes,
                 cellCenters,
                 samplePoints};
        }
    };

} // namespace hase::core
