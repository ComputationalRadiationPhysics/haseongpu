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
 *
 * @licence GPLv3
 **/

#pragma once
#include <alpaka/alpaka.hpp>

namespace hase::core
{

    struct TwoDimPoint
    {
        double x;
        double y;
    };

    typedef TwoDimPoint TwoDimDir;

    struct Point
    {
        double x;
        double y;
        double z;
    };

    typedef Point Vector;

    /**
     * @brief a Ray, defined by a startpoint, direction and length
     */
    struct Ray
    {
        Point p;
        Vector dir;
        double length;
    };

    struct NormalRay
    {
        TwoDimPoint p;
        TwoDimDir dir;
    };

    enum ReflectionPlane : int32_t
    {
        TOP_REFLECTION = 1,
        BOTTOM_REFLECTION = -1
    };

    ALPAKA_FN_HOST_ACC Vector direction(Point startPoint, Point endPoint);
    ALPAKA_FN_HOST_ACC double distance(Point startPoint, Point endPoint);
    ALPAKA_FN_HOST_ACC Ray generateRay(Point startPoint, Point endPoint);
    ALPAKA_FN_HOST_ACC Ray normalizeRay(Ray ray);

} // namespace hase::core
