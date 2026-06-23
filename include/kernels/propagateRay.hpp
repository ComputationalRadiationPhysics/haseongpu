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
 * @brief Propagates a ray through the explicit 3D prism-cell structure.
 *        On each step the next cell on the ray path is calculated from
 *        explicit face and neighbor lookup tables.
 *        length and startpoint of propagation is stored inside the
 *        ray struct. The propagation ends when the length of the ray
 *        is reduced to zero. It is possible to do propagation with
 *        or without reflection. In case of reflection, the rays will
 *        have a longer way.
 *
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once
#include <core/geometry.hpp>
#include <core/mesh.hpp>

namespace hase::kernels
{

    /**
     * @brief Direct ray propagation without reflection
     *
     * @param ray           The ray which will propagate through the explicit cells.
     * @param startCell     The cell where the startpoint of the ray is located.
     * @param mesh          Explicit 3D cell mesh.
     * @param sigmaA        Absorption value of the ray.
     * @param sigmaE        Emission value of the ray.
     *
     * @return gain         Integral of ray length multiplied by beta values.
     *                      See at the code or accompanying paper for more clarity.
     *
     */
    ALPAKA_FN_ACC double propagateRay(
        core::Ray ray,
        unsigned* startCell,
        core::DeviceMeshView const& mesh,
        double sigmaA,
        double sigmaE);
    /**
     * @brief Compatibility wrapper for no-reflection explicit-cell propagation.
     *
     * @param startPoint      Point where the ray should start from.
     * @param endPoint        Point where the will end.
     * @param reflections      Number of reflections the ray will do
     *                        from startPoint to endPoint.
     * @param reflectionPlane Plane of first reflection (upper or lower surface of gain medium)
     * @param startCell       The cell where the startpoint of the ray is located.
     * @param mesh            Explicit 3D cell mesh.
     * @param sigmaA          Absorption value of the ray.
     * @param sigmaE          Emission value of the ray.
     *
     * @return gain           Integral of ray length multiplied by beta values.
     *                        See at the code or accompanying paper for more clarity.
     *
     */
    ALPAKA_FN_ACC double propagateRayWithReflection(
        core::Point startPoint,
        core::Point endPoint,
        unsigned reflections,
        core::ReflectionPlane reflectionPlane,
        unsigned startCell,
        core::DeviceMeshView const& mesh,
        double sigmaA,
        double sigmaE);

} // namespace hase::kernels
