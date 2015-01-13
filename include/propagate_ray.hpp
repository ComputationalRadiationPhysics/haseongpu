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


/**
 * @brief Propagates an ray through the triangle/prism/crystal structure.
 *        On each step the next triangle on the ray path will be calculated 
 *        from the current triangle (startTriangle in the beginning).
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
#include <mesh.hpp>
#include <reflection.hpp> /* ReflectionPlane */

/**
 * @brief Direct ray propagation without reflection
 *
 * @param ray           The ray which will propagate through the prisms.
 * @param startLevel    The level where the startpoint of the ray is located.
 * @param startTriangle The triangle where the startpoint of the ray is locaed
 * @param mesh          All information about triangles, points, contants. 
 *                      See mesh.h for details.
 * @param sigmaA        Absorption value of the ray.
 * @param sigmaE        Emission value of the ray.
 *
 * @return gain         Integral of ray length multiplied by beta values.
 *                      See at the code or accompanying paper for more clarity.
 *
 */
__device__ double propagateRay(Ray ray, 
			       unsigned *startLevel, 
			       unsigned  *startTriangle, 
			       const Mesh &mesh,
			       const double sigmaA, 
			       const double sigmaE
			       );
/**
 * @brief Indirect ray propagation with reflections on upper and lower surface
 *
 * @param startPoint      Point where the ray should start from.
 * @param endPoint        Point where the will end.
 * @param reflections      Number of reflections the ray will do
 *                        from startPoint to endPoint.
 * @param reflectionPlane Plane of first reflection (upper or lower surface of gain medium)
 * @param startLevel      The level where the startpoint of the ray is located.
 * @param startTriangle   The triangle where the startpoint of the ray is locaed
 * @param mesh            All information about triangles, points, contants. 
 *                        See mesh.h for details.
 * @param sigmaA          Absorption value of the ray.
 * @param sigmaE          Emission value of the ray.
 *
 * @return gain           Integral of ray length multiplied by beta values.
 *                        See at the code or accompanying paper for more clarity.
 *
 */
__device__ double propagateRayWithReflection(Point startPoint, 
					     const Point endPoint, 
					     const unsigned reflections,
					     ReflectionPlane reflectionPlane,
					     unsigned startLevel, 
					     unsigned startTriangle, 
					     const Mesh &mesh, 
					     const double sigmaA, 
					     const double sigmaE);



