/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASENonGPU
 *
 * HASENonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASENonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASENonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @licence GPLv3
 *
 */

#ifndef Reflection_H
#define Reflection_H

#include <mesh.h>
#include <geometry.h>

/**
 * @brief Calculates the reflectionPoint and reflectionAngle with upper or lower surface
 *        of gain medium. Depending on the number of reflection the
 *        ray still has to do and the reflection surface. 
 *
 * @param startPoint      Point where the ray should start from.
 * @param endPoint        Point where the will end.
 * @param reflectionsLeft Number of reflections the ray will do
 *                        from startPoint to endPoint.
 * @param reflectionPlane Plane of first reflection (upper or lower surface of gain medium).
 * @param reflectionPoint Pointer to the point the ray will be reflected.
 * @param reflectionAngle Pointer to the value of the angle the ray will be reflected with.
 *
 * @return 0 if everything goes right
 *         1 otherwise
 *
 */
__device__ int calcNextReflection(const Point startPoint, 
				  const Point endPoint, 
				  const unsigned reflectionsLeft, 
				  const ReflectionPlane reflectionPlane, 
				  Point *reflectionPoint, 
				  double *reflectionAngle, 
				  const Mesh &mesh);

#endif /* Reflection_H */
