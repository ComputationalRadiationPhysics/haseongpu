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


#ifndef GEOMETRY_GPU_H
#define GEOMETRY_GPU_H

/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 *
 * @licence GPLv3
 **/

struct TwoDimPoint {
  double x;
  double y;
};

typedef TwoDimPoint TwoDimDir;

struct Point {
  double x;
  double y;
  double z;
};

typedef Point Vector;

/**
 * @brief a Ray, defined by a startpoint, direction and length
 */
struct Ray {
  Point p;
  Vector dir;
  float length;
};

struct NormalRay {
  TwoDimPoint p;
  TwoDimDir dir;
};

enum ReflectionPlane {TOP_REFLECTION = 1, BOTTOM_REFLECTION = -1};

__host__ __device__ Vector direction(Point startPoint, Point endPoint);
__host__ __device__ float distance(Point startPoint, Point endPoint);
__host__ __device__ Ray generateRay(Point startPoint, Point endPoint);
__host__ __device__ Ray normalizeRay(Ray ray);

#endif /* GEOMETRY_GPU_H */
