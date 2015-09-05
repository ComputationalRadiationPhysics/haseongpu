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



// ALPAKA
#include <alpaka/alpaka.hpp> /* ALPAKA_FN_HOST_ACC */

// HASEonGPU
#include <geometry.hpp> /* Vector, Point, Ray */

ALPAKA_FN_HOST_ACC Vector direction(Point startPoint, Point endPoint){
  Vector v = {endPoint.x - startPoint.x, endPoint.y - startPoint.y, endPoint.z - startPoint.z};
  return v;
}

ALPAKA_FN_HOST_ACC float distance(Point startPoint, Point endPoint){
  float x, y, z;
  x = endPoint.x - startPoint.x;
  y = endPoint.y - startPoint.y;
  z = endPoint.z - startPoint.z;
  float d = sqrt(x*x + y*y + z*z);
  return fabs(d);

}

ALPAKA_FN_HOST_ACC Ray generateRay(Point startPoint, Point endPoint){
  Ray ray;
  ray.p = startPoint;
  ray.dir = direction(startPoint, endPoint);
  ray.length = distance(startPoint, endPoint);

  return ray;

}

ALPAKA_FN_HOST_ACC Ray normalizeRay(Ray ray){
  ray.dir.x = ray.dir.x / ray.length;
  ray.dir.y = ray.dir.y / ray.length;
  ray.dir.z = ray.dir.z / ray.length;
  ray.length = 1;

  return ray;
}
