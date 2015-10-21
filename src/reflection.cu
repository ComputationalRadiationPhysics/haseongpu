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


#include <cassert>
#include <cmath> /* M_PI */

#include <reflection.hpp>
#include <mesh.hpp>
#include <geometry.hpp>

__device__ double calcIntersectionAngle(const Ray ray, double *reflectionAngle){
  // Calc intesection angle with z-plane
  double nominator = abs(ray.dir.z);
  double denominator = sqrt((ray.dir.x * ray.dir.x) + (ray.dir.y * ray.dir.y) + (ray.dir.z * ray.dir.z));
  if(denominator != 0.0){
    double radian = acos(nominator / denominator);
    *reflectionAngle = ((180. / M_PI) * radian);
    return 0;
  }
  return 1;
}

__device__ int calcPlaneIntersectionPoint(const Ray reflectionRay, const ReflectionPlane reflectionPlane, const Mesh &mesh, Point *intersectionPoint){
  // Assume that mesh is on x/y axis and parallel to x/y axis
  double planeZ = 0.0;
  if(reflectionPlane == TOP_REFLECTION){
    // Reflection on TOP plane
    planeZ = mesh.thickness * mesh.numberOfLevels;
  }
  double denominator = reflectionRay.dir.z;  
  if (denominator != 0.0){
    double nominator = planeZ - reflectionRay.p.z;
    double length = nominator / denominator;
    if(length > 0){
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
__device__ Ray generateReflectionRay(const Point startPoint, Point endPoint,  const int reflectionsLeft, const ReflectionPlane reflectionPlane, const Mesh &mesh){
  float mirrorPlaneZ = 0;
  if(reflectionsLeft % 2 == 0){
    // Even reflectionCount is postponement
    endPoint.z = endPoint.z + reflectionPlane * (reflectionsLeft * mesh.thickness * mesh.numberOfLevels); 
  }
  else {
    // Odd reflectionsCount is reflection

    if(reflectionPlane == TOP_REFLECTION){
      mirrorPlaneZ = ceil(reflectionsLeft/(double)2) * mesh.thickness * mesh.numberOfLevels;
    }
    else{
      mirrorPlaneZ = floor(reflectionsLeft/(double)2) * mesh.thickness * mesh.numberOfLevels * reflectionPlane;
    }

    endPoint.z = reflectionPlane * abs(( mirrorPlaneZ + mirrorPlaneZ - endPoint.z));
    
  }
  return generateRay(startPoint, endPoint);
}

__device__ int calcNextReflection(const Point startPoint, const Point endPoint, const unsigned reflectionsLeft, const ReflectionPlane reflectionPlane, Point *reflectionPoint, double *reflectionAngle, const Mesh &mesh){
  Ray reflectionRay = generateReflectionRay(startPoint, endPoint, reflectionsLeft, reflectionPlane, mesh);
  if(calcPlaneIntersectionPoint(reflectionRay, reflectionPlane, mesh, reflectionPoint)) return 1;
  if(calcIntersectionAngle(reflectionRay, reflectionAngle)) return 1;

  return 0;
}
