#ifndef GEOMETRY_GPU_H
#define GEOMETRY_GPU_H

#include "mesh.h"

/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 *
 * @licence GPLv3
 **/

__host__ __device__ Vector direction(Point startPoint, Point endPoint);
__host__ __device__ float distance(Point startPoint, Point endPoint);
__host__ __device__ Ray generateRay(Point startPoint, Point endPoint);
__host__ __device__ Ray normalizeRay(Ray ray);

#endif /* GEOMETRY_GPU_H */
