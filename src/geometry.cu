#include "mesh.h"

__host__ __device__ Vector direction(Point startPoint, Point endPoint){
  Vector v = {endPoint.x - startPoint.x, endPoint.y - startPoint.y, endPoint.z - startPoint.z};
  return v;
}

__host__ __device__ float distance(Point startPoint, Point endPoint){
  float x, y, z;
  x = endPoint.x - startPoint.x;
  y = endPoint.y - startPoint.y;
  z = endPoint.z - startPoint.z;
  float d = sqrt(x*x + y*y + z*z);
  return fabs(d);

}

__host__ __device__ Ray generateRay(Point startPoint, Point endPoint){
  Ray ray;
  ray.p = startPoint;
  ray.dir = direction(startPoint, endPoint);
  ray.length = distance(startPoint, endPoint);

  return ray;

}

__host__ __device__ Ray normalizeRay(Ray ray){
  ray.dir.x = ray.dir.x / ray.length;
  ray.dir.y = ray.dir.y / ray.length;
  ray.dir.z = ray.dir.z / ray.length;
  ray.length = 1;

  return ray;
}
