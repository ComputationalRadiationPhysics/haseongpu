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
