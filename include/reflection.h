#ifndef Reflection_H
#define Reflection_H

#include "mesh.h"

enum ReflectionPlane {TOP_REFLECTION = 1, BOTTOM_REFLECTION = -1};

__device__ int calcNextReflection(const Point startPoint, const Point endPoint, const unsigned reflectionsLeft, const ReflectionPlane reflectionPlane, Point *reflectionPoint, double *reflectionAngle, const Mesh &mesh);

#endif /* Reflection_H */
