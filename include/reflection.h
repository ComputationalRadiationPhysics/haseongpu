#ifndef Reflection_H
#define Reflection_H

#include "mesh.h"

enum ReflectionPlane {TOP_REFLECTION = 1, BOTTOM_REFLECTION = -1};

__device__ int calcNextReflection(Point startPoint, Point endPoint, unsigned reflectionsLeft, ReflectionPlane reflectionPlane,Point *reflectionPoint, double *reflectionAngle, Mesh *mesh);

#endif /* Reflection_H */
