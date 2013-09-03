#ifndef Reflection_H
#define Reflection_H

#include "mesh.h"
__device__ int calcNextReflection(Point startPoint, Point endPoint, unsigned reflectionsLeft, int reflectionPlane, Point *reflectionPoint, float *reflectionAngle, Mesh *mesh);

#endif /* Reflection_H */
