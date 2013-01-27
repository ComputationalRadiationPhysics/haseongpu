#ifndef GEOMETRIE_H
#define GEOMETRIE_H

#include "datatypes.h"

PointCu  collide_triangle(TriangleCu t, RayCu r);
float  collide_prism(PrismCu pr, RayCu r);
float4 to_barycentric(TriangleCu t, RayCu r);
PointCu intersection(PlaneCu p, RayCu r);
float distance(PointCu a, PointCu b);
VectorCu crossproduct(VectorCu a, VectorCu b);
float skalar_mul(VectorCu a, VectorCu b);

#endif /* GEOMETRIE_H */
