#ifndef GEOMETRY_GPU_H
#define GEOMETRY_GPU_H

#include "datatypes.h"
#include "curand_kernel.h"

__device__ PointCu addVectorToPoint(PointCu p, VectorCu v);
__device__ float    collide_prism_gpu(PrismCu pr, RayCu r, VectorCu rayDirection, double absRayDistance);
__device__ float    distance_gpu(PointCu a, PointCu b);
__device__ VectorCu crossproduct_gpu(VectorCu a, VectorCu b);
__device__ float    skalar_mul_gpu(VectorCu a, VectorCu b);
__device__ double  intersectionRayTriangleGPU(PointCu rayOrigin, VectorCu rayDirection, PointCu p1, PointCu p2,PointCu p3);
__device__ RayCu    generateRayGpu(PointCu sample, PrismCu startPrism, curandState localState);
__device__ PrismCu  selectPrism(int id, PrismCu *prisms, int totalNumberOfPrisms);


#endif /* GEOMETRY_GPU_H */
