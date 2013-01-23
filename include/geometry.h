#ifndef GEOMETRY_H
#define GEOMETRY_H 

__device__ PointCu createPoint(float x, float y, float z, float w);
__device__ VectorCu createVector(PointCu a, PointCu b);
__device__ TriangleCu createTriangle(PointCu a, PointCu b, PointCu c);

__device__ float distance_gpu(PointCu a, PointCu b);
__device__ VectorCu crossproduct_gpu(VectorCu a, VectorCu b);
__device__ float skalar_mul_gpu(VectorCu a, VectorCu b);
__device__ PointCu intersectionRayTriangleGPU(PointCu rayOrigin, PointCu rayObjective, PointCu p1, PointCu p2,PointCu p3);
__device__ float collide_prism_gpu(PrismCu pr, RayCu r);
__global__ void trace_on_prisms(PrismCu* prisms, const unsigned max_prisms, RayCu* rays, const unsigned max_rays_per_sample, PointCu *samples, const unsigned blocks_per_sample);

#endif /* GEOMETRY_H */
