#ifndef MESH_H
#define MESH_H 

typedef float4 Point;

struct Vector : float4 {
  __host__ __device__ float length();
  __host__ __device__ void normalize();
};

struct Ray : Point {
  Vector dir;

  __host__ __device__ float length();
  __host__ __device__ void normalize();
};

struct Edge {
  Ray normal;
  Triangle *neighbor; // index in global triangle list
};

struct Triangle {
  Point A;
  Point B;
  Point C;

  double *betaValues; // array of betaValues, one for each prism/level above the triangle
  Edge edges[3]; // edges of the triangle size: 3
};

#endif /* MESH_H */
