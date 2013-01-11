#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include <vector>

//----------------------------------------------------
// Structures
//----------------------------------------------------
typedef struct point {
  float x;
  float y;
  float z;
} POINT;

typedef struct vector {
  float x;
  float y;
  float z;
} VECTOR;

typedef struct ray {
  point start;
  vector direction;
} RAY;

typedef struct triangle {
  point a;
  point b;
  point c;
} TRIANGLE;

typedef struct plane {
  point start;
  vector normal;
  
} PLANE;

//------------------------------------------
typedef float4 point_cu;
typedef float4 vector_cu;

typedef struct triangle_cu{
  point_cu A;
  point_cu B;
  point_cu C;
} TRIANGLE_CU;

typedef struct prism_cu{
  triangle_cu t1;
  triangle_cu t2;
} PRISM_CU;

typedef struct plane_cu {
  point_cu P;
  vector_cu normal;
} PLANE_CU;

typedef struct ray_par {
  float wavelength, b, c, d;
} RAY_PAR;

typedef struct ray_cu {
  point_cu P;
  vector_cu direction;
  ray_par data;
} RAY_CU;

//----------------------------------------------------
// Auxillary function declaration
//----------------------------------------------------

float distance(point a, point b);
void  print_point(point p);

// New functions
bool  collide(triangle_cu t, point_cu p);
bool  collide(triangle_cu t, ray_cu r);
float4 to_barycentric(triangle_cu t, point_cu p);
point_cu intersection(plane_cu p, ray_cu r);
std::vector<triangle_cu> generate_triangles(int height, int width, float level);
std::vector<ray_cu> generate_rays(int height, int width, float level, unsigned max_rays);
ray_cu   generate_ray(int height, int weight, float level);

//----------------------------------------------------
//  Calculations
//----------------------------------------------------
__device__ float4 to_barycentric_gpu(triangle_cu t, point_cu p){
  float x1,x2,x3, y1,y2,y3, x,y;
  float4 b;

  x1 = t.A.x;
  x2 = t.B.x;
  x3 = t.C.x;

  y1 = t.A.y;
  y2 = t.B.y;
  y3 = t.C.y;
    
  x = p.x;
  y = p.y;

  b.x = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3)) / ((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
  b.y = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3)) / ((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
  b.z = 1 - b.x - b.y;
  b.w = 0;
  return b;
}

__device__ bool collide_gpu(triangle_cu t, point_cu p){
  float4 b = to_barycentric_gpu(t, p);
  bool has_collided = true;
  if(b.x < 0 || b.x > 1 || b.y < 0 || b.y > 1 || b.z < 0 || b.z > 1) has_collided = false;
  return has_collided;
}

__device__ point_cu intersection_gpu(plane_cu pl, ray_cu r){
  point_cu intersection_point = {0.0,0.0,0.0};

  float t, d;

  // vector coordinates
  float n1, n2, n3, x1, x2, x3, p1, p2, p3, a1, a2, a3;
  
  // just get the coordinates from the structs
  n1 = pl.normal.x;
  n2 = pl.normal.y;
  n3 = pl.normal.z;

  a1 = pl.P.x;
  a2 = pl.P.y;
  a3 = pl.P.z;

  x1 = r.P.x;
  x2 = r.P.y;
  x3 = r.P.z;

  p1 = r.direction.x;
  p2 = r.direction.y;
  p3 = r.direction.z;

  // calculation of intersection
  d = n1*a1 + n2*a2 + n3*a3;
  t = (d - n1*x1 - n2*x2 - n3*x3) / (n1*p1 + n2*p2 + n3*p3);

  intersection_point.x = x1 + t * p1;
  intersection_point.y = x2 + t * p2;
  intersection_point.z = x3 + t * p3;

  return intersection_point;

}


__global__ void trace_on_gpu(triangle_cu* triangles, const unsigned max_triangles, ray_cu* rays, const unsigned max_rays, float4 *collisions)
{
  // Cuda ids
  unsigned tid = threadIdx.x;
  unsigned bid = blockIdx.x + blockIdx.y * gridDim.x;
  unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;

  // Local data
  const plane_cu plane = {
    {0, 0, 0, 1},
    {0, 0, 1, 0}};
  point_cu intersection_point;
  triangle_cu triangle = triangles[gid];
  unsigned ray_i;

  __syncthreads();
  // Calculation
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
    intersection_point = intersection_gpu(plane, rays[ray_i]);
    if(collide_gpu(triangle, intersection_point)){
      collisions[gid].x++;
    }
    __syncthreads();
    
  }
  
}

int main(){
  const unsigned max_rays = 2;
  const unsigned max_triangles =  100;
  unsigned length = ceil(sqrt(max_triangles / 2));
  unsigned i;
  unsigned ray_i;
  unsigned triangle_i;
  ray_cu r;
  point_cu p;
  float runtime_gpu = 0.0;
  float runtime_cpu = 0.0;
  cudaEvent_t start, stop;
  const plane_cu pl = {
    {0, 0, 0, 1},
    {0, 0, 1, 0}};

  // Generate testdata
  std::vector<triangle_cu> triangles = generate_triangles(length, length, 0);
  std::vector<ray_cu> rays = generate_rays(length, length, 0, max_rays);
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // CPU Raytracing
  cudaEventRecord(start, 0);
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
    p = intersection(pl, rays[ray_i]);
    for(triangle_i = 0; triangle_i < triangles.size(); ++triangle_i){
      if(collide(triangles[triangle_i], p)){
  	fprintf(stdout, "CPU: Ray %d hits on Triangle %d\n", ray_i, triangle_i);
  	// Do something on collision
      }

      }

  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&runtime_cpu, start, stop);

  // GPU Raytracing
  triangle_cu* h_triangles, *d_triangles;
  ray_cu* h_rays, *d_rays;
  float4* h_collisions, *d_collisions;
  int threads = 32;
  int blocks = ceil(max_triangles / threads);

  // Memory allocation on host
  cudaHostAlloc( (void**)&h_triangles, max_triangles * sizeof(triangle_cu), cudaHostAllocDefault);
  cudaHostAlloc( (void**)&h_rays, max_rays * sizeof(ray_cu), cudaHostAllocDefault);
  cudaHostAlloc( (void**)&h_collisions, max_triangles * sizeof(float4), cudaHostAllocDefault);

  // Memory initialisation on host
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
    h_rays[i] = rays[i];
  }

  for(triangle_i = 0; triangle_i < max_triangles; ++triangle_i){
    h_collisions[triangle_i].x = 0;
    h_triangles[triangle_i] = triangles[triangle_i];
  }

  // Memory allocation on device
  cudaMalloc(&d_triangles, max_triangles * sizeof(triangle_cu));
  cudaMalloc(&d_rays, max_rays * sizeof(ray_cu));
  cudaMalloc(&d_collisions, max_triangles * sizeof(float4));

  // Copy data from host to device
  cudaEventRecord(start, 0);
  cudaMemcpy(d_triangles, h_triangles, max_triangles * sizeof(triangle_cu), cudaMemcpyHostToDevice);
  cudaMemcpy(d_rays, h_rays, max_rays * sizeof(ray_cu), cudaMemcpyHostToDevice);
  cudaMemcpy(d_collisions, h_collisions, max_triangles * sizeof(float4), cudaMemcpyHostToDevice);

  // Start kernel
  trace_on_gpu<<<threads, blocks>>>(d_triangles, max_triangles, d_rays, max_rays, d_collisions);

  // Copy data from device to host
  cudaMemcpy(h_collisions, d_collisions, max_triangles * sizeof(float4), cudaMemcpyDeviceToHost);

  // Evaluate device data
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&runtime_gpu, start, stop);
  for(triangle_i = 0; triangle_i < max_triangles; ++triangle_i){
    if(h_collisions[triangle_i].x > 0){
      fprintf(stderr, "GPU: %d collisions on triangle %d\n", (int)h_collisions[triangle_i].x, triangle_i );
    }
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "Triangles   : %d\n", max_triangles);
  fprintf(stderr, "Rays        : %d\n", max_rays);
  fprintf(stderr, "GPU Blocks  : %d\n", blocks);
  fprintf(stderr, "GPU Threads : %d\n", threads);
  fprintf(stderr, "Runtime_GPU : %f\n", runtime_gpu);
  fprintf(stderr, "Runtime_CPU : %f\n", runtime_cpu);
  fprintf(stderr, "\n");

  // Cleanup
  cudaFreeHost(h_triangles);
  cudaFreeHost(h_rays);
  cudaFreeHost(h_collisions);

  return 0;
}

//----------------------------------------------------
// Auxillary function definition
//----------------------------------------------------

float4 to_barycentric(triangle_cu t, point_cu p){
  float x1,x2,x3, y1,y2,y3, x,y;
  float4 b;

  x1 = t.A.x;
  x2 = t.B.x;
  x3 = t.C.x;

  y1 = t.A.y;
  y2 = t.B.y;
  y3 = t.C.y;
    
  x = p.x;
  y = p.y;

  b.x = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3)) / ((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
  b.y = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3)) / ((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
  b.z = 1 - b.x - b.y;
  b.w = 0;
  assert(fabs((b.x + b.y + b.z) - 1) < 0.0001);
  return b;
}

/**
   @brief Detects collisions of triangle and point with
   precondition, that the point is on the same 
   plane as the point.
**/
bool collide(triangle_cu t, point_cu p){
  float4 b = to_barycentric(t, p);
  bool has_collided = true;
  if(b.x < 0 || b.x > 1) has_collided = false;
  if(b.y < 0 || b.y > 1) has_collided = false;
  if(b.z < 0 || b.z > 1) has_collided = false;

  return has_collided;
}

/**
   @brief Detects collisions of a triangle and a ray without
   a precondition.
**/
bool collide(triangle_cu t, ray_cu r){
  plane_cu pl;
  float b1, b2, b3, c1, c2, c3;

  b1 = t.B.x;
  b2 = t.B.y;
  b3 = t.B.z;

  c1 = t.C.x;
  c2 = t.C.y;
  c3 = t.C.z;

  pl.P = t.A;
  pl.normal.x = (b2*c3 - b3*c2);
  pl.normal.y = (b3*c1 - b1*c3);
  pl.normal.z = (b1*c2 - b2*c1);

  return collide(t, intersection(pl, r));
}

/**
   @brief Intersection calculates the intersection between a plane p
   and a ray r. There is no detection for rays in the plane
   or for parallel plane. 

   It uses the normal of the plane to derive the coordinate form 
   of the plane. With the help of a coordinate form it is very
   easy to get the intersection point between a ray and a plane.

   ray   g: y~ = x~ + t*p~
   plane E: y~ = a~ + r*b~ + s*c~
   d  = n1*(x1+t*p1) + n2*(x2+t*p2) + n3*(x3+t*p3)
   d  = n~ * a~
**/
point_cu intersection(plane_cu pl, ray_cu r){
  point_cu intersection_point = {0.0,0.0,0.0};

  float t, d;

  // vector coordinates
  float n1, n2, n3, x1, x2, x3, p1, p2, p3, a1, a2, a3;
  
  // just get the coordinates from the structs
  n1 = pl.normal.x;
  n2 = pl.normal.y;
  n3 = pl.normal.z;

  a1 = pl.P.x;
  a2 = pl.P.y;
  a3 = pl.P.z;

  x1 = r.P.x;
  x2 = r.P.y;
  x3 = r.P.z;

  p1 = r.direction.x;
  p2 = r.direction.y;
  p3 = r.direction.z;

  // calculation of intersection
  d = n1*a1 + n2*a2 + n3*a3;
  t = (d - n1*x1 - n2*x2 - n3*x3) / (n1*p1 + n2*p2 + n3*p3);

  intersection_point.x = x1 + t * p1;
  intersection_point.y = x2 + t * p2;
  intersection_point.z = x3 + t * p3;

  return intersection_point;

}

float distance(point a, point b){
  float d = sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y),2) + pow((b.z - a.z),2));
  return fabs(d);
}

std::vector<triangle_cu> generate_triangles(int height, int weight, float level){
  int h,w;
  std::vector<triangle_cu> triangles;
  for(h = 0; h < height; ++h){
    for(w = 0; w < weight; ++w){
      triangle_cu t1 = {
	{float(h), float(w), level, 1},
	{float(h), float(w+1), level, 1},
	{float(h+1), float(w), level, 1}};
      triangle_cu t2 = {
	{float(h), float(w+1), level, 1},
	{float(h+1), float(w+1), level, 1},
	{float(h+1), float(w), level, 1}};
      triangles.push_back(t1);
      triangles.push_back(t2);

    }

  }

  return triangles;
}

ray_cu generate_ray(const int height, const int width, const float level){
  float random1 = float(rand() % height) + (rand() / (float) RAND_MAX);
  float random2 = float(rand() % width) + (rand() / (float) RAND_MAX);
  ray_cu r = {
    {random1, random2, level, 1},
    {0,0,1, 0}};
  return r;
}

std::vector<ray_cu> generate_rays(const int height, const int width, const float level, const unsigned max_rays){
  std::vector<ray_cu> rays;
  unsigned ray_i;
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
    ray_cu ray = generate_ray(height, width, level);
    rays.push_back(ray);
  }
  return rays;
}

void print_point(point p){
  fprintf(stdout, "Point\n");
  fprintf(stdout, "x: %f\n", p.x);
  fprintf(stdout, "y: %f\n", p.y);
  fprintf(stdout, "z: %f\n", p.z);

}
