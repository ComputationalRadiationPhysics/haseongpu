#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include <vector>

#define CUDA_CHECK_RETURN(value) {				\
  cudaError_t _m_cudaStat = value;				\
  if (_m_cudaStat != cudaSuccess) {				\
    fprintf(stderr, "Error %s at line %d in file %s\n",	\
    cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);	\
    exit(1);							\
  }								\
}

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
} PRISM_CU;

typedef struct plane_cu {
  point_cu P;
  vector_cu normal;
} PLANE_CU;

typedef struct ray_cu {
  point_cu P;
  vector_cu direction;
} RAY_CU;

//----------------------------------------------------
// Auxillary function declaration
//----------------------------------------------------

float distance(point a, point b);
void  print_point(point p);

// New functions
bool  collide(triangle_cu t, point_cu p);
bool  collide(triangle_cu t, ray_cu r);
bool  collide(prism_cu pr, ray_cu r);
float4 to_barycentric(triangle_cu t, point_cu p);
point_cu intersection(plane_cu p, ray_cu r);
std::vector<triangle_cu> generate_triangles(int height, int width, float level);
std::vector<prism_cu> generate_prisms(int height, int width, float level);
std::vector<ray_cu> generate_rays(int height, int width, int level, unsigned max_rays);
ray_cu   generate_ray(int height, int weight, int level);

//----------------------------------------------------
// Device Code
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

__device__ point_cu intersection_gpu(plane_cu pl, ray_cu r){
  point_cu intersection_point = {0.0,0.0,0.0};

  float t = 0;
  float d = 0;

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

__device__ bool collide_gpu(triangle_cu t, point_cu p){
  float4 b = to_barycentric_gpu(t, p);
  return  (b.x > 0) && (b.x < 1) && (b.y > 0) && (b.y < 1) && (b.z > 0) && (b.z < 1) && (b.z == b.z);
}

__device__ bool collide_gpu(triangle_cu t, ray_cu r){
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

  return collide_gpu(t, intersection_gpu(pl, r));
}

__device__ bool collide_gpu(prism_cu pr, ray_cu r){
  bool has_collide = false;
  point_cu A1 = pr.t1.A;
  point_cu B1 = pr.t1.B;
  point_cu C1 = pr.t1.C;
  point_cu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.w, 1};
  point_cu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.w, 1};
  point_cu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.w, 1};

  triangle_cu triangles[8] = {
    pr.t1,
    {A2, B2, C2},
    {A1, B1, A2},
    {B1, B2, A2},
    {B1, C1, C2},
    {B1, B2, C2},
    {A1, C1, C2},
    {A1, A2, C2}};

  has_collide = has_collide || collide_gpu(triangles[0], r);
  has_collide = has_collide || collide_gpu(triangles[1], r);
  has_collide = has_collide || collide_gpu(triangles[2], r); 
  has_collide = has_collide || collide_gpu(triangles[3], r);
  has_collide = has_collide || collide_gpu(triangles[4], r); 
  has_collide = has_collide || collide_gpu(triangles[5], r); 
  has_collide = has_collide || collide_gpu(triangles[6], r); 
  has_collide = has_collide || collide_gpu(triangles[7], r);

  return has_collide;
}

__global__ void trace_on_prisms(prism_cu* prisms, const unsigned max_prisms, ray_cu* rays, const unsigned max_rays, float4 *collisions){
  // Cuda ids
  //unsigned tid = threadIdx.x;
  //unsigned bid = blockIdx.x + blockIdx.y * gridDim.x;
  
  unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;

  // Local data
  prism_cu prism = prisms[gid];
  unsigned local_collisions = 0;
  unsigned ray_i;

  __syncthreads();
  // Calculation
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
     local_collisions += (1 * (int)collide_gpu(prism, rays[ray_i]));
  
	       
  }
  __syncthreads();
  collisions[gid].x = local_collisions;
  
}
//*/


//----------------------------------------------------
// Host Code
//----------------------------------------------------
int main(){
  const unsigned max_rays = 1000;
  const unsigned max_triangles = 1000;
  const unsigned length = ceil(sqrt(max_triangles / 2));
  const unsigned depth  = 1;
  const unsigned max_prisms = length * length * depth * 2;
  unsigned ray_i, prism_i;
  float runtime_gpu = 0.0;
  float runtime_cpu = 0.0;
  cudaEvent_t start, stop;
  bool use_cpu = true;
  bool use_gpu = true;

  // Generate testdata
  std::vector<prism_cu> prisms = generate_prisms(length, length, depth);
  std::vector<ray_cu> rays = generate_rays(length, length, depth, max_rays);
  std::vector<float> collisions(max_prisms, 0);
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // CPU Raytracing
  cudaEventRecord(start, 0);
  if(use_cpu){
    for(ray_i = 0; ray_i < rays.size(); ++ray_i){
      for(prism_i = 0; prism_i < prisms.size(); ++prism_i){
	if(collide(prisms[prism_i], rays[ray_i])){
	  fprintf(stdout, "CPU: Ray %d hits on prism %d\n", ray_i, prism_i);
	  collisions[prism_i]++;
	}

      }
    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&runtime_cpu, start, stop);
  }

  // GPU Raytracing
  ray_cu* h_rays, *d_rays;
  prism_cu* h_prisms, *d_prisms;
  float4* h_collisions, *d_collisions;
  int threads = 256;
  int blocks = ceil(max_prisms / threads);
  if(use_gpu){

    // Memory allocation on host
    CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_prisms, max_prisms * sizeof(prism_cu), cudaHostAllocDefault));
    CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_rays, max_rays * sizeof(ray_cu), cudaHostAllocDefault));
    CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_collisions, max_prisms * sizeof(float4), cudaHostAllocDefault));

    // Memory initialisation on host
    for(ray_i = 0; ray_i < max_rays; ++ray_i){
      h_rays[ray_i] = rays[ray_i];
    }
    for(prism_i = 0; prism_i < max_prisms; ++prism_i){
      h_collisions[prism_i].x = 0;
      h_prisms[prism_i] = prisms[prism_i];
    }

    // Memory allocation on device
    CUDA_CHECK_RETURN(cudaMalloc(&d_rays, max_rays * sizeof(ray_cu)));
    CUDA_CHECK_RETURN(cudaMalloc(&d_prisms, max_prisms * sizeof(prism_cu)));
    CUDA_CHECK_RETURN(cudaMalloc(&d_collisions, max_prisms * sizeof(float4)));

    // Copy data from host to device
    cudaEventRecord(start, 0);
    CUDA_CHECK_RETURN(cudaMemcpy(d_rays, h_rays, max_rays * sizeof(ray_cu), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(d_prisms, h_prisms, max_prisms * sizeof(prism_cu), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(d_collisions, h_collisions, max_prisms * sizeof(float4), cudaMemcpyHostToDevice));

    // Start kernel
    trace_on_prisms<<<threads, blocks>>>(d_prisms, max_prisms, d_rays, max_rays, d_collisions);

    // Copy data from device to host
    CUDA_CHECK_RETURN(cudaMemcpy(h_collisions, d_collisions, max_prisms * sizeof(float4), cudaMemcpyDeviceToHost));

    // Evaluate device data
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&runtime_gpu, start, stop);
    for(prism_i = 0; prism_i < max_prisms; ++prism_i){
      if(h_collisions[prism_i].x > 0)
	fprintf(stderr, "GPU: (%f, %f, %f, %f) collission on prism %d\n", h_collisions[prism_i].x, h_collisions[prism_i].y, h_collisions[prism_i].z, h_collisions[prism_i].w, prism_i);

    }
    for(prism_i = 0; prism_i < max_prisms; ++prism_i){
      if((h_collisions[prism_i].x != collisions[prism_i]) && use_cpu && use_gpu){
	fprintf(stderr, "\033[31;1m[Error]\033[m CPU(%.0f) != GPU(%.0f) on prism %d\n",collisions[prism_i], h_collisions[prism_i].x, prism_i);
      }
    }
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Prism       : %d\n", max_prisms);
  fprintf(stderr, "Triangles   : %d\n", max_prisms * 8);
  fprintf(stderr, "Rays        : %d\n", max_rays);
  fprintf(stderr, "GPU Blocks  : %d\n", blocks);
  fprintf(stderr, "GPU Threads : %d\n", threads);
  fprintf(stderr, "Runtime_GPU : %f s\n", runtime_gpu / 1000.0);
  fprintf(stderr, "Runtime_CPU : %f s\n", runtime_cpu / 1000.0);
  fprintf(stderr, "\n");

  // Cleanup
  cudaFreeHost(h_rays);
  cudaFreeHost(h_prisms);
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

  // In case of division by 0 --> nan
  if((fabs((b.x + b.y + b.z) - 1)) != (fabs((b.x + b.y + b.z) - 1)))
    b.z = 2;
  return b;
}

/**
   @brief Detects collisions of triangle and point with
   precondition, that the point is on the same 
   plane as the point.
**/
bool collide(triangle_cu t, point_cu p){
  float4 b = to_barycentric(t, p);
  return (b.x > 0) && (b.x < 1) && (b.y > 0) && (b.y < 1) && (b.z > 0) && (b.z < 1) && (b.z == b.z);
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

bool collide(prism_cu pr, ray_cu r){
  bool has_collide;
  point_cu A1 = pr.t1.A;
  point_cu B1 = pr.t1.B;
  point_cu C1 = pr.t1.C;
  point_cu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.w, 1};
  point_cu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.w, 1};
  point_cu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.w, 1};

  triangle_cu triangles[8] = {
    pr.t1,
    {A2, B2, C2},
    {A1, B1, A2},
    {B1, B2, A2},
    {B1, C1, C2},
    {B1, B2, C2},
    {A1, C1, C2},
    {A1, A2, C2}};

  has_collide = 
    collide(triangles[0], r)
    || collide(triangles[1], r)
    || collide(triangles[2], r) 
    || collide(triangles[3], r)
    || collide(triangles[4], r) 
    || collide(triangles[5], r) 
    || collide(triangles[6], r) 
    || collide(triangles[7], r);

  return has_collide;
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

std::vector<prism_cu> generate_prisms(int height, int weight, float level){
  int h,w,l;
  std::vector<prism_cu> prisms;
  for(l = 0; l < level; ++l){
    for(h = 0; h < height; ++h){
      for(w = 0; w < weight; ++w){
	triangle_cu a1 = {
	  {float(h), float(w), l, l+1},
	  {float(h), float(w+1), l, l+1},
	  {float(h+1), float(w), l, l+1}};
	triangle_cu b1 = {
	  {float(h), float(w+1), l, 1+1},
	  {float(h+1), float(w+1), l, 1+1},
	  {float(h+1), float(w), l, 1+1}};
      
	prism_cu pr1 = {a1};
	prism_cu pr2 = {b1};

	prisms.push_back(pr1);
	prisms.push_back(pr2);

      }

    }

  }

  return prisms;
}

ray_cu generate_ray(const int heigth, const int width, const int level){
  float rand_heigth = float(rand() % heigth) + (rand() / (float) RAND_MAX);
  float rand_width  = float(rand() % width ) + (rand() / (float) RAND_MAX);
  float rand_level  = float(rand() % level ) + (rand() / (float) RAND_MAX);

  float dir_x = (rand() / (float) RAND_MAX);
  float dir_y = (rand() / (float) RAND_MAX);
  float dir_z = (rand() / (float) RAND_MAX);

  ray_cu r = {
    {rand_heigth, rand_width, rand_level, 1},
    {dir_x, dir_y, dir_z, 0}};
  return r;
}

std::vector<ray_cu> generate_rays(const int height, const int width, const int level, const unsigned max_rays){
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
