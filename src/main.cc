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
ray_cu   generate_ray(int height, int weight, float level);

//----------------------------------------------------
//  Calculations
//----------------------------------------------------
__global__ void trace_on_gpu(int max_triangles, triangle* triangles)
{
  //Thread ID
  //int thread_i = blockDim.y * blockIdx.y + threadIdx.y;
  //int block_j = blockDim.x * blockIdx.x + threadIdx.x;

}

int main(){
  const unsigned triangle_cnt = 2;
  const unsigned max_rays = 1000;
  const unsigned max_triangles = 20000;
  unsigned length = 20000 / (100 * 2);
  unsigned ray_cnt = 0;
  unsigned i;
  unsigned ray_i;
  ray_cu r;
  point_cu p;

  plane_cu pl = {
    {0, 0, 0, 1},
    {0, 0, 1, 0}};

  std::vector<triangle_cu> triangles = generate_triangles(length, length, 0);

  // CPU Raytracing
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
    r = generate_ray(length, length, 0);
    p = intersection(pl, r);
    for(i = 0; i < triangles.size(); ++i){
      if(collide(triangles[i], p)){
  	fprintf(stdout, "Ray %d on Triangle %d: Ahh collision, don't panic\n", ray_cnt, i);
  	// Do something on collision
      }

    }
    ray_cnt++;

  }

  // GPU Raytracing
  triangle_cu** h_triangles;
  triangle_cu** d_triangles;
  int threads = 0;
  int blocks = 0;

  // Memory allocation on host
  cudaHostAlloc( (void**)&h_triangles, max_triangles * sizeof(triangle_cu), cudaHostAllocDefault);

  // Memory allocation on device
  cudaMalloc(&d_triangles, max_triangles * sizeof(triangle_cu));

  // Copy data from host to device
  cudaMemcpy(d_triangles, h_triangles, max_triangles * sizeof(triangle_cu), cudaMemcpyHostToDevice);

  // Start kernel
  trace_on_gpu<<<threads, blocks>>>(d_triangles, max_triangles);

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
  //printf("%4.20f\n", (b.x + b.y + b.z));
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

ray_cu generate_ray(int height, int width, float level){
  float random1 = float(rand() % height) + (rand() / (float) RAND_MAX);
  float random2 = float(rand() % width) + (rand() / (float) RAND_MAX);
  ray_cu r = {
    {random1, random2, level, 1},
    {0,0,1, 0}};
  return r;
}

void print_point(point p){
  fprintf(stdout, "Point\n");
  fprintf(stdout, "x: %f\n", p.x);
  fprintf(stdout, "y: %f\n", p.y);
  fprintf(stdout, "z: %f\n", p.z);

}
