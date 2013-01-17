#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include <vector>
#include "curand_kernel.h"

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
typedef float4 point_cu;
typedef float4 vector_cu;

typedef struct triangle_cu{
  point_cu A;
  point_cu B;
  point_cu C;
} TRIANGLE_CU;

typedef struct prism_cu{
  triangle_cu t1; // height in w coordinate
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
// Device Code
//----------------------------------------------------
__device__ float distance_gpu(point_cu a, point_cu b){
  float d = sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y),2) + pow((b.z - a.z),2));
  return fabs(d);
}
__device__ vector_cu crossproduct_gpu(vector_cu a, vector_cu b){
  vector_cu c = {
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
  };
  return c;
}

__device__ float skalar_mul_gpu(vector_cu a, vector_cu b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

__device__ float4 to_barycentric_gpu(triangle_cu tr, ray_cu ray){
  float4 b = {0,0,0,0};
  vector_cu e1, e2, q, s, r;
  point_cu p0, p1, p2;
  float a, f, u, v, t;

  p0 = tr.A;
  p1 = tr.B;
  p2 = tr.C;

  e1.x = p1.x - p0.x;
  e1.y = p1.y - p0.y;
  e1.z = p1.z - p0.z;

  e2.x = p2.x - p0.x;
  e2.y = p2.y - p0.y;
  e2.z = p2.z - p0.z;

  q = crossproduct_gpu(ray.direction, e2);
  a = skalar_mul_gpu(e1, q);
  
  // a is to close to 0
  if(fabs(a) < 0.0000001)
    return b;

  f = 1 / a;
  
  s.x = ray.P.x - p0.x;
  s.y = ray.P.y - p0.y;
  s.z = ray.P.z - p0.z;

  u = f * skalar_mul_gpu(s, q);

  if(u < 0.0)
    return b;

  r = crossproduct_gpu(s, e1);
  v = f * skalar_mul_gpu(ray.direction, r);
  if( v < 0.0 || (u + v) > 1)
    return b;
  
  t = f * skalar_mul_gpu(e2, q);
  
  b.x = u;
  b.y = v;
  b.z = t;
  b.w = 1;

  return b;
}

__device__ point_cu intersection_gpu(plane_cu pl, ray_cu r){
  point_cu intersection_point = {0, 0, 0, 0};

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
  // this case for parallel rays, will be ignored for easier calculations
  float tmp = (n1*p1 + n2*p2 + n3*p3);
  if(tmp == 0)
    return intersection_point;

  d = n1*a1 + n2*a2 + n3*a3;
  t = (d - n1*x1 - n2*x2 - n3*x3) / tmp;

  intersection_point.x = x1 + t * p1;
  intersection_point.y = x2 + t * p2;
  intersection_point.z = x3 + t * p3;
  intersection_point.w = 1;
  return intersection_point;

}

/**
   @brief Detects collisions of a triangle and a ray without
   a precondition.
**/
__device__ point_cu collide_triangle_gpu(triangle_cu t, ray_cu r){
  plane_cu pl;
  float b1, b2, b3, c1, c2, c3;

  b1 = t.B.x - t.A.x;
  b2 = t.B.y - t.A.y;
  b3 = t.B.z - t.A.z;

  c1 = t.C.x - t.A.x;
  c2 = t.C.y - t.A.y;
  c3 = t.C.z - t.A.z;

  pl.P = t.A;
  pl.normal.x = (b2*c3 - b3*c2);
  pl.normal.y = (b3*c1 - b1*c3);
  pl.normal.z = (b1*c2 - b2*c1);

  float4 b = to_barycentric_gpu(t, r);
  // Maybe we can calculate intersection be barycentric coords
  point_cu p = intersection_gpu(pl, r);
  // Can be reduced to return p;
  if(b.w == 1){
    return p;
  }
  else{
    point_cu no_inter = {0,0,0,0}; 
    return no_inter;
  }



}

__device__ float collide_prism_gpu(prism_cu pr, ray_cu r){
  //bool has_collide;
  point_cu intersections[2];
  point_cu A1 = pr.t1.A;
  point_cu B1 = pr.t1.B;
  point_cu C1 = pr.t1.C;
  point_cu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.z + pr.t1.A.w, 1};
  point_cu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.z + pr.t1.B.w, 1};
  point_cu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.z + pr.t1.C.w, 1};

  triangle_cu triangles[8] = {
    pr.t1,
    {A2, B2, C2},
    {A1, B1, A2},
    {B1, B2, A2},
    {B1, C1, C2},
    {B1, B2, C2},
    {A1, C1, C2},
    {A1, A2, C2}};

  unsigned i; 
  unsigned j = 0;
  // test for collision on all triangles of an prism
  for(i = 0; i < 8; ++i){
    point_cu p = collide_triangle_gpu(triangles[i], r);
    if(p.w == 0)
    // No Collision for this triangle
      continue;
    // Filter double Collision on edges or vertices
    if(j != 0){
      if(intersections[j-1].x != p.x || intersections[j-1].y != p.y || intersections[j-1].z != p.z){
	intersections[j++] = p;
      }
    }
    else{
      intersections[j++] = p;
    }

  }
  if(j > 1)
    return distance_gpu(intersections[0], intersections[1]);
  else
    return 0;
}

__global__ void trace_on_prisms(prism_cu* prisms, const unsigned max_prisms, ray_cu* rays, const unsigned max_rays_per_sample, point_cu *samples, const unsigned blocks_per_sample){
  // Cuda ids
  unsigned tid = threadIdx.x;
  unsigned bid = blockIdx.x + blockIdx.y * gridDim.x;
  unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;

  // Local data
  unsigned prism_i;
  unsigned sample_i = bid / blocks_per_sample;
  point_cu sample_point = samples[sample_i];
  ray_cu ray = rays[gid];
  
  __syncthreads();
  for(prism_i = 0; prism_i < max_prisms; ++prism_i){
    float distance = fabs(collide_prism_gpu(prisms[prism_i], ray));
    if(distance > 0){
      ray.P.w += distance;

    }
	
  }
  rays[gid].P.w = ray.P.w;
  atomicAdd(&(samples[sample_i].w), ray.P.w);
  __syncthreads();


}

//----------------------------------------------------
// Auxillary function declaration
//----------------------------------------------------
// Geometrie functions
point_cu  collide_triangle(triangle_cu t, ray_cu r);
float  collide_prism(prism_cu pr, ray_cu r);
float4 to_barycentric(triangle_cu t, ray_cu r);
point_cu intersection(plane_cu p, ray_cu r);
float distance(point_cu a, point_cu b);
vector_cu crossproduct(vector_cu a, vector_cu b);
float skalar_mul(vector_cu a, vector_cu b);

// Testdata generation
std::vector<triangle_cu> generate_triangles(int height, int width, float level);
std::vector<prism_cu> generate_prisms(int height, int width, float level);
std::vector<ray_cu> generate_rays(int height, int width, int level, unsigned max_rays);
std::vector<point_cu> generate_samples(int height, int width, int level);
std::vector<ray_cu> generate_sample_rays(int height, int width, int level, unsigned max_rays, point_cu sample);
ray_cu   generate_ray(int height, int weight, int level);

// Debug functions
void print_point(point_cu p);
void print_vector(vector_cu v);
void print_plane(plane_cu pl);

//----------------------------------------------------
// Host Code
//----------------------------------------------------
int main(){
  const unsigned max_rays = 256;
  const unsigned max_triangles = 8;
  const unsigned length = ceil(sqrt(max_triangles / 2));
  const unsigned depth  = 1;
  unsigned ray_i, prism_i, sample_i;
  float runtime_gpu = 0.0;
  float runtime_cpu = 0.0;
  cudaEvent_t start, stop;
  bool use_cpu = true;
  bool use_gpu = true;

  // Generate testdata
  std::vector<prism_cu> prisms = generate_prisms(length, length, depth);
  std::vector<point_cu> samples = generate_samples(length, length, depth);
  std::vector<ray_cu> rays;
  std::vector<float> sample_data(samples.size(), 0);
  std::vector<float> ray_data(samples.size() * max_rays, 0);
  for(sample_i = 0; sample_i < samples.size(); ++sample_i){
    std::vector<ray_cu> sample_ray = generate_sample_rays(length, length, depth, max_rays, samples[sample_i]);
    for(ray_i = 0; ray_i < sample_ray.size(); ++ray_i)
      rays.push_back(sample_ray[ray_i]);
  }
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // CPU Raytracing
  cudaEventRecord(start, 0);
  if(use_cpu){
    for(sample_i = 0; sample_i < samples.size(); ++sample_i){
      for(ray_i = 0; ray_i < max_rays; ++ray_i){
	for(prism_i = 0; prism_i < prisms.size(); ++prism_i){
	  float distance = fabs(collide_prism(prisms[prism_i], rays[sample_i * max_rays + ray_i]));
	  if(distance > 0){
	    fprintf(stdout, "CPU: Sample %d Ray %d hits on prism %d with distance %f\n", sample_i, ray_i, prism_i, distance);
	    sample_data[sample_i] += distance;
	    ray_data[(sample_i * max_rays) + ray_i] += distance;
	    
	  }

	}

      }

    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&runtime_cpu, start, stop);
  }

  // GPU Raytracing
  ray_cu *h_rays, *d_rays;
  prism_cu *h_prisms, *d_prisms;
  point_cu *h_samples, *d_samples;

  int threads = 256;
  int blocks_per_sample = ceil(max_rays / threads);
  int blocks = blocks_per_sample * samples.size();
  
  if(use_gpu){

    // Memory allocation on host
    CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_prisms, prisms.size() * sizeof(prism_cu), cudaHostAllocDefault));
    CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_rays, samples.size() * max_rays * sizeof(ray_cu), cudaHostAllocDefault));
    CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_samples, samples.size() * sizeof(point_cu), cudaHostAllocDefault));

    // Memory initialisation on host
    for(ray_i = 0; ray_i < max_rays * samples.size(); ++ray_i){
      h_rays[ray_i] = rays[ray_i];
    }

    for(prism_i = 0; prism_i < prisms.size(); ++prism_i){
      h_prisms[prism_i] = prisms[prism_i];
    }

    for(sample_i = 0; sample_i < samples.size(); ++sample_i){
      h_samples[sample_i] = samples[sample_i];
    }

    // Memory allocation on device
    CUDA_CHECK_RETURN(cudaMalloc(&d_rays, samples.size() * max_rays * sizeof(ray_cu)));
    CUDA_CHECK_RETURN(cudaMalloc(&d_prisms, prisms.size() * sizeof(prism_cu)));
    CUDA_CHECK_RETURN(cudaMalloc(&d_samples, samples.size() * sizeof(point_cu)));

    // Copy data from host to device
    cudaEventRecord(start, 0);
    CUDA_CHECK_RETURN(cudaMemcpy(d_rays, h_rays, samples.size() * max_rays * sizeof(ray_cu), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(d_prisms, h_prisms, prisms.size() * sizeof(prism_cu), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(d_samples, h_samples, samples.size() * sizeof(point_cu), cudaMemcpyHostToDevice));

    // Start kernel
    //trace_on_prisms<<<threads, blocks>>>(d_prisms, max_prisms, d_rays, max_rays, d_samples, blocks_per_sample);
    trace_on_prisms<<<blocks, threads>>>(d_prisms, prisms.size(), d_rays, max_rays, d_samples, blocks_per_sample);

    // Copy data from device to host
    CUDA_CHECK_RETURN(cudaMemcpy(h_samples, d_samples, samples.size() * sizeof(point_cu), cudaMemcpyDeviceToHost));
    CUDA_CHECK_RETURN(cudaMemcpy(h_rays, d_rays, samples.size() * max_rays * sizeof(ray_cu), cudaMemcpyDeviceToHost));
    
    // Evaluate device data
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&runtime_gpu, start, stop);
    }    
  if(use_gpu && use_cpu){
    for(sample_i = 0; sample_i < samples.size(); ++sample_i){
      if(fabs(sample_data[sample_i] - h_samples[sample_i].w) < 0.1 && use_cpu && use_gpu){
	fprintf(stderr, "CPU == GPU: Sample %d with value %f \n", sample_i, sample_data[sample_i]);

      }
      else{
	fprintf(stderr, "\033[31;1m[Error]\033[m CPU(%f) != GPU(%f) on sample %d\n", sample_data[sample_i], h_samples[sample_i].w, sample_i);

      }

    }
    


    for(ray_i = 0; ray_i < rays.size(); ++ray_i){
      if(fabs(ray_data[ray_i] - h_rays[ray_i].P.w) < 0.00001){
	//fprintf(stderr, "CPU == GPU: Ray %d with value %f \n", ray_i, ray_data[ray_i]);

      }
      else{
	fprintf(stderr, "\033[31;1m[Error]\033[m CPU(%f) != GPU(%f) on ray %d\n", ray_data[ray_i], h_rays[ray_i].P.w, ray_i);

      }
    }
  }

  // Print statistics
  fprintf(stderr, "\n");
  fprintf(stderr, "Prism             : %d\n", (int)prisms.size());
  fprintf(stderr, "Triangles         : %d\n", (int)prisms.size() * 8);
  fprintf(stderr, "Samples           : %d\n", (int)samples.size());
  fprintf(stderr, "Rays/Sample       : %d\n", max_rays);
  fprintf(stderr, "GPU Blocks        : %d\n", blocks);
  fprintf(stderr, "GPU Threads       : %d\n", threads);
  fprintf(stderr, "GPU Blocks/Sample : %d\n", blocks_per_sample);
  fprintf(stderr, "Runtime_GPU       : %f s\n", runtime_gpu / 1000.0);
  fprintf(stderr, "Runtime_CPU       : %f s\n", runtime_cpu / 1000.0);
  fprintf(stderr, "Speedup CPU/GPU   : %.1f\n", runtime_cpu / runtime_gpu);
  fprintf(stderr, "\n");

  // Cleanup
  cudaFreeHost(h_rays);
  cudaFreeHost(h_prisms);
  cudaFreeHost(h_samples);

  return 0;
}

//----------------------------------------------------
// Auxillary function definition
//----------------------------------------------------

/**
   @brief Algorithm based on a paper. exact copy of pseudo code.
 **/
float4 to_barycentric(triangle_cu tr, ray_cu ray){
  float4 b = {0,0,0,0};
  vector_cu e1, e2, q, s, r;
  point_cu p0, p1, p2;
  float a, f, u, v, t;

  p0 = tr.A;
  p1 = tr.B;
  p2 = tr.C;

  e1.x = p1.x - p0.x;
  e1.y = p1.y - p0.y;
  e1.z = p1.z - p0.z;

  e2.x = p2.x - p0.x;
  e2.y = p2.y - p0.y;
  e2.z = p2.z - p0.z;

  q = crossproduct(ray.direction, e2);
  a = skalar_mul(e1, q);
  
  // a is to close to 0
  if(fabs(a) < 0.000001)
    return b;

  f = 1 / a;
  
  s.x = ray.P.x - p0.x;
  s.y = ray.P.y - p0.y;
  s.z = ray.P.z - p0.z;

  u = f * skalar_mul(s, q);

  if(u < 0.0)
    return b;

  r = crossproduct(s, e1);
  v = f * skalar_mul(ray.direction, r);
  if( v < 0.0 || (u + v) > 1)
    return b;
  
  t = f * skalar_mul(e2, q);
  
  b.x = u;
  b.y = v;
  b.z = t;
  b.w = 1;

  return b;
}

/**
   @brief Detects collisions of a triangle and a ray without
   a precondition.
**/
point_cu collide_triangle(triangle_cu t, ray_cu r){
  plane_cu pl;
  float b1, b2, b3, c1, c2, c3;

  b1 = t.B.x - t.A.x;
  b2 = t.B.y - t.A.y;
  b3 = t.B.z - t.A.z;

  c1 = t.C.x - t.A.x;
  c2 = t.C.y - t.A.y;
  c3 = t.C.z - t.A.z;

  pl.P = t.A;
  pl.normal.x = (b2*c3 - b3*c2);
  pl.normal.y = (b3*c1 - b1*c3);
  pl.normal.z = (b1*c2 - b2*c1);

  float4 b = to_barycentric(t, r);
  // Maybe we can calculate intersection be barycentric coords
  point_cu p = intersection(pl, r);
  // Can be reduced to return p;
  if(b.w == 1){
    return p;
  }
  else{
    point_cu no_inter = {0,0,0,0}; 
    return no_inter;
  }

}

float collide_prism(prism_cu pr, ray_cu r){
  //bool has_collide;
  point_cu intersections[2];
  point_cu A1 = pr.t1.A;
  point_cu B1 = pr.t1.B;
  point_cu C1 = pr.t1.C;
  point_cu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.z + pr.t1.A.w, 1};
  point_cu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.z + pr.t1.B.w, 1};
  point_cu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.z + pr.t1.C.w, 1};

  triangle_cu triangles[8] = {
    pr.t1,
    {A2, B2, C2},
    {A1, B1, A2},
    {B1, B2, A2},
    {B1, C1, C2},
    {B1, B2, C2},
    {A1, C1, C2},
    {A1, A2, C2}};

  unsigned i; 
  unsigned j = 0;
  // test for collision on all triangles of an prism
  for(i = 0; i < 8; ++i){
    point_cu p = collide_triangle(triangles[i], r);
    if(p.w == 0)
    // No Collision for this triangle
      continue;
    // Filter double Collision on edges or vertices
    if(j != 0){
      if(intersections[j-1].x != p.x || intersections[j-1].y != p.y || intersections[j-1].z != p.z){
	intersections[j++] = p;
      }
    }
    else{
      intersections[j++] = p;
    }

  }
  if(j > 1)
    return distance(intersections[0], intersections[1]);
  else
    return 0;
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
  point_cu intersection_point = {0, 0, 0, 0};

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
  // this case for parallel rays, will be ignored for easier calculations
  float tmp = (n1*p1 + n2*p2 + n3*p3);
  if(tmp == 0)
    return intersection_point;

  d = n1*a1 + n2*a2 + n3*a3;
  t = (d - n1*x1 - n2*x2 - n3*x3) / tmp;

  intersection_point.x = x1 + t * p1;
  intersection_point.y = x2 + t * p2;
  intersection_point.z = x3 + t * p3;
  intersection_point.w = 1;
  return intersection_point;

}

float distance(point_cu a, point_cu b){
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

std::vector<prism_cu> generate_prisms(int height, int width, float level){
  int h,w,l;
  std::vector<prism_cu> prisms;
  for(l = 0; l < level; ++l){
    for(h = 0; h < height; ++h){
      for(w = 0; w < width; ++w){
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
    {rand_heigth, rand_width, rand_level, 0},
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

std::vector<point_cu> generate_samples(int height, int width, int level){
  std::vector<point_cu> sample_points;
  int h,w,l;
  for(l = 0; l <= level; ++l){
    for(h = 0; h <= height; ++h){
      for(w = 0; w <= width; ++w){
	point_cu p = {float(h), float(w), float(l), 0};
	sample_points.push_back(p);
      }
    }
  }
  return sample_points;
}

std::vector<ray_cu> generate_sample_rays(int height, int width, int level, unsigned max_rays, point_cu sample){
  std::vector<ray_cu> rays;
  unsigned ray_i;
  for(ray_i = 0; ray_i < max_rays; ++ray_i){
    ray_cu ray = generate_ray(height, width, level);
    ray.direction.x  =  sample.x - ray.P.x;
    ray.direction.y  =  sample.y - ray.P.y;
    ray.direction.z  =  sample.z - ray.P.z;
    ray.P = sample;
    rays.push_back(ray);
  }
  return rays;
}

vector_cu crossproduct(vector_cu a, vector_cu b){
  vector_cu c = {
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
  };
  return c;
}

float skalar_mul(vector_cu a, vector_cu b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

void print_point(point_cu p){
  fprintf(stderr, "Point (%f, %f, %f)\n",p.x, p.y, p.z);

}

void print_vector(vector_cu v){
  fprintf(stderr, "Vector (%f, %f, %f)\n",v.x, v.y, v.z);

}

void print_plane(plane_cu pl){
  fprintf(stderr, "Plane: \n\t");
  print_point(pl.P);
  fprintf(stderr, "\t");
  print_vector(pl.normal);
}
