// Libraies
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include <vector>

// User header files
#include "geometrie.h"
#include "datatypes.h"
#include "generate_testdata.h"

#define CUDA_CHECK_RETURN(value) {				\
  cudaError_t _m_cudaStat = value;				\
  if (_m_cudaStat != cudaSuccess) {				\
    fprintf(stderr, "Error %s at line %d in file %s\n",	\
    cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);	\
    exit(1);							\
  }								\
}

//----------------------------------------------------
// Device Code
//----------------------------------------------------

__device__ float distance_gpu(PointCu a, PointCu b){
  float d = sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y),2) + pow((b.z - a.z),2));
  return fabs(d);
}
__device__ VectorCu crossproduct_gpu(VectorCu a, VectorCu b){
  VectorCu c = {
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
  };
  return c;
}


__device__ float skalar_mul_gpu(VectorCu a, VectorCu b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

/**
   @brief Calculates the barycentric coordinates of the triangle
          and the intersectionpoint of the ray and the triangle.

   @return PointCu {0,0,0,0} if there is no intersection triangle/ray
   @return PointCu {x,y,z,1} barycentric coordinates of intersection triangle/ray
 **/
__device__ float4 to_barycentric_gpu(TriangleCu tr, RayCu ray){
  float4 b = {0,0,0,0};
  VectorCu e1, e2, q, s, r, ray_direction;
  PointCu p0, p1, p2;
  float a, f, u, v, t;

  p0 = tr.A;
  p1 = tr.B;
  p2 = tr.C;

  ray_direction.x = ray.direction.x - ray.P.x;
  ray_direction.y = ray.direction.y - ray.P.y;
  ray_direction.z = ray.direction.z - ray.P.z;

  e1.x = p1.x - p0.x;
  e1.y = p1.y - p0.y;
  e1.z = p1.z - p0.z;

  e2.x = p2.x - p0.x;
  e2.y = p2.y - p0.y;
  e2.z = p2.z - p0.z;

  q = crossproduct_gpu(ray_direction, e2);
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
  v = f * skalar_mul_gpu(ray_direction, r);
  if( v < 0.0 || (u + v) > 1)
    return b;
  
  t = f * skalar_mul_gpu(e2, q);
  
  b.x = u;
  b.y = v;
  b.z = t;
  b.w = 1;

  return b;
}

__device__ PointCu intersection_gpu(PlaneCu pl, RayCu r){
  PointCu intersection_point = {0, 0, 0, 0};

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

  p1 = r.direction.x - r.P.x;
  p2 = r.direction.y - r.P.y;
  p3 = r.direction.z - r.P.z;

  // calculation of intersection
  // this case for parallel rays, will be ignored for easier calculations
  float denominator = (n1*p1 + n2*p2 + n3*p3);
  if(abs(denominator) <= 0.000001)
    return intersection_point;

  d = n1*a1 + n2*a2 + n3*a3;
  t = (d - n1*x1 - n2*x2 - n3*x3) / denominator;

  // ignore intersections before the ray 
  if(t < 0)
    return intersection_point;

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
__device__ PointCu collide_triangle_gpu(TriangleCu t, RayCu r){
  PlaneCu pl;
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
  PointCu p = intersection_gpu(pl, r);
  if(b.w == 1 && p.w == 1){
    return p;
  }
  else{
    PointCu no_inter = {0,0,0,0}; 
    return no_inter;
  }

}

/* __device__ float collide_prism_gpu(PrismCu pr, RayCu r){ */
/*   //bool has_collide; */
/*   PointCu intersections[2]; */
/*   PointCu A1 = pr.t1.A; */
/*   PointCu B1 = pr.t1.B; */
/*   PointCu C1 = pr.t1.C; */
/*   PointCu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.z + pr.t1.A.w, 1}; */
/*   PointCu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.z + pr.t1.B.w, 1}; */
/*   PointCu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.z + pr.t1.C.w, 1}; */

/*   TriangleCu triangles[8] = { */
/*     pr.t1, */
/*     {A2, B2, C2}, */
/*     {A1, B1, A2}, */
/*     {B1, B2, A2}, */
/*     {B1, C1, C2}, */
/*     {B1, B2, C2}, */
/*     {A1, C1, C2}, */
/*     {A1, A2, C2}}; */

/*   unsigned i;  */
/*   unsigned j = 0; */
/*   // test for collision on all triangles of an prism */
/*   for(i = 0, j = 0; i < 8; ++i){ */
/*     PointCu p = collide_triangle_gpu(triangles[i], r); */
/*     if(p.w == 0) */
/*     // No Collision for this triangle */
/*       continue; */
/*     // Filter double Collision on edges or vertices */
/*     if(j != 0){ */
/*       if(intersections[j-1].x != p.x || intersections[j-1].y != p.y || intersections[j-1].z != p.z){ */
/* 	intersections[j++] = p; */
/*       } */
/*     } */
/*     else{ */
/*       intersections[j++] = p; */
/*     } */

/*   } */
/*   //if(j > 1) */
/*     return (j/2) * distance_gpu(intersections[0], intersections[1]); */
/*     //else */
/*     //return 0; */
/* } */


/**
   @attention slower than commented collide_prism_gpu but more pretty
**/
__device__ float collide_prism_gpu(PrismCu pr, RayCu r){
  //bool has_collide;
  PointCu first_intersection = {0, 0, 0, 0};
  PointCu intersection_point = {0, 0, 0, 0};
  PointCu A1 = pr.t1.A;
  PointCu B1 = pr.t1.B;
  PointCu C1 = pr.t1.C;
  PointCu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.z + pr.t1.A.w, 1};
  PointCu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.z + pr.t1.B.w, 1};
  PointCu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.z + pr.t1.C.w, 1};
  float ray_distance = distance_gpu(r.P, r.direction);
  TriangleCu triangles[8] = {
    pr.t1,
    {A2, B2, C2},
    {A1, B1, A2},
    {B1, B2, A2},
    {B1, C1, C2},
    {B1, B2, C2},
    {A1, C1, C2},
    {A1, A2, C2}};

  // test for collision on all triangles of an prism
  unsigned i; 
  for(i = 0; i < 8; ++i){
    intersection_point = collide_triangle_gpu(triangles[i], r);
    if(intersection_point.w != 0){
      if(first_intersection.w == 0){
	first_intersection = intersection_point;
      }
      else{
	// Filter double collisions
	if(first_intersection.x != intersection_point.x || first_intersection.y != intersection_point.y || first_intersection.z != intersection_point.z){
	  /*
	  if(distance_gpu(r.P, first_intersection) <= ray_distance && distance_gpu(r.P, intersection_point) > ray_distance)
	    return distance_gpu(r.direction, first_intersection);

	  if(distance_gpu(r.P, first_intersection) >= ray_distance && distance_gpu(r.P, intersection_point) < ray_distance)
	    return distance_gpu(r.direction, intersection_point);
	  
	  if(distance_gpu(r.P, first_intersection) > ray_distance || distance_gpu(r.direction, first_intersection) > ray_distance)
	    return 0;
	  */
	  return distance_gpu(first_intersection, intersection_point);
	}

      }

    }

  }

  return 0;
}

__global__ void trace_on_prisms(PrismCu* prisms, const unsigned max_prisms, RayCu* rays, const unsigned max_rays_per_sample, PointCu *samples, const unsigned blocks_per_sample){
  // Cuda ids
  unsigned tid = threadIdx.x;
  unsigned bid = blockIdx.x + blockIdx.y * gridDim.x;
  unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;

  // Local data
  unsigned prism_i;
  unsigned sample_i = bid / blocks_per_sample;
  PointCu sample_point = samples[sample_i];
  RayCu ray = rays[gid];
  
  __syncthreads();
  for(prism_i = 0; prism_i < max_prisms; ++prism_i){
    float distance = fabs(collide_prism_gpu(prisms[prism_i], ray));
    if(distance > 0){
      ray.P.w += distance;
    }
    __syncthreads();	
  }
  rays[gid].P.w = ray.P.w;
  atomicAdd(&(samples[sample_i].w), ray.P.w);



}

//----------------------------------------------------
// Auxillary function declaration
//----------------------------------------------------
// Debug functions
void print_point(PointCu p);
void print_vector(VectorCu v);
void print_plane(PlaneCu pl);

//----------------------------------------------------
// Host Code
//----------------------------------------------------
int main(){
  const unsigned max_rays = 256;
  const unsigned max_triangles = 32;
  const unsigned length = ceil(sqrt(max_triangles / 2));
  const unsigned depth  = 2;
  unsigned ray_i, prism_i, sample_i;
  float runtime_gpu = 0.0;
  float runtime_cpu = 0.0;
  cudaEvent_t start, stop;
  bool use_cpu = true;
  bool use_gpu = true;

  // Generate testdata
  std::vector<PrismCu> prisms = generate_prisms(length, length, depth);
  std::vector<PointCu> samples = generate_samples(length, length, depth);
  std::vector<RayCu> rays;
  std::vector<float> sample_data(samples.size(), 0);
  std::vector<float> ray_data(samples.size() * max_rays, 0);
  for(sample_i = 0; sample_i < samples.size(); ++sample_i){
    std::vector<RayCu> sample_ray = generate_sample_rays(length, length, depth, max_rays, samples[sample_i]);
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
	float d1 = distance(rays[sample_i * max_rays + ray_i].P, rays[sample_i * max_rays + ray_i].direction);
	float d2 = ray_data[(sample_i * max_rays) + ray_i];
	if(fabs(d1-d2) > 0.000001){
	  fprintf(stderr, "\033[31;1m[Error]\033[m Sample %d Ray %d with wrong distance real_distance(%f) != sum_distance(%f)\n", sample_i, ray_i, d1, d2);
	  
	  fprintf(stderr, "Sample: ");print_point(rays[sample_i * max_rays + ray_i].P);
	  fprintf(stderr, "Objective: ");print_point(rays[sample_i * max_rays + ray_i].direction);
	}
	//assert(fabs(d1 - d2) < 0.000001);
      }

    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&runtime_cpu, start, stop);
  }

  // GPU Raytracing
  RayCu *h_rays, *d_rays;
  PrismCu *h_prisms, *d_prisms;
  PointCu *h_samples, *d_samples;

  int threads = 256;
  int blocks_per_sample = ceil(max_rays / threads);
  int blocks = blocks_per_sample * samples.size();
  
  if(use_gpu){

    // Memory allocation on host
    CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_prisms, prisms.size() * sizeof(PrismCu), cudaHostAllocDefault));
    CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_rays, samples.size() * max_rays * sizeof(RayCu), cudaHostAllocDefault));
    CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_samples, samples.size() * sizeof(PointCu), cudaHostAllocDefault));

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
    CUDA_CHECK_RETURN(cudaMalloc(&d_rays, samples.size() * max_rays * sizeof(RayCu)));
    CUDA_CHECK_RETURN(cudaMalloc(&d_prisms, prisms.size() * sizeof(PrismCu)));
    CUDA_CHECK_RETURN(cudaMalloc(&d_samples, samples.size() * sizeof(PointCu)));

    // Copy data from host to device
    cudaEventRecord(start, 0);
    CUDA_CHECK_RETURN(cudaMemcpy(d_rays, h_rays, samples.size() * max_rays * sizeof(RayCu), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(d_prisms, h_prisms, prisms.size() * sizeof(PrismCu), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(d_samples, h_samples, samples.size() * sizeof(PointCu), cudaMemcpyHostToDevice));

    // Start kernel
    trace_on_prisms<<<blocks, threads>>>(d_prisms, prisms.size(), d_rays, max_rays, d_samples, blocks_per_sample);

    // Copy data from device to host
    CUDA_CHECK_RETURN(cudaMemcpy(h_samples, d_samples, samples.size() * sizeof(PointCu), cudaMemcpyDeviceToHost));
    CUDA_CHECK_RETURN(cudaMemcpy(h_rays, d_rays, samples.size() * max_rays * sizeof(RayCu), cudaMemcpyDeviceToHost));
    
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
void print_point(PointCu p){
  fprintf(stderr, "Point (%f, %f, %f)\n",p.x, p.y, p.z);

}

void print_vector(VectorCu v){
  fprintf(stderr, "Vector (%f, %f, %f)\n",v.x, v.y, v.z);

}

void print_plane(PlaneCu pl){
  fprintf(stderr, "Plane: \n\t");
  print_point(pl.P);
  fprintf(stderr, "\t");
  print_vector(pl.normal);
}
