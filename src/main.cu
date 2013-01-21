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
#include "print.h"

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

__device__ PointCu intersectionRayTriangleGPU(PointCu rayOrigin, //Ursprung des Strahls
				   PointCu rayObjective, //Richtungsvektor des Strahls
				   PointCu p1, //1.Punkt des Dreiecks
				   PointCu p2, //2.Punkt des Dreiecks
				   PointCu p3) //3.Punkt des Dreiecks
{
  double s2; //2.barizentrische Koordinate des Dreiecks
  double s3; //3.barizentrische Koordinate des Dreiecks
  //1.barizentrische Koordinate des Dreiecks ergibt sich mit 1.-s2-s3
  double t; //Geradenparameter
  PointCu intersectionPoint = {0, 0, 0, 0};

  //Grenzwert fuer numerische Stabilitaet
  const double eps = 1e-6; //empirischer Wert, bei Moeller/Trumbore 1e-6

  //Variable fuer Determinante
  double determinante;

  //side12 und side13 sind Vektoren der Seiten
  //cross ist eine Hilfsvariable
  VectorCu side12, side13, rayDirection, cross;

  //Berechnung von Vektoren der Seiten:
  //1.Seite side12 von p1 nach p2
  side12.x = p2.x - p1.x;
  side12.y = p2.y - p1.y;
  side12.z = p2.z - p1.z;

  //2.Seite side13 von p1 nach p3
  side13.x = p3.x - p1.x;
  side13.y = p3.y - p1.y;
  side13.z = p3.z - p1.z;

  rayDirection.x = rayObjective.x - rayOrigin.x;
  rayDirection.y = rayObjective.y - rayOrigin.y;
  rayDirection.z = rayObjective.z - rayOrigin.z;

  //Gleichsetzen von Gereadengleichung und Ebenengleichung
  //Modell:  sp=p1+s2*(p2-p1)+s3*(p3-p1)=rayOrigin+t*rayDirection
  //Berechnung mit Cramerscher Regel zum Loesen des Gleichungssystems:
  //s2*(p2-p1)+s3*(p3-p1)-t*rayDirection = rayOrigin-p1 -> zu bestimmende Parameter: s2,s3,-t 

  //Kreuzprodukt von side13 und rayDirection
  cross.x = side13.y * rayDirection.z - side13.z * rayDirection.y;
  cross.y = side13.z * rayDirection.x - side13.x * rayDirection.z;
  cross.z = side13.x * rayDirection.y - side13.y * rayDirection.x;

  //Berechnung der Determinante mit Skalarprodukt
  determinante = cross.x * side12.x + cross.y * side12.y + cross.z * side12.z;


  //Test auf Parallelitaet
  //numerische Stabilitaet!!!

  if (determinante > -eps && determinante < eps){
    return intersectionPoint;
  }

  //Abstand Ursprung des Strahls zu p1
  VectorCu p1_rayOrigin; //=rayOrigin-p1;
  p1_rayOrigin.x = rayOrigin.x - p1.x;
  p1_rayOrigin.y = rayOrigin.y - p1.y;
  p1_rayOrigin.z = rayOrigin.z - p1.z;

  //barizentrische Koordinaten
  // sp=s1*p1+s2*p2+s3*p3
  //2. barizentrische Koordinate s2
  //=Skalarprodukt von p1_rayOrigin und cross
  s2 = cross.x * p1_rayOrigin.x + cross.y * p1_rayOrigin.y + cross.z * p1_rayOrigin.z;

  //Hilfsvariable
  VectorCu tempcross;
  //zunaenaehst Kreuzprodukt von rayDirection und side12
  tempcross.x = rayDirection.y * side12.z - rayDirection.z * side12.y;
  tempcross.y = rayDirection.z * side12.x - rayDirection.x * side12.z;
  tempcross.z = rayDirection.x * side12.y - rayDirection.y * side12.x;

  //s3=Skalarprodukt von rayDirection und side12
  //s3=(rayDirection x side12) *p1_rayOrigin
  s3 = tempcross.x * p1_rayOrigin.x + tempcross.y * p1_rayOrigin.y + tempcross.z * p1_rayOrigin.z;

  //Cramersche Regel -> Division durchfuehren
  double invdet = 1. / determinante;

  s2 = invdet*s2;
  s3 = invdet*s3;

  //weitere Verwendung der Hilfsvariable fuer Berechnung des Geradenparameters t
  //zunaechst Kreuzprodukt von side13 und side12
  tempcross.x = side13.y * side12.z - side13.z * side12.y;
  tempcross.y = side13.z * side12.x - side13.x * side12.z;
  tempcross.z = side13.x * side12.y - side13.y * side12.x;

  //t ist dann das Skalarprodukt von tempcross und p1_rayOrigin
  //t=(seite13,seite12) *p1_rayOrigin = -(seite12,seite13) *p1_rayOrigin
  t = tempcross.x * p1_rayOrigin.x + tempcross.y * p1_rayOrigin.y + tempcross.z * p1_rayOrigin.z;

  t = invdet*t;

  //Test,ob der Schnittpunkt innerhalb des Dreiecks liegt:

  //Ueberschereitungstest fuer barizentrische Koordinaten
  if (s2 < 0. || s2 > 1.) return intersectionPoint;

  //Ueberschereitungstest fuer barizentrische Koordinaten
  if (s3 < 0. || s3 > 1.) return intersectionPoint;

  //0 <= s1=1-s2-s3 <= 1 -> s2+s3<1   (s2+s3>0 schon durchgefuehrt,da s2>0 s3>0)
  if (s2 + s3 > 1.) return intersectionPoint;

  //Test, ob Strahl in Richtung des Dreiecks zeigt:
  if (t < 0.) return intersectionPoint;

  //Schnittpunktberechnung

  intersectionPoint.x = rayOrigin.x + t * rayDirection.x;
  intersectionPoint.y = rayOrigin.y + t * rayDirection.y;
  intersectionPoint.z = rayOrigin.z + t * rayDirection.z;
  intersectionPoint.w = 1;

  return intersectionPoint;
 
}

/**
   @brief Detects collisions of a triangle and a ray without
   a precondition.
**/
__device__ PointCu collide_triangle_gpu(TriangleCu t, RayCu r){
 return intersectionRayTriangleGPU(r.P, r.direction, t.A, t.B, t.C);

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
	  
	  if(distance_gpu(r.P, first_intersection) <= ray_distance && distance_gpu(r.P, intersection_point) > ray_distance)
	    return distance_gpu(r.direction, first_intersection);

	  if(distance_gpu(r.P, first_intersection) >= ray_distance && distance_gpu(r.P, intersection_point) < ray_distance)
	    return distance_gpu(r.direction, intersection_point);
	  
	  if(distance_gpu(r.P, first_intersection) > ray_distance || distance_gpu(r.direction, first_intersection) > ray_distance)
	    return 0;

	  return distance_gpu(first_intersection, intersection_point);
	}

      }

    }

  }

  return 0;
}

__global__ void trace_on_prisms(PrismCu* prisms, const unsigned max_prisms, RayCu* rays, const unsigned max_rays_per_sample, PointCu *samples, const unsigned blocks_per_sample){
  // Cuda ids
  //unsigned tid = threadIdx.x;
  unsigned bid = blockIdx.x + blockIdx.y * gridDim.x;
  unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;

  // Local data
  unsigned prism_i;
  unsigned sample_i = bid / blocks_per_sample;
  RayCu ray = rays[gid];
  unsigned beta_per_ray = 1;
  unsigned importance_per_prism = 1;
  
  // Calculations
  __syncthreads();
  for(prism_i = 0; prism_i < max_prisms; ++prism_i){
    float distance = fabs(collide_prism_gpu(prisms[prism_i], ray));

    ray.P.w += distance * beta_per_ray;
    __syncthreads();	
  }
  __syncthreads();

  // Check Solution
  if(fabs(ray.P.w - distance_gpu(ray.P, ray.direction)) > 0.00001){
	  printf("\033[31;1m[Error]\033[m Sample %d Ray %d with wrong distance real_distance(%f) != sum_distance(%f)\n", sample_i, gid, distance_gpu(ray.P, ray.direction), ray.P.w);
	  return;
  }

  // Copy data to global
  rays[gid].P.w = ray.P.w;
  atomicAdd(&(samples[sample_i].w), (ray.P.w * importance_per_prism));

}
//----------------------------------------------------
// Host Code
//----------------------------------------------------
int main(){
  const unsigned max_rays = 256;
  const unsigned max_triangles = 8;
  const unsigned length = ceil(sqrt(max_triangles / 2));
  const unsigned depth  = 4;
  unsigned ray_i, prism_i, sample_i;
  float runtime_gpu = 0.0;
  float runtime_cpu = 0.0;
  cudaEvent_t start, stop;
  bool use_cpu = true;
  bool use_gpu = true;

  // Generate testdata
  fprintf(stderr, "C Generate Testdata\n");
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
    fprintf(stderr, "C Start CPU-Raytracing\n");
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
    fprintf(stderr, "C Start GPU Raytracing\n");
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

