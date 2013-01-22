#include "geometry.h"
#include "datatypes.h"
#include <vector>

#ifndef ASE_BRUTEFORCE_CPU
#define ASE_BRUTEFORCE_CPU

float runAseBruteforceCpu(std::vector<PointCu> *samples, std::vector<PrismCu> *prisms, std::vector<RayCu> *rays, std::vector<float> *ase){
  unsigned ray_i, prism_i, sample_i;
  unsigned rays_per_sample = rays->size() / samples->size();
  float runtime_cpu = 0.0;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  // Start Calculations
  fprintf(stderr, "C Start CPU-Raytracing\n");
  for(sample_i = 0; sample_i < samples->size(); ++sample_i){
    for(ray_i = 0; ray_i < rays_per_sample; ++ray_i){
      for(prism_i = 0; prism_i < prisms->size(); ++prism_i){
	float distance = fabs(collide_prism(prisms->at(prism_i), rays->at(sample_i * rays_per_sample + ray_i)));
	if(distance > 0){
	  //fprintf(stdout, "CPU: Sample %d Ray %d hits on prism %d with distance %f\n", sample_i, ray_i, prism_i, distance);
	  ase->at(sample_i) += distance;
	  rays->at(sample_i * rays_per_sample + ray_i).P.w += distance;
	    
	}

      }
      // Evaluate solution by length test
      float d1 = distance(rays->at(sample_i * rays_per_sample + ray_i).P, rays->at(sample_i * rays_per_sample + ray_i).direction);
      float d2 = rays->at(sample_i * rays_per_sample + ray_i).P.w;
      if(fabs(d1-d2) > 0.000001){
      	fprintf(stderr, "\033[31;1m[Error]\033[m Sample %d Ray %d with wrong distance real_distance(%f) != sum_distance(%f)\n", sample_i, ray_i, d1, d2);
	  
      }
      
    }

  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&runtime_cpu, start, stop);
  return runtime_cpu;

}

#endif /* ASE_BRUTEFORCE_CPU */
