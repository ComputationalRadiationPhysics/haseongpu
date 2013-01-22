#include "geometry_gpu.h"
#include "datatypes.h"

#ifndef ASE_BRUTEFORCE_KERNEL_H
#define ASE_BRUTEFORCE_KERNEL_H

#define CUDA_CHECK_RETURN(value) {				\
  cudaError_t _m_cudaStat = value;				\
  if (_m_cudaStat != cudaSuccess) {				\
    fprintf(stderr, "Error %s at line %d in file %s\n",	\
    cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);	\
    exit(1);							\
  }								\
}


__global__ void ase_bruteforce_kernel(PrismCu* prisms, const unsigned max_prisms, RayCu* rays, const unsigned max_rays_per_sample, PointCu *samples, const unsigned blocks_per_sample){
  // Cuda ids
  unsigned tid = threadIdx.x;
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

float runAseBruteforceGpu(std::vector<PointCu> *samples, std::vector<PrismCu> *prisms, std::vector<RayCu> *rays, std::vector<float> *ase, unsigned threads){
  RayCu *h_rays, *d_rays;
  PrismCu *h_prisms, *d_prisms;
  PointCu *h_samples, *d_samples;
  float runtime_gpu = 0.0;
  cudaEvent_t start, stop;

  int blocks_per_sample = ceil(rays->size() / (threads * samples->size()));
  int blocks = blocks_per_sample * samples->size();
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  fprintf(stderr, "C Start GPU Raytracing\n");
  // Memory allocation on host
  CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_prisms, prisms->size() * sizeof(PrismCu), cudaHostAllocDefault));
  CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_rays, rays->size() * sizeof(RayCu), cudaHostAllocDefault));
  CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_samples, samples->size() * sizeof(PointCu), cudaHostAllocDefault));

  // Memory initialisation on host
  unsigned ray_i;
  for(ray_i = 0; ray_i < rays->size(); ++ray_i){
    h_rays[ray_i] = rays->at(ray_i);
  }

  unsigned prism_i;
  for(prism_i = 0; prism_i < prisms->size(); ++prism_i){
    h_prisms[prism_i] = prisms->at(prism_i);
  }

  unsigned sample_i;
  for(sample_i = 0; sample_i < samples->size(); ++sample_i){
    h_samples[sample_i] = samples->at(sample_i);
  }

  // Memory allocation on device
  CUDA_CHECK_RETURN(cudaMalloc(&d_rays, rays->size() * sizeof(RayCu)));
  CUDA_CHECK_RETURN(cudaMalloc(&d_prisms, prisms->size() * sizeof(PrismCu)));
  CUDA_CHECK_RETURN(cudaMalloc(&d_samples, samples->size() * sizeof(PointCu)));

  // Copy data from host to device
  CUDA_CHECK_RETURN(cudaMemcpy(d_rays, h_rays, rays->size() * sizeof(RayCu), cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(d_prisms, h_prisms, prisms->size() * sizeof(PrismCu), cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(d_samples, h_samples, samples->size() * sizeof(PointCu), cudaMemcpyHostToDevice));
  cudaEventRecord(start, 0);
  
  // Start kernel
  ase_bruteforce_kernel<<<blocks, threads>>>(d_prisms, prisms->size(), d_rays, rays->size() / samples->size(), d_samples, blocks_per_sample);

  // Copy data from device to host
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&runtime_gpu, start, stop);

  CUDA_CHECK_RETURN(cudaMemcpy(h_samples, d_samples, samples->size() * sizeof(PointCu), cudaMemcpyDeviceToHost));
  CUDA_CHECK_RETURN(cudaMemcpy(h_rays, d_rays, rays->size() * sizeof(RayCu), cudaMemcpyDeviceToHost));

 
  // Copy data to vectors
  for(sample_i = 0; sample_i < samples->size(); ++sample_i){
    ase->at(sample_i) = h_samples[sample_i].w;

  }

  // Cleanup data  
  cudaFreeHost(h_rays);
  cudaFreeHost(h_prisms);
  cudaFreeHost(h_samples);

  
  return runtime_gpu;
}

#endif /* ASE_BRUTEFORCE_KERNEL_H */
