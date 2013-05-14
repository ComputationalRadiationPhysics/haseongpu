#include "geometry_gpu.h"
#include "datatypes.h"
#include "curand_kernel.h"
#include "buildgrid.h"

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

__global__ void random_setup_kernel ( curandState * state, unsigned long seed )
{

  unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init ( seed, gid, 0, &state[gid] );
} 

/**
 * @brief Calculates ase phi for every sample point in samples[i].w
 * 
 * @param prisms            Array of prisms for intersection
 *        max_prisms        Total number of prisms
 *        samples           Array of sample points (.w defined as ase phi for this sample)
 *        blocks_per_sample How much blocks calculating on one smaple (a block has 256 threads)
 *        betas             Array of constant values for each prism
 *        globalState       Value for random number generation
 *
 * Every Thread generates a ray and tests this ray against every prism for
 * intersection. From the intersection distances of the rays will be
 * calculated the ase phi.
 *
 */
__global__ void ase_bruteforce_kernel(PrismCu* prisms, const unsigned max_prisms, PointCu *samples, const unsigned blocks_per_sample, const double *betas, curandState* globalState){
  // Cuda ids
  //unsigned tid = threadIdx.x;
  unsigned bid = blockIdx.x + blockIdx.y * gridDim.x;
  unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;

  // Local data
  const double sigma_e = 2.4E-20;
  const double sigma_a = 1.16E-21;
  const double N_tot   = 2.76E20;
  //const double clad_abs = 5.5;
  //const int clad_num = 3;
  unsigned prism_i;
  unsigned sample_i = bid / blocks_per_sample;
  double sumRayDistance = 0.f;
  double gain = 1;
  double distance = 0;
  double absRayDistance = 0;
  VectorCu rayDirection;

  // Random number generator
  curandState localState = globalState[gid];
  unsigned  raySourcePrism_i = selectPrism(gid, prisms, max_prisms);
  PointCu     sample = samples[sample_i];
  RayCu       ray = generateRayGpu(sample, prisms[raySourcePrism_i], localState);
  globalState[gid] = localState;

  // Precalculation of ray direction and total ray distance
  absRayDistance = distance_gpu(ray.P, ray.direction);
  rayDirection.x = ray.direction.x - ray.P.x;
  rayDirection.y = ray.direction.y - ray.P.y;
  rayDirection.z = ray.direction.z - ray.P.z;

  // Calculations
  __syncthreads();
  for(prism_i = 0; prism_i < max_prisms; ++prism_i){
    // In case of intersection, distance will be the intersection distance
    distance = collide_prism_gpu(prisms[prism_i], ray, rayDirection, absRayDistance);
    // Sum up distance to check the result
    sumRayDistance += distance;
    // Calculate the gain (has to be done for every cycle ?) 
    gain *= exp(N_tot * (betas[prism_i] * (sigma_e + sigma_a) - sigma_a) * distance);
    __syncthreads();
  }
  gain /= (absRayDistance * absRayDistance);
  __syncthreads();



  // Check Solution 
  if(fabs(sumRayDistance - absRayDistance) > 0.00001){
    printf("\033[31;1m[Error]\033[m Sample %d Ray %d with wrong distance real_distance(%f) != sum_distance(%f)\n", sample_i, gid, absRayDistance, sumRayDistance);
    return;
  }

  // Copy data to global
  // Cant explain multiplication by 0.007 but is done in original code
  // ==> has been removed, dont know why

  atomicAdd(&(samples[sample_i].w), (gain * betas[raySourcePrism_i]));

}

/**
 *  @brief Precalculatitions, preparing data for kernel execution and execute kernel that
 *         calculates ase phi of every sample point.
 *  
 *  @params samples    Sample points for ase phi
 *          prisms     3d prisms used for intersection
 *          betas      Constant values for each prism
 *          ase        Solutions of kernel calculations for every sample point
 *          rays_total Total number of rays
 *  
 *  @return runtime    Time the kernel spent for calculations
 *          threads    Number of threads per block (will be changed by reference)
 *          blocks     Number of total blocks (will be changed by reference)
 *          rays_total Can change because it has to fit onto hardware
 *
 * First the function calculates the number of blocks and threads for the kernel with
 * help of total number of rays and number of sample points. 
 * Problem data will be transfered from vectors to arrays and moved to the
 * device Memory.
 * The Kernel calculates the ase phi values for every sample point and
 * finally just some postcalculations.
 *
 */
float runAseBruteforceGpu(std::vector<PointCu> *samples, std::vector<PrismCu> *prisms, std::vector<double> *betas, std::vector<double> *ase, unsigned &threads, unsigned &blocks, unsigned &rays_total){
  PrismCu *h_prisms, *d_prisms;
  PointCu *h_samples, *d_samples;
  double *h_betas, *d_betas;
  float runtime_gpu = 0.0;
  cudaEvent_t start, stop;
  curandState *devStates;

  // Calculation of problem size (threads, blocks, blocks_per_sample)
  threads = 256;
  unsigned rays_per_sample = ceil(rays_total / (float)samples->size());
  const int blocks_per_sample = ceil(rays_per_sample / (float)threads);
  blocks = blocks_per_sample * samples->size();
  rays_total = blocks * threads;
  rays_per_sample = rays_total / samples->size();

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // Memory allocation on host
  fprintf(stderr, "C Prepare data for GPU-Calculations\n");
  CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_prisms, prisms->size() * sizeof(PrismCu), cudaHostAllocDefault));
  CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_samples, samples->size() * sizeof(PointCu), cudaHostAllocDefault));
  CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&h_betas, betas->size() * sizeof(double), cudaHostAllocDefault));

  // Memory initialisation on host
  unsigned prism_i;
  for(prism_i = 0; prism_i < prisms->size(); ++prism_i){
    h_prisms[prism_i] = prisms->at(prism_i);
  }

  unsigned sample_i;
  for(sample_i = 0; sample_i < samples->size(); ++sample_i){
    h_samples[sample_i] = samples->at(sample_i);
  }

  unsigned beta_i;
  for(beta_i = 0; beta_i < betas->size(); ++beta_i){
    h_betas[beta_i] = betas->at(beta_i);
  }

  /* Grid grid; */
  /* int4 dimGrid = {1, 1, 1, 0}; */
  /* prepareGrid(&grid, dimGrid, h_prisms, prisms->size()); */

  // Memory allocation on device
  CUDA_CHECK_RETURN(cudaMalloc(&d_prisms, prisms->size() * sizeof(PrismCu)));
  CUDA_CHECK_RETURN(cudaMalloc(&d_samples, samples->size() * sizeof(PointCu)));
  CUDA_CHECK_RETURN(cudaMalloc(&d_betas, betas->size() * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&devStates, threads * blocks * sizeof(curandState)));


  // Copy data from host to device
  CUDA_CHECK_RETURN(cudaMemcpy(d_prisms, h_prisms, prisms->size() * sizeof(PrismCu), cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(d_samples, h_samples, samples->size() * sizeof(PointCu), cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(d_betas, h_betas, betas->size() * sizeof(double), cudaMemcpyHostToDevice));
  cudaEventRecord(start, 0);
  
  // Start random setup kernel, is needed later to generate random numbers
  fprintf(stderr, "C Init GPU Random number\n");
  random_setup_kernel <<< blocks, threads >>> ( devStates, time(NULL) );

  // Start kernel
  fprintf(stderr, "C Start GPU Raytracing\n");
  ase_bruteforce_kernel<<< blocks, threads >>>(d_prisms, prisms->size(), d_samples, blocks_per_sample, d_betas, devStates);

  // Copy data from device to host
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&runtime_gpu, start, stop);

  CUDA_CHECK_RETURN(cudaMemcpy(h_samples, d_samples, samples->size() * sizeof(PointCu), cudaMemcpyDeviceToHost));

 
  // Copy data to vectors and scale by number of rays per sample
  for(sample_i = 0; sample_i < samples->size(); ++sample_i){
    ase->at(sample_i) = (h_samples[sample_i].w  / (rays_per_sample * 4.0f * 3.14159)) ; 

  }

  // Cleanup data  
  cudaFreeHost(h_prisms);
  cudaFreeHost(h_samples);
  cudaFreeHost(h_betas);
  cudaDeviceReset();
  
  return runtime_gpu;
}

#endif /* ASE_BRUTEFORCE_KERNEL_H */
