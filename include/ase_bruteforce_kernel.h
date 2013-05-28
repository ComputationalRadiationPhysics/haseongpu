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

__global__ void random_setup_kernel ( curandState * state, unsigned long seed );


__global__ void ase_bruteforce_kernel(PrismCu* prisms, 
				      const unsigned max_prisms, 
				      PointCu *samples, 
				      const unsigned blocks_per_sample, 
				      const double *betas, 
				      const float cladAbsorption,
				      const float cladNumber,
				      const float nTot,
				      const float sigmaA,
				      const float sigmaE,
				      curandState* globalState);

float runAseBruteforceGpu(std::vector<PointCu> *samples, 
			  std::vector<PrismCu> *prisms, 
			  std::vector<double> *betas, 
			  std::vector<double> *ase,
			  float cladAbsorption,
			  float cladNumber,
			  float nTot,
			  float sigmaA,
			  float sigmaE,
			  unsigned &threads, 
			  unsigned &blocks, 
			  unsigned &rays_total);

#endif /* ASE_BRUTEFORCE_KERNEL_H */
