#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector_types.h>
#include <assert.h>
#include <vector>
#include <curand_kernel.h>
#include <cudachecks.h>
#include <importance_sampling.h>
#include <test_functions.h>
#include <calc_sample_phi_ase.h>
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>
#include <mesh.h>

#define SEED 1234

float calcDndtAse (unsigned &threads, 
		   unsigned &blocks, 
		   unsigned &hostRaysPerSample,
		   Mesh mesh,
		   Mesh hostMesh,
		   std::vector<double> *betaCellsVector,
		   float nTot,
		   float sigmaA,
		   float sigmaE,
		   float crystalFluorescence,
		   std::vector<double> *dndtAse){

  // Variable declaration
  // CPU
  double *hostImportance;
  unsigned *hostRaysPerPrism;
  cudaEvent_t start, stop;
  float runtimeGpu;
  unsigned *hostIndicesOfPrisms;
  float *hostPhiAse;

  // GPU
  float *phiAse;
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  double *importance;
  unsigned *indicesOfPrisms;
  float *sumPhi;
  unsigned *raysDump;
  unsigned *raysPerPrism;

  //OPTIMIZE: find perfect number of threads - MUST be the same as the size of shared memory in kernel
  threads = 256; 
  blocks = 200;

  hostPhiAse          = (float*)    malloc (hostMesh.numberOfSamples * sizeof(float));
  hostImportance      = (double*)   malloc (hostMesh.numberOfPrisms  * sizeof(double));
  hostRaysPerPrism    = (unsigned*) malloc (hostMesh.numberOfPrisms  * sizeof(unsigned));
  hostIndicesOfPrisms = (unsigned*) malloc (hostRaysPerSample        * sizeof(unsigned));
  runtimeGpu = 0.0;

  for(unsigned i=0; i < hostRaysPerSample; ++i) hostIndicesOfPrisms[i] = 0;
  for(unsigned i=0; i < hostMesh.numberOfSamples; ++i) hostPhiAse[i] = 0.f;
  for(unsigned i=0; i < hostMesh.numberOfPrisms; ++i) hostRaysPerPrism[i] = 1;
  for(unsigned i=0; i < hostMesh.numberOfPrisms; ++i) hostImportance[i] = 1.0;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // Init mersenne twister PRNG
  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, blocks * sizeof(curandStateMtgp32)));
  CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));
  CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));
  CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, blocks, SEED));

  // Memory allocation on device
  CUDA_CHECK_RETURN(cudaMalloc(&phiAse, hostMesh.numberOfSamples * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&importance, hostMesh.numberOfPrisms * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&indicesOfPrisms, hostRaysPerSample * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&raysPerPrism, hostMesh.numberOfPrisms * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&sumPhi, sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&raysDump, sizeof(unsigned)));

  // Copy host to device
  CUDA_CHECK_RETURN(cudaMemcpy(phiAse, hostPhiAse, hostMesh.numberOfSamples * sizeof(float), cudaMemcpyHostToDevice));

  // Calculate Phi Ase foreach sample
  fprintf(stderr, "\nC Start Phi Ase calculation\n");
  cudaEventRecord(start, 0);
  for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){
    if(sample_i % 200 == 0) fprintf(stderr, "C Sampling point %d/%d done\n", sample_i, hostMesh.numberOfSamples);

    importanceSampling(sample_i, mesh, hostRaysPerSample, sigmaA, sigmaE, nTot, importance, sumPhi, raysPerPrism, indicesOfPrisms, raysDump, threads, blocks);
    CUDA_CHECK_RETURN(cudaMemcpy(hostRaysPerPrism,raysPerPrism, hostMesh.numberOfPrisms*sizeof(unsigned),cudaMemcpyDeviceToHost));

    // Prism scheduling for gpu threads
    for(unsigned prism_i=0, absoluteRay = 0; prism_i < hostMesh.numberOfPrisms; ++prism_i){
      for(unsigned ray_i=0; ray_i < hostRaysPerPrism[prism_i]; ++ray_i){
        hostIndicesOfPrisms[absoluteRay++] = prism_i;
        assert(absoluteRay <= hostRaysPerSample);
      }
    }

    // Copy dynamic sample data to device
    CUDA_CHECK_RETURN(cudaMemcpy(indicesOfPrisms, hostIndicesOfPrisms, hostRaysPerSample * sizeof(unsigned), cudaMemcpyHostToDevice));

    // Start Kernel
    calcSamplePhiAse<<< blocks, threads >>>(devMTGPStates, mesh, indicesOfPrisms, importance, hostRaysPerSample, phiAse, sample_i, sigmaA, sigmaE, nTot);

  }
  // Copy solution back to host
  CUDA_CHECK_RETURN(cudaMemcpy(hostPhiAse, phiAse, hostMesh.numberOfSamples * sizeof(float), cudaMemcpyDeviceToHost));

  // Stop time
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&runtimeGpu, start, stop);

  // Calculate dndt Ase
  for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){
    hostPhiAse[sample_i] = float( (double(hostPhiAse[sample_i]) / (hostRaysPerSample * 4.0f * 3.14159)));
    double gain_local = double(nTot) * (betaCellsVector->at(sample_i)) * double(sigmaE + sigmaA) - double(nTot * sigmaA);
    dndtAse->at(sample_i) = gain_local * hostPhiAse[sample_i] / crystalFluorescence;
        
  }

  // Free Memory
  free(hostPhiAse);
  free(hostImportance);
  free(hostRaysPerPrism);
  free(hostIndicesOfPrisms);
  cudaFree(phiAse);
  cudaFree(importance);
  cudaFree(indicesOfPrisms);
  cudaFree(raysPerPrism);
  cudaFree(sumPhi);
  cudaFree(raysDump);
  cudaDeviceReset();

  return runtimeGpu;

}
