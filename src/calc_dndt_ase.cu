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
#include <ctime> /* progressBar */
#include <progressbar.h> /*progressBar */

#define SEED 4321

float calcDndtAse (unsigned &threads, 
		   unsigned &blocks,
		   unsigned &hostRaysPerSample,
		   Mesh mesh,
		   Mesh hostMesh,
		   std::vector<double> *betaCellsVector,
		   float nTot,
		   std::vector<float> *hostSigmaA,
		   std::vector<float> *hostSigmaE,
		   float crystalFluorescence,
		   std::vector<double> *dndtAse
		   ){

  // Variable declaration
  // CPU
  double *hostImportance;
  unsigned *hostRaysPerPrism;
  float runtime;
  unsigned *hostIndicesOfPrisms;
  float *hostPhiAse;
  time_t starttime,progressStartTime;

  // GPU
  float *phiAse;
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  double *importance;
  unsigned *indicesOfPrisms;
  float *sumPhi;
  unsigned *raysDump;
  unsigned *raysPerPrism;
  unsigned *cumulativeSums;
  float * sigmaA;
  float * sigmaE;

  //OPTIMIZE: find perfect number of threads - MUST be the same as the size of shared memory in kernel
  dim3 blockDim(256);
  dim3 gridDim(200, hostSigmaE->size());
  threads = blockDim.x;
  blocks = gridDim.x * gridDim.y;
    
  starttime = time(0);

  hostPhiAse          = (float*)    malloc (hostMesh.numberOfSamples * hostSigmaE->size() * sizeof(float));
  hostImportance      = (double*)   malloc (hostMesh.numberOfPrisms  * hostSigmaE->size() * sizeof(double));
  hostRaysPerPrism    = (unsigned*) malloc (hostMesh.numberOfPrisms  * hostSigmaE->size() * sizeof(unsigned));
  hostIndicesOfPrisms = (unsigned*) malloc (hostRaysPerSample        * hostSigmaE->size() * sizeof(unsigned));

  for(unsigned i=0; i < hostRaysPerSample * hostSigmaE->size(); ++i) hostIndicesOfPrisms[i] = 0;
  for(unsigned i=0; i < hostMesh.numberOfSamples * hostSigmaE->size(); ++i) hostPhiAse[i] = 0.f;
  for(unsigned i=0; i < hostMesh.numberOfPrisms * hostSigmaE->size(); ++i) hostRaysPerPrism[i] = 1;
  for(unsigned i=0; i < hostMesh.numberOfPrisms * hostSigmaE->size(); ++i) hostImportance[i] = 1.0;

  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.x * gridDim.y * sizeof(curandStateMtgp32)));

  // TODO RUN for each hostSigmaE->(i) ???
  CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));
  CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));
  CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED));

  // Memory allocation on device
  CUDA_CHECK_RETURN(cudaMalloc(&phiAse, hostMesh.numberOfSamples * hostSigmaE->size() * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&importance, hostMesh.numberOfPrisms * hostSigmaE->size() * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&indicesOfPrisms, hostRaysPerSample * hostSigmaE->size() * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&raysPerPrism, hostMesh.numberOfPrisms * hostSigmaE->size() * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&sumPhi, hostSigmaE->size() * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&raysDump, hostSigmaE->size() * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&cumulativeSums,  hostMesh.numberOfPrisms * hostSigmaE->size() * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&sigmaA, hostSigmaE->size() * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&sigmaE, hostSigmaE->size() * sizeof(float)));

  // Copy host to device
  CUDA_CHECK_RETURN(cudaMemcpy(phiAse, hostPhiAse, hostMesh.numberOfSamples * hostSigmaE->size() * sizeof(float), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(sigmaA, &(hostSigmaA->at(0)), hostSigmaA->size() * sizeof(float), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(sigmaE, &(hostSigmaE->at(0)), hostSigmaE->size() * sizeof(float), cudaMemcpyHostToDevice));
  

  // Calculate Phi Ase foreach sample
  fprintf(stderr, "\nC Start Phi Ase calculation\n");
  progressStartTime = time(0);
  for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){

    importanceSampling(sample_i, mesh, hostRaysPerSample, sigmaA, sigmaE, nTot, importance, sumPhi, raysPerPrism, indicesOfPrisms, raysDump, cumulativeSums, blockDim, gridDim);

    CUDA_CHECK_RETURN(cudaMemcpy(hostRaysPerPrism, raysPerPrism, hostMesh.numberOfPrisms * sizeof(unsigned),cudaMemcpyDeviceToHost));

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
    calcSamplePhiAse<<< gridDim, blockDim >>>(devMTGPStates, mesh, indicesOfPrisms, importance, hostRaysPerSample, phiAse, sample_i, hostSigmaA->at(0), hostSigmaE->at(0), nTot);

    // update progressbar
    if((sample_i+1) % 10 == 0) fancyProgressBar(sample_i,hostMesh.numberOfSamples,60,progressStartTime);

  }
  // Copy solution back to host
  CUDA_CHECK_RETURN(cudaMemcpy(hostPhiAse, phiAse, hostMesh.numberOfSamples * sizeof(float), cudaMemcpyDeviceToHost));

  // Calculate dndt Ase
  for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){
    hostPhiAse[sample_i] = float( (double(hostPhiAse[sample_i]) / (hostRaysPerSample * 4.0f * 3.14159)));
    double gain_local = double(nTot) * (betaCellsVector->at(sample_i)) * double(hostSigmaE->at(0) + hostSigmaA->at(0)) - double(nTot * hostSigmaA->at(0));
    dndtAse->at(sample_i) = gain_local * hostPhiAse[sample_i] / crystalFluorescence;

  }
  // Stop time
  runtime = difftime(time(0),starttime);


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

  return runtime;

}
