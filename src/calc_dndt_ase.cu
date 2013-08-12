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
		   std::vector<double> *hostSigmaA,
		   std::vector<double> *hostSigmaE,
		   std::vector<double> *dndtAse,
		   std::vector<float> *hostPhiAse,
		   std::vector<double> *expectation
		   ){

  // Variable declaration
  // CPU
  double *hostImportance;
  unsigned *hostRaysPerPrism;
  float runtime;
  unsigned *hostIndicesOfPrisms;
  float *hostPhiAseSquare;
  time_t starttime,progressStartTime;

  // GPU
  float *phiAse;
  float *phiAseSquare;
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  double *importance;
  unsigned *indicesOfPrisms;
  float *sumPhi;
  unsigned *raysDump;
  unsigned *raysPerPrism;
  unsigned *cumulativeSums;
  double * sigmaA;
  double * sigmaE;

  dim3 blockDim(256);
  dim3 gridDim(200, hostSigmaE->size());
  threads = blockDim.x;
  blocks = gridDim.x * gridDim.y;
    
  starttime = time(0);

  //hostPhiAse          = (float*)    malloc (hostMesh.numberOfSamples * gridDim.y * sizeof(float));
  hostPhiAseSquare    = (float*)    malloc (hostMesh.numberOfSamples * gridDim.y * sizeof(float));
  hostImportance      = (double*)   malloc (hostMesh.numberOfPrisms  * gridDim.y * sizeof(double));
  hostRaysPerPrism    = (unsigned*) malloc (hostMesh.numberOfPrisms  * gridDim.y * sizeof(unsigned));
  hostIndicesOfPrisms = (unsigned*) malloc (hostRaysPerSample        * gridDim.y * sizeof(unsigned));

  for(unsigned i=0; i < hostRaysPerSample * gridDim.y; ++i) hostIndicesOfPrisms[i] = 0;
  for(unsigned i=0; i < hostMesh.numberOfSamples * gridDim.y; ++i) hostPhiAseSquare[i] = 0.f;
  for(unsigned i=0; i < hostMesh.numberOfPrisms * gridDim.y; ++i) hostRaysPerPrism[i] = 1;
  for(unsigned i=0; i < hostMesh.numberOfPrisms * gridDim.y; ++i) hostImportance[i] = 1.0;

  // CUDA Mersenne twister for more than 200 blocks
  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.y * gridDim.x  * sizeof(curandStateMtgp32)));
  CUDA_CALL(cudaMalloc((void**)&devKernelParams, gridDim.y * sizeof(mtgp32_kernel_params)));
  for(unsigned mersenne_i = 0; mersenne_i < gridDim.y; ++mersenne_i){
    CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, &(devKernelParams[mersenne_i])));
    CURAND_CALL(curandMakeMTGP32KernelState(&(devMTGPStates[gridDim.x * mersenne_i]), mtgp32dc_params_fast_11213, &(devKernelParams[mersenne_i]), gridDim.x, SEED + mersenne_i));
  }

  // Memory allocation on device
  CUDA_CHECK_RETURN(cudaMalloc(&phiAse, hostMesh.numberOfSamples * gridDim.y * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&phiAseSquare, hostMesh.numberOfSamples * gridDim.y * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&importance, hostMesh.numberOfPrisms * gridDim.y * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&indicesOfPrisms, hostRaysPerSample * gridDim.y * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&raysPerPrism, hostMesh.numberOfPrisms * gridDim.y * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&sumPhi, gridDim.y * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&raysDump, gridDim.y * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&cumulativeSums,  hostMesh.numberOfPrisms * gridDim.y * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&sigmaA, gridDim.y * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&sigmaE, gridDim.y * sizeof(double)));

  // Copy host to device
  CUDA_CHECK_RETURN(cudaMemcpy(phiAse, &(hostPhiAse->at(0)), hostMesh.numberOfSamples * gridDim.y * sizeof(float), cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(phiAseSquare, hostPhiAseSquare, hostMesh.numberOfSamples * gridDim.y * sizeof(float), cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(sigmaA, &(hostSigmaA->at(0)), hostSigmaA->size() * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(sigmaE, &(hostSigmaE->at(0)), gridDim.y * sizeof(double), cudaMemcpyHostToDevice));
  
  // Calculate Phi Ase foreach sample
  fprintf(stderr, "\nC Start Phi Ase calculation\n");
  progressStartTime = time(0);
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){

    hostRaysPerSample = importanceSampling(sample_i, mesh, hostRaysPerSample, sigmaA, sigmaE, importance, sumPhi, raysPerPrism, indicesOfPrisms, raysDump, cumulativeSums, blockDim, gridDim);

    CUDA_CHECK_RETURN(cudaMemcpy(hostRaysPerPrism, raysPerPrism, hostMesh.numberOfPrisms * gridDim.y * sizeof(unsigned),cudaMemcpyDeviceToHost));

    // Prism scheduling for gpu threads
    for(unsigned wave_i=0; wave_i < gridDim.y; ++wave_i){
      for(unsigned prism_i=0, absoluteRay = 0; prism_i < hostMesh.numberOfPrisms; ++prism_i){
    	for(unsigned ray_i=0; ray_i < hostRaysPerPrism[prism_i + hostMesh.numberOfPrisms * wave_i]; ++ray_i){
    	  hostIndicesOfPrisms[absoluteRay + hostRaysPerSample * wave_i] = prism_i;
    	  absoluteRay++;
    	  assert(absoluteRay <= hostRaysPerSample);
    	}
      }
    }

    // Copy dynamic sample data to device
    CUDA_CHECK_RETURN(cudaMemcpy(indicesOfPrisms, hostIndicesOfPrisms, hostRaysPerSample * gridDim.y * sizeof(unsigned), cudaMemcpyHostToDevice));

    // Start Kernel
    calcSamplePhiAse<<< gridDim, blockDim , blockDim.x * sizeof(double)>>>(devMTGPStates, mesh, indicesOfPrisms, importance, hostRaysPerSample, phiAse, phiAseSquare, sample_i, sigmaA, sigmaE);

    // update progressbar
    if((sample_i+1) % 10 == 0) fancyProgressBar(sample_i,hostMesh.numberOfSamples,60,progressStartTime);

  }
  // Copy solution back to host
  CUDA_CHECK_RETURN(cudaMemcpy(&(hostPhiAse->at(0)), phiAse, hostMesh.numberOfSamples * gridDim.y * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CHECK_RETURN(cudaMemcpy(hostPhiAseSquare, phiAseSquare, hostMesh.numberOfSamples * gridDim.y * sizeof(float), cudaMemcpyDeviceToHost));

  // Calculate expectations
  for(unsigned wave_i = 0; wave_i < gridDim.y; ++wave_i){
    for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){
      int sampleOffset = sample_i + hostMesh.numberOfSamples * wave_i;
      double a = hostPhiAseSquare[sampleOffset] / hostRaysPerSample;
      double b = (hostPhiAse->at(sampleOffset) / hostRaysPerSample) * (hostPhiAse->at(sampleOffset) / hostRaysPerSample);

      expectation->at(sampleOffset) =  sqrt(abs((a - b) / hostRaysPerSample));
      
    }

  }

  // Calculate dndt Ase
  for(unsigned wave_i = 0; wave_i < gridDim.y; ++wave_i){
    for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){
      int sampleOffset = sample_i + hostMesh.numberOfSamples * wave_i;
      hostPhiAse->at(sampleOffset) = float((double(hostPhiAse->at(sampleOffset)) / (hostRaysPerSample * 4.0f * 3.14159)));
      double gain_local = double(hostMesh.nTot) * hostMesh.betaCells[sample_i] * double(hostSigmaE->at(wave_i) + hostSigmaA->at(wave_i)) - double(hostMesh.nTot * hostSigmaA->at(wave_i));
      dndtAse->at(sampleOffset) = gain_local * hostPhiAse->at(sampleOffset) / hostMesh.crystalFluorescence;

    }
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
  cudaFree(sigmaA);
  cudaFree(sigmaE);
  cudaDeviceReset();

  return runtime;

}
