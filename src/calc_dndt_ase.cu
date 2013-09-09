#include "calc_dndt_ase.h"
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
#include "calc_sample_phi_ase.h"
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>

/* include MTGP pre-computed parameter sets */
/* include <curand_mtgp32dc_p_11213.h> */

#include <cuda_runtime_api.h>
#include <mesh.h>
#include <ctime> /* progressBar */
#include <progressbar.h> /*progressBar */

#define SEED 1234

/**
 * @brief Calculates which ray should start in which prism. Thus
 *        every thread in on gpu knows the index of the prism
 *        where its rays starts.
 *
 **/
void calcIndicesOfPrism(unsigned *indicesOfPrisms, unsigned *numberOfReflections, unsigned* raysPerPrism, unsigned reflectionSlices, unsigned raysPerSample, Mesh mesh, dim3 gridDim){
  for(unsigned reflection_i =0; reflection_i < reflectionSlices; ++reflection_i){
    for(unsigned prism_i=0, absoluteRay = 0; prism_i < mesh.numberOfPrisms; ++prism_i){
      unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;
      //unsigned wavelengthOffset = wave_i * mesh.numberOfPrisms;
      int reflectionPlane  = (reflection_i % 2 == 0)? -1 : 1;
      unsigned reflections = (reflection_i + 1) / 2;

      for(unsigned ray_i=0; ray_i < raysPerPrism[prism_i + reflectionOffset]; ++ray_i){
        indicesOfPrisms[absoluteRay] = prism_i;
        numberOfReflections[absoluteRay] = reflection_i;
        absoluteRay++;
        assert(absoluteRay <= raysPerSample);
      }

    }

  }
}


/**
 * @brief Gives every 200 blocks an index to the sigma_a/_e array or -1
 *        if this wavelength will be ignored.
 **/
void calcIndicesOfWavelengths(int *indicesOfWavelength, dim3 gridDim, std::vector<bool> ignoreWavelength){
  for(unsigned wave_i=0; wave_i < gridDim.y; ++wave_i){
    if(ignoreWavelength[wave_i]){
      indicesOfWavelength[wave_i] = -1;
    }
    else{
      indicesOfWavelength[wave_i] = wave_i;

    }

  }

}

double calcExpectation(double phiAse, double phiAseSquare, unsigned raysPerSample){
  double a = phiAseSquare / raysPerSample;
  double b = (phiAse / raysPerSample) * (phiAse / raysPerSample);

  return sqrt(abs((a - b) / raysPerSample));
}

float calcDndtAse (unsigned &threads, 
    unsigned &blocks,
    unsigned &hostRaysPerSample,
    Mesh mesh,
    Mesh hostMesh,
    std::vector<double> hostSigmaA,
    std::vector<double> hostSigmaE,
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
  int *hostNumberOfReflections;
  float *hostPhiAseSquare;
  time_t starttime,progressStartTime;
  unsigned hostRaysPerSampleSave;
  float expectationThreshold;
  unsigned maxRaysPerSample;
  unsigned maxReflections;
  unsigned reflectionSlices;
  unsigned *hostRaysDump;
  bool distributeRandomly;

  // GPU
  float *phiAse;
  float *phiAseSquare;
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  double *importance;
  unsigned *indicesOfPrisms;
  unsigned *numberOfReflections;
  unsigned *raysPerPrism;
  unsigned *cumulativeSums;

  // Variable Definitions
  dim3 blockDim(256);
  dim3 gridDim(200, hostSigmaE.size());
  threads = blockDim.x;
  blocks = gridDim.x;

  starttime = time(0);
  hostRaysPerSampleSave = hostRaysPerSample;
  expectationThreshold = 0.003;
  maxRaysPerSample = 100000000; // 100M
  maxReflections = 0;
  reflectionSlices = 1 + 2 * maxReflections;
  distributeRandomly = true;

  // Memory allocation on host
  hostPhiAseSquare         = (float*)    malloc (hostMesh.numberOfSamples * gridDim.y * sizeof(float));
  hostImportance           = (double*)   malloc (hostMesh.numberOfPrisms  * reflectionSlices * sizeof(double));
  hostRaysPerPrism         = (unsigned*) malloc (hostMesh.numberOfPrisms  * reflectionSlices * sizeof(unsigned));
  hostIndicesOfPrisms      = (unsigned*) malloc (maxRaysPerSample         * sizeof(unsigned));
  hostNumberOfReflections  = (int*)      malloc (maxRaysPerSample         * sizeof(int));
  hostRaysDump             = (unsigned*) malloc (1                        * sizeof(unsigned));

  for(unsigned i=0; i < maxRaysPerSample ; ++i) hostIndicesOfPrisms[i] = 0;
  for(unsigned i=0; i < maxRaysPerSample ; ++i) hostNumberOfReflections[i] = 0;
  for(unsigned i=0; i < hostMesh.numberOfSamples * gridDim.y; ++i) hostPhiAseSquare[i] = 0.f;
  for(unsigned i=0; i < hostMesh.numberOfPrisms * reflectionSlices; ++i) hostRaysPerPrism[i] = 1;
  for(unsigned i=0; i < hostMesh.numberOfPrisms * reflectionSlices; ++i) hostImportance[i] = 1.0;
  *hostRaysDump = 0;

  // CUDA Mersenne twister for more than 200 blocks (for every wavelength)
  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.x  * sizeof(curandStateMtgp32)));
  CUDA_CALL(cudaMalloc((void**)&devKernelParams,sizeof(mtgp32_kernel_params)));

  // TODO remove unused states (if using only 1 wavelength at a time...)
  for(unsigned wave_i = 0; wave_i < gridDim.y; ++wave_i){
    CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, &(devKernelParams[wave_i])));
    CURAND_CALL(curandMakeMTGP32KernelState(&(devMTGPStates[gridDim.x * wave_i]), mtgp32dc_params_fast_11213, &(devKernelParams[wave_i]), gridDim.x, SEED + wave_i));
  }

  // Memory allocation on device
  CUDA_CHECK_RETURN(cudaMalloc(&indicesOfPrisms, maxRaysPerSample * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&numberOfReflections, maxRaysPerSample * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&phiAse, hostMesh.numberOfSamples * gridDim.y * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&phiAseSquare, hostMesh.numberOfSamples * gridDim.y * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&importance, hostMesh.numberOfPrisms * reflectionSlices * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&raysPerPrism, hostMesh.numberOfPrisms * reflectionSlices * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&cumulativeSums,  hostMesh.numberOfPrisms * sizeof(unsigned)));

  // Copy host to device
  CUDA_CHECK_RETURN(cudaMemcpy(phiAse, &(hostPhiAse->at(0)), hostMesh.numberOfSamples * gridDim.y * sizeof(float), cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(phiAseSquare, hostPhiAseSquare, hostMesh.numberOfSamples * gridDim.y * sizeof(float), cudaMemcpyHostToDevice));

  // Calculate Phi Ase foreach sample
  fprintf(stderr, "\nC Start Phi Ase calculation\n");
  progressStartTime = time(0);
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

  for(unsigned wave_i = 0; wave_i < gridDim.y; ++wave_i){
    for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){
      int sampleOffset = sample_i + hostMesh.numberOfSamples * wave_i;
      hostRaysPerSample = hostRaysPerSampleSave;

      while(true){

        for(unsigned i=0; i < hostRaysPerSample && i < maxRaysPerSample ; ++i) hostIndicesOfPrisms[i] = 0;
        for(unsigned i=0; i < hostRaysPerSample && i < maxRaysPerSample ; ++i) hostNumberOfReflections[i] = 0;

        importanceSampling(sample_i, reflectionSlices, mesh, hostRaysPerSample, hostSigmaA[wave_i], hostSigmaE[wave_i], importance, raysPerPrism, hostRaysDump, distributeRandomly, blockDim, gridDim);
        CUDA_CHECK_RETURN(cudaMemcpy(hostRaysPerPrism, raysPerPrism, hostMesh.numberOfPrisms * reflectionSlices * sizeof(unsigned),cudaMemcpyDeviceToHost));

        // Prism scheduling for gpu threads
        calcIndicesOfPrism(hostIndicesOfPrisms, hostNumberOfReflections, hostRaysPerPrism, reflectionSlices, hostRaysPerSample, hostMesh, gridDim);
        CUDA_CHECK_RETURN(cudaMemcpy(indicesOfPrisms, hostIndicesOfPrisms, hostRaysPerSample * sizeof(unsigned), cudaMemcpyHostToDevice));
        CUDA_CHECK_RETURN(cudaMemcpy(numberOfReflections, hostNumberOfReflections, hostRaysPerSample * sizeof(unsigned), cudaMemcpyHostToDevice));

        // Start Kernel
        calcSamplePhiAse<<< 200, blockDim >>>(devMTGPStates, mesh, indicesOfPrisms, wave_i, numberOfReflections, importance, hostRaysPerSample, phiAse, phiAseSquare, sample_i, hostSigmaA[wave_i], hostSigmaE[wave_i]);


        // Copy solution (for this samplepoint) back to host
        CUDA_CHECK_RETURN(cudaMemcpy(&(hostPhiAse->at(sampleOffset)), &(phiAse[sampleOffset]), sizeof(float), cudaMemcpyDeviceToHost));
        CUDA_CHECK_RETURN(cudaMemcpy(&(hostPhiAseSquare[sampleOffset]), &(phiAseSquare[sampleOffset]), sizeof(float), cudaMemcpyDeviceToHost));

        // Check square error
        expectation->at(sampleOffset) =  calcExpectation(hostPhiAse->at(sampleOffset), hostPhiAseSquare[sampleOffset], hostRaysPerSample);

        if(expectation->at(sampleOffset) <= expectationThreshold) break;
        if((hostRaysPerSample * 10) > maxRaysPerSample)          break;

        // fprintf(stderr,"increasing from %d to %d\n",hostRaysPerSample, hostRaysPerSample*10);
        // If the threshold is still too high, increase the number of rays and reset the previously calculated value
        hostRaysPerSample *= 10;
        hostPhiAse->at(sampleOffset) = 0;
        hostPhiAseSquare[sampleOffset] = 0;
        CUDA_CHECK_RETURN( cudaMemcpy(&(phiAse[sampleOffset]), &(hostPhiAse->at(sampleOffset)), sizeof(float), cudaMemcpyHostToDevice));
        CUDA_CHECK_RETURN( cudaMemcpy(&(phiAseSquare[sampleOffset]), &(hostPhiAseSquare[sampleOffset]), sizeof(float), cudaMemcpyHostToDevice));

      }
      // Update progressbar
      if((sample_i+1) % 10 == 0) fancyProgressBar(sample_i,hostMesh.numberOfSamples,60,progressStartTime);

      // Calculate dndt Ase, after one point is completely sampled
      hostPhiAse->at(sampleOffset) = float((double(hostPhiAse->at(sampleOffset)) / (hostRaysPerSample * 4.0f * 3.14159)));
      double gain_local = double(hostMesh.nTot) * hostMesh.betaCells[sample_i] * double(hostSigmaE[wave_i] + hostSigmaA[wave_i]) - double(hostMesh.nTot * hostSigmaA[wave_i]);
      dndtAse->at(sampleOffset) = gain_local * hostPhiAse->at(sampleOffset) / hostMesh.crystalFluorescence;


    }
  }


  // Stop time
  runtime = difftime(time(0),starttime);

  // Free Memory
  // HINT Don't free importance if we return value to main
  free(hostImportance);
  free(hostRaysPerPrism);
  free(hostIndicesOfPrisms);
  free(hostRaysDump);
  cudaFree(phiAse);
  cudaFree(importance);
  cudaFree(indicesOfPrisms);
  cudaFree(raysPerPrism);
  cudaDeviceReset();

  return runtime;
}

