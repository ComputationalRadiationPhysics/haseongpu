#include "calc_dndt_ase.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <assert.h>
#include <vector>
#include <curand_kernel.h>
#include <cudachecks.h>
#include <importance_sampling.h>
#include <test_functions.h>
#include <cuda_utils.h> /* copyToDevice, copyFromDevice */
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
void calcIndicesOfPrism(std::vector<unsigned> &indicesOfPrisms, std::vector<unsigned> &numberOfReflections, std::vector<unsigned> raysPerPrism, unsigned reflectionSlices, unsigned raysPerSample, Mesh mesh){
  // Init vectors with zero (slow and not needed anymore)
  // for(unsigned i=0;  i < indicesOfPrisms.size() ; ++i) indicesOfPrisms[i] = 0;
  // for(unsigned i=0;  i < numberOfReflections.size() ; ++i) numberOfReflections[i] = 0;

  // Calc new values
  unsigned absoluteRay = 0;
  for(unsigned reflection_i =0; reflection_i < reflectionSlices; ++reflection_i){
    for(unsigned prism_i=0; prism_i < mesh.numberOfPrisms; ++prism_i){
      unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;
      for(unsigned ray_i=0; ray_i < raysPerPrism[prism_i + reflectionOffset]; ++ray_i){
        indicesOfPrisms[absoluteRay] = prism_i;
        numberOfReflections[absoluteRay] = reflection_i;
        absoluteRay++;
        assert(absoluteRay <= raysPerSample);

      }

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
		   unsigned maxRaysPerSample,
		   Mesh mesh,
		   Mesh hostMesh,
		   std::vector<double> hostSigmaA,
		   std::vector<double> hostSigmaE,
		   float expectationThreshold,
		   bool useReflections,
		   std::vector<double> &dndtAse,
		   std::vector<float> &hostPhiAse,
		   std::vector<double> &expectation
		   ){

  // Variable declaration
  // CPU
  float runtime;
  time_t starttime,progressStartTime;
  unsigned hostRaysPerSampleSave;
  unsigned maxReflections;
  unsigned reflectionSlices;
  bool distributeRandomly;

  std::cout << hostRaysPerSample << std::endl;
  std::cout << maxRaysPerSample << std::endl;

  // GPU
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;

  // Variable Definitions
  dim3 blockDim(256);
  dim3 gridDim(200, hostSigmaE.size());
  threads = blockDim.x;
  blocks = gridDim.x;

  starttime = time(0);
  hostRaysPerSampleSave = hostRaysPerSample;

  if(useReflections){
    maxReflections = hostMesh.getMaxReflections(); 
  }
  else {
    maxReflections = 0;
  }

  reflectionSlices = 1 + 2 * maxReflections;
  distributeRandomly = true;

  // Memory allocation on host
  std::vector<unsigned> hostIndicesOfPrisms(maxRaysPerSample, 0);
  std::vector<unsigned> hostNumberOfReflections(maxRaysPerSample, 0);
  std::vector<double>   hostImportance(hostMesh.numberOfPrisms * reflectionSlices, 0);
  std::vector<unsigned> hostRaysPerPrism(hostMesh.numberOfPrisms * reflectionSlices, 1);
  std::vector<float>    hostPhiAseSquare(hostMesh.numberOfSamples * gridDim.y, 0);

  // Memory allocation/init and copy for device
  unsigned *indicesOfPrisms     = copyToDevice(hostIndicesOfPrisms);
  unsigned *numberOfReflections = copyToDevice(hostNumberOfReflections);
  unsigned *raysPerPrism        = copyToDevice(hostRaysPerPrism);
  double   *importance          = copyToDevice(hostImportance);
  float    *phiAseSquare        = copyToDevice(hostPhiAseSquare);
  float    *phiAse              = copyToDevice(hostPhiAse);
  
  // CUDA Mersenne twister for more than 200 blocks (for every wavelength)
  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.x  * sizeof(curandStateMtgp32)));
  CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));

  // TODO remove unused states (if using only 1 wavelength at a time...)
  for(unsigned wave_i = 0; wave_i < gridDim.y; ++wave_i){
    CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, &(devKernelParams[wave_i])));
    CURAND_CALL(curandMakeMTGP32KernelState(&(devMTGPStates[gridDim.x * wave_i]), mtgp32dc_params_fast_11213, &(devKernelParams[wave_i]), gridDim.x, SEED + wave_i));
  }

  // Calculate Phi Ase foreach sample
  fprintf(stderr, "\nC Start Phi Ase calculation\n");
  progressStartTime = time(0);
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

  std::vector<unsigned> centerSample(expectation.size(), 0);

  for(unsigned wave_i = 0; wave_i < gridDim.y; ++wave_i){
    for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){
      int sampleOffset = sample_i + hostMesh.numberOfSamples * wave_i;
      hostRaysPerSample = hostRaysPerSampleSave;

      while(true){
        importanceSampling(sample_i, reflectionSlices, mesh, hostRaysPerSample, hostSigmaA[wave_i], hostSigmaE[wave_i], importance, raysPerPrism, distributeRandomly, blockDim, gridDim);
	copyFromDevice(hostRaysPerPrism, raysPerPrism);

        // Prism scheduling for gpu threads
        calcIndicesOfPrism(hostIndicesOfPrisms, hostNumberOfReflections, hostRaysPerPrism, reflectionSlices, hostRaysPerSample, hostMesh);
	copyToDevice(hostIndicesOfPrisms, indicesOfPrisms);
	copyToDevice(hostNumberOfReflections, numberOfReflections);

	// TESTING OUTPUT
	 if(sample_i == 1386)
	   centerSample.assign(hostRaysPerPrism.begin(), hostRaysPerPrism.end());

        // Start Kernel
        calcSamplePhiAse<<< 200, blockDim >>>(devMTGPStates, mesh, indicesOfPrisms, wave_i, numberOfReflections, importance, hostRaysPerSample, phiAse, phiAseSquare, sample_i, hostSigmaA[wave_i], hostSigmaE[wave_i]);

        // Copy solution (for this samplepoint) back to host
	hostPhiAse[sampleOffset]       = copyFromDevice(&(phiAse[sampleOffset]));
	hostPhiAseSquare[sampleOffset] = copyFromDevice(&(phiAseSquare[sampleOffset]));

        // Check square error
        expectation.at(sampleOffset) = calcExpectation(hostPhiAse.at(sampleOffset), hostPhiAseSquare[sampleOffset], hostRaysPerSample);

        if(expectation.at(sampleOffset) < expectationThreshold) break;
        if((hostRaysPerSample * 10) > maxRaysPerSample)         break;

        // fprintf(stderr,"increasing from %d to %d\n",hostRaysPerSample, hostRaysPerSample*10);
        // If the threshold is still too high, increase the number of rays and reset the previously calculated value
        hostRaysPerSample *= 10;
        hostPhiAse.at(sampleOffset) = 0;
        hostPhiAseSquare[sampleOffset] = 0;
	copyToDevice(hostPhiAse[sampleOffset], &(phiAse[sampleOffset]));
	copyToDevice(hostPhiAseSquare[sampleOffset], &(phiAseSquare[sampleOffset]));

      }
      // Update progressbar
      if((sample_i+1) % 10 == 0) fancyProgressBar(sample_i,hostMesh.numberOfSamples,60,progressStartTime);

      // Calculate dndt Ase, after one point is completely sampled
      hostPhiAse.at(sampleOffset) = float((double(hostPhiAse.at(sampleOffset)) / (hostRaysPerSample * 4.0f * 3.14159)));
      double gain_local = double(hostMesh.nTot) * hostMesh.betaCells[sample_i] * double(hostSigmaE[wave_i] + hostSigmaA[wave_i]) - double(hostMesh.nTot * hostSigmaA[wave_i]);
      dndtAse.at(sampleOffset) = gain_local * hostPhiAse.at(sampleOffset) / hostMesh.crystalFluorescence;


    }
  }

  // Stop time
  runtime = difftime(time(0),starttime);

  // TESTING OUTPUT
   expectation.assign(centerSample.begin(), centerSample.end());

  // Free Memory
  cudaFree(phiAse);
  cudaFree(importance);
  cudaFree(indicesOfPrisms);
  cudaFree(raysPerPrism);
  cudaFree(numberOfReflections);
  cudaFree(phiAseSquare);
  cudaDeviceReset();

  return runtime;
}

