#include "calc_dndt_ase.h"
#include "map_rays_to_prisms.h"
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
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>


#define SEED 1234

/**
 * @brief Calculates which ray should start in which prism. Thus
 *        every thread in on gpu knows the index of the prism
 *        where its rays starts.
 *
 **/

double calcExpectation(double phiAse, double phiAseSquare, unsigned raysPerSample){
  double a = phiAseSquare / raysPerSample;
  double b = (phiAse / raysPerSample) * (phiAse / raysPerSample);

  return sqrt(abs((a - b) / raysPerSample));
}



double getDndtAse(const Mesh& mesh, const double sigmaA, const double sigmaE, const std::vector<float>& phiAse, const unsigned sampleOffset){
  unsigned sample_i = sampleOffset % mesh.numberOfSamples;
  double gain_local = mesh.nTot * mesh.betaCells[sample_i] * (sigmaE + sigmaA) - double(mesh.nTot * sigmaA);

  return gain_local * phiAse[sampleOffset] / mesh.crystalFluorescence;
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
  time_t starttime,progressStartTime;
  unsigned hostRaysPerSampleSave;
  unsigned maxReflections;
  unsigned reflectionSlices;
  bool distributeRandomly;

  // GPU
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;

  // Variable Definitions
  dim3 blockDim(256);
  dim3 gridDim(200, hostSigmaE.size());

  // give back to calling function
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


  // Memory allocation/init and copy for device
  thrust::device_vector<unsigned> numberOfReflections(maxRaysPerSample,0);
  thrust::device_vector<unsigned> indicesOfPrisms(maxRaysPerSample,0);
  thrust::device_vector<unsigned> raysPerPrism(hostMesh.numberOfPrisms * reflectionSlices, 1);
  thrust::device_vector<unsigned> prefixSum   (hostMesh.numberOfPrisms * reflectionSlices, 0);
  thrust::device_vector<double>   importance  (hostMesh.numberOfPrisms * reflectionSlices, 0);

  //OPTIMIZE: use only 1 value for phiAseSquare
  thrust::device_vector<float> phiAseSquare(hostMesh.numberOfSamples * gridDim.y, 0);
  thrust::device_vector<float> phiAse(hostPhiAse);

  // CUDA Mersenne twister for more than 200 blocks (for every wavelength)
  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.x  * sizeof(curandStateMtgp32)));
  CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));
  CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));
  CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED));

  // Calculate Phi Ase foreach sample
  fprintf(stderr, "\nC Start Phi Ase calculation\n");
  progressStartTime = time(0);
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

  // std::vector<unsigned> centerSample(expectation.size(), 0);

  for(unsigned wave_i = 0; wave_i < gridDim.y; ++wave_i){
    for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){
      int sampleOffset = sample_i + hostMesh.numberOfSamples * wave_i;
      hostRaysPerSample = hostRaysPerSampleSave;

      while(true){
        importanceSampling(
            sample_i, reflectionSlices, mesh, hostRaysPerSample, hostSigmaA[wave_i], hostSigmaE[wave_i],
            thrust::raw_pointer_cast(&importance[0]), thrust::raw_pointer_cast(&raysPerPrism[0]),
            distributeRandomly, blockDim, gridDim);

        // Prism scheduling for gpu threads
        mapRaysToPrisms(indicesOfPrisms,numberOfReflections,raysPerPrism,prefixSum,reflectionSlices,hostRaysPerSample,hostMesh.numberOfPrisms);

        // TESTING OUTPUT
        //if(sample_i == 1386)
        //  centerSample.assign(hostRaysPerPrism.begin(), hostRaysPerPrism.end());

        // Start Kernel
        calcSamplePhiAse<<< 200, blockDim >>>(
            devMTGPStates, 
            mesh, 
            thrust::raw_pointer_cast(&indicesOfPrisms[0]), 
            wave_i, 
            thrust::raw_pointer_cast(&numberOfReflections[0]), 
            thrust::raw_pointer_cast(&importance[0]),
            hostRaysPerSample, 
            thrust::raw_pointer_cast(&phiAse[0]), 
            thrust::raw_pointer_cast(&phiAseSquare[0]),
            sample_i, 
            hostSigmaA[wave_i], 
            hostSigmaE[wave_i]
            );

        // Check square error
        expectation.at(sampleOffset) = calcExpectation(phiAse[sampleOffset], phiAseSquare[sampleOffset], hostRaysPerSample);

        if(expectation[sampleOffset] < expectationThreshold) break;
        if(hostRaysPerSample * 10 > maxRaysPerSample)         break;

        // If the threshold is still too high, increase the number of rays and reset the previously calculated value
        hostRaysPerSample *= 10;
        phiAse[sampleOffset] = 0;
        phiAseSquare[sampleOffset] = 0;

      }
      // Update progressbar
      if((sample_i+1) % 10 == 0) fancyProgressBar(sample_i,hostMesh.numberOfSamples,60,progressStartTime);

      // get phiASE, normalize it and get dndtAse
      hostPhiAse.at(sampleOffset) = phiAse[sampleOffset];
      hostPhiAse[sampleOffset] /= hostRaysPerSample * 4.0f * 3.14159;
      dndtAse[sampleOffset] = getDndtAse(hostMesh,hostSigmaA[wave_i],hostSigmaE[wave_i],hostPhiAse,sampleOffset);

    }
  }

  // TESTING OUTPUT
  // expectation.assign(centerSample.begin(), centerSample.end());

  // Free Memory
  cudaFree(devMTGPStates);
  cudaFree(devKernelParams);

  return difftime(time(0),starttime);
}

