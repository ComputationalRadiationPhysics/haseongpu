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

double calcExpectation(const double phiAse, const double phiAseSquare, const unsigned raysPerSample){
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
    const unsigned maxRaysPerSample,
    const Mesh& mesh,
    const Mesh& hostMesh,
    const std::vector<double>& hostSigmaA,
    const std::vector<double>& hostSigmaE,
    const float expectationThreshold,
    const bool useReflections,
    std::vector<double> &dndtAse,
    std::vector<float> &hostPhiAse,
    std::vector<double> &expectation
    ){

  // Optimization to use more L1 cache
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

  using thrust::device_vector;
  using  thrust::raw_pointer_cast;

  // variable Definitions CPU
  time_t starttime                = time(0);
  unsigned hostRaysPerSampleSave  = hostRaysPerSample;
  unsigned maxReflections         = useReflections ? hostMesh.getMaxReflections() : 0;
  unsigned reflectionSlices       = 1 + (2 * maxReflections);
  unsigned numberOfWavelengths    = hostSigmaE.size();
  bool distributeRandomly         = true;
  dim3 blockDim(256);             //OPTIMIZE: find perfect number of threads
  dim3 gridDim(200);              //can't be more than 200 due to restrictions from the Mersenne Twister
  threads                         = blockDim.x; // give back to calling function
  blocks                          = gridDim.x;
  // Memory allocation/init and copy for device memory
  device_vector<unsigned> numberOfReflections(maxRaysPerSample,  0);
  device_vector<unsigned> indicesOfPrisms    (maxRaysPerSample,  0);
  device_vector<float>    phiAse             (hostPhiAse.size(), 0);
  device_vector<unsigned> raysPerPrism       (hostMesh.numberOfPrisms * reflectionSlices, 1);
  device_vector<unsigned> prefixSum          (hostMesh.numberOfPrisms * reflectionSlices, 0);
  device_vector<double>   importance         (hostMesh.numberOfPrisms * reflectionSlices, 0);
  device_vector<float>    phiAseSquare       (hostMesh.numberOfSamples * numberOfWavelengths, 0); //OPTIMIZE: use only 1 value
 
 // CUDA Mersenne twister (can not have more than 200 blocks!)
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.x  * sizeof(curandStateMtgp32)));
  CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));
  CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));
  CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED));
  

  // Calculate Phi Ase for each wavelength
  for(unsigned wave_i = 0; wave_i < numberOfWavelengths; ++wave_i){
    time_t progressStartTime = time(0);
    fprintf(stderr, "\nC Phi_ASE calculation for wavelength %d\n",wave_i);

    //calculation for each sample point
    for(unsigned sample_i = 0; sample_i < hostMesh.numberOfSamples; ++sample_i){
      unsigned sampleOffset  = sample_i + hostMesh.numberOfSamples * wave_i;
      hostRaysPerSample = hostRaysPerSampleSave;

      while(true){
        importanceSampling(
            sample_i, reflectionSlices, mesh, hostRaysPerSample, hostSigmaA[wave_i], hostSigmaE[wave_i],
            raw_pointer_cast(&importance[0]), raw_pointer_cast(&raysPerPrism[0]),
            distributeRandomly, blockDim, gridDim
            );

        // Prism scheduling for gpu threads
        mapRaysToPrisms(indicesOfPrisms,numberOfReflections,raysPerPrism,prefixSum,reflectionSlices,hostRaysPerSample,hostMesh.numberOfPrisms);

        // Start Kernel
        calcSamplePhiAse<<< gridDim, blockDim >>>(
            devMTGPStates,
            mesh, 
            raw_pointer_cast(&indicesOfPrisms[0]), 
            wave_i, 
            raw_pointer_cast(&numberOfReflections[0]), 
            raw_pointer_cast(&importance[0]),
            hostRaysPerSample, 
            raw_pointer_cast(&phiAse[0]), 
            raw_pointer_cast(&phiAseSquare[0]),
            sample_i, 
            hostSigmaA[wave_i], hostSigmaE[wave_i]
            );

        // Check square error
        expectation.at(sampleOffset) = calcExpectation(phiAse[sampleOffset], phiAseSquare[sampleOffset], hostRaysPerSample);

        if(expectation[sampleOffset] < expectationThreshold) break;
        if(hostRaysPerSample * 10 > maxRaysPerSample)        break;

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

  // Free Memory
  cudaFree(devMTGPStates);
  cudaFree(devKernelParams);

  return difftime(time(0),starttime);
}
