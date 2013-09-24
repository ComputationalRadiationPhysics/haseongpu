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
    unsigned &hRaysPerSample,
    const unsigned maxRaysPerSample,
    const Mesh& dMesh,
    const Mesh& hMesh,
    const std::vector<double>& hSigmaA,
    const std::vector<double>& hSigmaE,
    const float expectationThreshold,
    const bool useReflections,
    std::vector<double> &dndtAse,
    std::vector<float> &hPhiAse,
    std::vector<double> &expectation
    ){

  // Optimization to use more L1 cache
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

  using thrust::device_vector;
  using  thrust::raw_pointer_cast;

  // variable Definitions CPU
  time_t starttime                = time(0);
  unsigned hRaysPerSampleSave     = hRaysPerSample;
  unsigned maxReflections         = useReflections ? hMesh.getMaxReflections() : 0;
  unsigned reflectionSlices       = 1 + (2 * maxReflections);
  unsigned numberOfWavelengths    = hSigmaE.size();
  bool distributeRandomly         = true;
  dim3 blockDim(128);             
  dim3 gridDim(200);              //can't be more than 200 due to restrictions from the Mersenne Twister
  threads                         = blockDim.x; // give back to calling function
  blocks                          = gridDim.x;

  // Memory allocation/init and copy for device memory
  device_vector<unsigned> dNumberOfReflections(maxRaysPerSample,  0);
  device_vector<unsigned> dIndicesOfPrisms    (maxRaysPerSample,  0);
  device_vector<float>    dPhiAse             (hPhiAse.size(), 0);
  device_vector<unsigned> dRaysPerPrism       (hMesh.numberOfPrisms * reflectionSlices, 1);
  device_vector<unsigned> dPrefixSum          (hMesh.numberOfPrisms * reflectionSlices, 0);
  device_vector<double>   dImportance         (hMesh.numberOfPrisms * reflectionSlices, 0);
  device_vector<float>    dPhiAseSquare       (hMesh.numberOfSamples * numberOfWavelengths, 0); //OPTIMIZE: use only 1 value
 
  // CUDA Mersenne twister (can not have more than 200 blocks!)
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.x  * sizeof(curandStateMtgp32)));
  CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));
  CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));
  CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED));
  

  thrust::host_vector<unsigned> hNumberOfReflections(maxRaysPerSample,0);
  thrust::host_vector<unsigned> hIndicesOfPrisms(maxRaysPerSample,0);
  unsigned midRaysPerSample=0;
  // Calculate Phi Ase for each wavelength
  for(unsigned wave_i = 0; wave_i < numberOfWavelengths; ++wave_i){
    time_t progressStartTime = time(0);
    fprintf(stderr, "\nC Phi_ASE calculation for wavelength %d\n",wave_i);

    //calculation for each sample point
    for(unsigned sample_i = 0; sample_i < hMesh.numberOfSamples; ++sample_i){
      unsigned sampleOffset  = sample_i + hMesh.numberOfSamples * wave_i;
      hRaysPerSample = hRaysPerSampleSave;

      unsigned hRaysPerSampleDump = 0;
      while(true){
	hRaysPerSampleDump = importanceSampling(
            sample_i, reflectionSlices, dMesh, hRaysPerSample, hSigmaA[wave_i], hSigmaE[wave_i],
            raw_pointer_cast(&dImportance[0]), raw_pointer_cast(&dRaysPerPrism[0]),
            distributeRandomly, blockDim, gridDim
            );

        // Prism scheduling for gpu threads
        mapRaysToPrisms(dIndicesOfPrisms, dNumberOfReflections, dRaysPerPrism, dPrefixSum, reflectionSlices, hRaysPerSampleDump, hMesh.numberOfPrisms);

        if(sample_i == 1386){
          thrust::copy(dNumberOfReflections.begin(),dNumberOfReflections.end(),hNumberOfReflections.begin());
          thrust::copy(dIndicesOfPrisms.begin(),dIndicesOfPrisms.end(),hIndicesOfPrisms.begin());
          midRaysPerSample=hRaysPerSample;
        }
        // Start Kernel
        calcSamplePhiAse<<< gridDim, blockDim >>>(
            devMTGPStates,
            dMesh, 
            raw_pointer_cast(&dIndicesOfPrisms[0]), 
            wave_i, 
            raw_pointer_cast(&dNumberOfReflections[0]), 
            raw_pointer_cast(&dImportance[0]),
            hRaysPerSampleDump, 
            raw_pointer_cast(&dPhiAse[0]), 
            raw_pointer_cast(&dPhiAseSquare[0]),
            sample_i, 
            hSigmaA[wave_i], hSigmaE[wave_i]
            );

        // Check square error
        expectation.at(sampleOffset) = calcExpectation(dPhiAse[sampleOffset], dPhiAseSquare[sampleOffset], hRaysPerSampleDump);

        if(expectation[sampleOffset] < expectationThreshold) break;
        if(hRaysPerSample * 10 > maxRaysPerSample)           break;

        // If the threshold is still too high, increase the number of rays and reset the previously calculated value
        hRaysPerSample             *= 10;
        dPhiAse[sampleOffset]       = 0;
        dPhiAseSquare[sampleOffset] = 0;

      }
      // Update progressbar
      if((sample_i+1) % 10 == 0) fancyProgressBar(sample_i, hMesh.numberOfSamples, 60, progressStartTime);

      // get phiASE, normalize it and get dndtAse
      hPhiAse.at(sampleOffset) = dPhiAse[sampleOffset];
      hPhiAse[sampleOffset]   /= hRaysPerSampleDump * 4.0f * 3.14159;
      dndtAse[sampleOffset]    = getDndtAse(hMesh, hSigmaA[wave_i], hSigmaE[wave_i], hPhiAse, sampleOffset);

    }
  }

  std::vector<unsigned> reflectionsPerPrism(hMesh.numberOfPrisms,0);

  
  for(unsigned i=0 ; i<20 ; ++i){
    fprintf(stderr, "%d IndicesOfPrisms: %d\n",i, hIndicesOfPrisms[i]);
    fprintf(stderr, "    reflections:    %d\n",hNumberOfReflections[i]);
  }
  for(unsigned i=0; i<midRaysPerSample; ++i){
    unsigned index = hIndicesOfPrisms[i];
    reflectionsPerPrism[index] = max(reflectionsPerPrism[index],hNumberOfReflections[i]);
  }

  for(unsigned i=0; i<reflectionsPerPrism.size();++i){
    dndtAse[i] = reflectionsPerPrism[i];
  }


  // Free Memory
  cudaFree(devMTGPStates);
  cudaFree(devKernelParams);

  return difftime(time(0),starttime);
}
