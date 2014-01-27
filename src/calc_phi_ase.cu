#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <curand_kernel.h>
#include <curand_mtgp32_host.h>
#include <cuda_runtime_api.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <write_to_vtk.h>
#include <calc_phi_ase.h>
#include <map_rays_to_prisms.h>
#include <cudachecks.h>
#include <importance_sampling.h>
#include <calc_sample_phi_ase.h>
#include <mesh.h>
#include <progressbar.h> /*progressBar */
#include <logging.h>
#include <types.h>

#define SEED 4321
#define RAY_STEPS 5

double calcMSE(const double phiAse, const double phiAseSquare, const unsigned raysPerSample){
  double a = phiAseSquare / raysPerSample;
  double b = (phiAse / raysPerSample) * (phiAse / raysPerSample);

  return sqrt(abs((a - b) / raysPerSample));
}

std::vector<int> generateRaysPerSampleExpList(int minRaysPerSample, int maxRaysPerSample, int steps){
  std::vector<int> raysPerSample;

  if((minRaysPerSample == maxRaysPerSample) || steps < 2){
    raysPerSample.push_back(minRaysPerSample);
    return raysPerSample;
  }

  for(int i = 0; i < steps; ++i){
    int step_val = minRaysPerSample * pow((maxRaysPerSample / minRaysPerSample), (i / (float)(steps - 1)));
    raysPerSample.push_back(step_val);

  }
  
  return raysPerSample;

}

float calcPhiAse (const unsigned hMinRaysPerSample,
		  const unsigned maxRaysPerSample,
		  const unsigned maxRepetitions,
		  const Mesh& dMesh,
		  const Mesh& hMesh,
		  const std::vector<double>& hSigmaA,
		  const std::vector<double>& hSigmaE,
		  const std::vector<float>& mseThreshold,
		  const bool useReflections,
		  std::vector<float> &phiAse,
		  std::vector<double> &mse,
		  std::vector<unsigned> &totalRays,
		  const unsigned gpu_i,
		  const unsigned minSample_i,
		  const unsigned maxSample_i,
		  float &runtime){

  // Optimization to use more L1 cache
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  cudaSetDevice(gpu_i);

  using thrust::device_vector;
  using thrust::raw_pointer_cast;

  // variable Definitions CPU
  time_t starttime                = time(0);
  unsigned maxReflections         = useReflections ? hMesh.getMaxReflections() : 0;
  unsigned reflectionSlices       = 1 + (2 * maxReflections);
  unsigned numberOfWavelengths    = hSigmaE.size();
  // In some cases distributeRandomly has to be true !
  // Otherwise bad or no ray distribution possible.
  bool distributeRandomly         = true;
  dim3 blockDim(128);             //can't be more than 256 due to restrictions from the Mersenne Twister
  dim3 gridDim(200);              //can't be more than 200 due to restrictions from the Mersenne Twister

  // Divide RaysPerSample range into steps
  std::vector<int>  raysPerSampleList = generateRaysPerSampleExpList(hMinRaysPerSample, maxRaysPerSample, RAY_STEPS);
  std::vector<int>::iterator raysPerSampleIter = raysPerSampleList.begin();

  // Memory allocation/init and copy for device memory
  device_vector<unsigned> dNumberOfReflectionSlices(maxRaysPerSample, 0);
  device_vector<float>    dGainSum            (1, 0);
  device_vector<float>    dGainSumSquare      (1, 0);
  device_vector<unsigned> dRaysPerPrism       (hMesh.numberOfPrisms * reflectionSlices, 1);
  device_vector<unsigned> dPrefixSum          (hMesh.numberOfPrisms * reflectionSlices, 0);
  device_vector<double>   dImportance         (hMesh.numberOfPrisms * reflectionSlices, 0);
  device_vector<double>   dPreImportance      (hMesh.numberOfPrisms * reflectionSlices, 0);
  device_vector<unsigned> dIndicesOfPrisms    (maxRaysPerSample,  0);
  device_vector<double>   dSigmaA             (hSigmaA.begin(),hSigmaA.end());
  device_vector<double>   dSigmaE             (hSigmaE.begin(),hSigmaE.end());

  // CUDA Mersenne twister (can not have more than 200 blocks!)
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.x  * sizeof(curandStateMtgp32)));
  CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));
  CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));
  CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED + minSample_i));

  // Calculate Phi Ase for each wavelength
  for(unsigned wave_i = 0; wave_i < numberOfWavelengths; ++wave_i){

    // Calculation for each sample point
    for(unsigned sample_i = minSample_i; sample_i < maxSample_i; ++sample_i){
      unsigned sampleOffset  = sample_i + hMesh.numberOfSamples * wave_i;
      unsigned hRaysPerSampleDump = 0; 
      raysPerSampleIter = raysPerSampleList.begin();
      bool mseTooHigh=true;

      importanceSamplingPropagation(sample_i,
				    reflectionSlices,
				    dMesh,
				    hSigmaA[wave_i],
				    hSigmaE[wave_i],
				    raw_pointer_cast(&dPreImportance[0]), 
				    blockDim,
				    gridDim);

      float hSumPhi = thrust::reduce(dPreImportance.begin(), dPreImportance.end(),0.);

      while(mseTooHigh){
        CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED + sample_i));
        unsigned run = 0;
        while(run < maxRepetitions && mseTooHigh){
          run++;

	  thrust::copy(dPreImportance.begin(),dPreImportance.end(),dImportance.begin());
          hRaysPerSampleDump = importanceSamplingDistribution(reflectionSlices,
							      dMesh,
							      *raysPerSampleIter,
							      raw_pointer_cast(&dImportance[0]), 
							      raw_pointer_cast(&dRaysPerPrism[0]),
							      hSumPhi,
							      distributeRandomly,
							      blockDim,
							      gridDim);
          
          // Prism scheduling for gpu threads
          mapRaysToPrisms(dIndicesOfPrisms, dNumberOfReflectionSlices, dRaysPerPrism, dPrefixSum, reflectionSlices, hRaysPerSampleDump, hMesh.numberOfPrisms);

          // Start Kernel
          dGainSum[0]       = 0;
          dGainSumSquare[0] = 0;

          if(useReflections){
            calcSampleGainSumWithReflection<<< gridDim, blockDim >>>(devMTGPStates,
								     dMesh, 
								     raw_pointer_cast(&dIndicesOfPrisms[0]), 
								     wave_i, 
								     raw_pointer_cast(&dNumberOfReflectionSlices[0]), 
								     raw_pointer_cast(&dImportance[0]),
								     hRaysPerSampleDump, 
								     raw_pointer_cast(&dGainSum[0]), 
								     raw_pointer_cast(&dGainSumSquare[0]),
								     sample_i, 
								     //hSigmaA[wave_i], 
								     //hSigmaE[wave_i],
                     raw_pointer_cast(&dSigmaA[0]),
                     raw_pointer_cast(&dSigmaE[0]),
                     hSigmaA.size(),
								     raw_pointer_cast(&(device_vector<unsigned> (1,0))[0]));
          }
          else{
            calcSampleGainSum<<< gridDim, blockDim >>>(devMTGPStates,
						       dMesh, 
						       raw_pointer_cast(&dIndicesOfPrisms[0]), 
						       wave_i, 
						       raw_pointer_cast(&dImportance[0]),
						       hRaysPerSampleDump, 
						       raw_pointer_cast(&dGainSum[0]), 
						       raw_pointer_cast(&dGainSumSquare[0]),
						       sample_i, 
						       //hSigmaA[wave_i], 
						       //hSigmaE[wave_i],
                   raw_pointer_cast(&dSigmaA[0]),
                   raw_pointer_cast(&dSigmaE[0]),
                   hSigmaA.size(),
						       raw_pointer_cast(&(device_vector<unsigned> (1,0))[0]));
          }

          float mseTmp = calcMSE(dGainSum[0], dGainSumSquare[0], hRaysPerSampleDump);

          assert(!isnan(dGainSum[0]));
          assert(!isnan(dGainSumSquare[0]));
          assert(!isnan(mseTmp));

          if(mse.at(sampleOffset) > mseTmp){
	    mse.at(sampleOffset) = mseTmp;
	    phiAse.at(sampleOffset) = dGainSum[0]; 
	    phiAse.at(sampleOffset)   /= *raysPerSampleIter * 4.0f * M_PI;
            totalRays.at(sampleOffset) = *raysPerSampleIter;
          }
          if(mse.at(sampleOffset) < mseThreshold.at(wave_i)) mseTooHigh = false;
        }

	// Increase rays per sample or break, when mseThreshold was not met
	raysPerSampleIter++;
	if(raysPerSampleIter == raysPerSampleList.end())
	  break;
	  
      }
      // Update progressbar
      if(verbosity & V_PROGRESS){
        fancyProgressBar(hMesh.numberOfSamples);
      }
    }

  }
  
  // Free Memory
  cudaFree(devMTGPStates);
  cudaFree(devKernelParams);


  runtime = difftime(time(0),starttime);
  return runtime;
}
