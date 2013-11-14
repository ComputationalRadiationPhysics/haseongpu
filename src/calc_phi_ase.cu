#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <assert.h>
#include <vector>
#include <curand_kernel.h>
#include <curand_mtgp32_host.h>
#include <cuda_runtime_api.h>
#include <ctime> /* progressBar */
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

#define SEED 1234

double calcExpectation(const double phiAse, const double phiAseSquare, const unsigned raysPerSample){
  double a = phiAseSquare / raysPerSample;
  double b = (phiAse / raysPerSample) * (phiAse / raysPerSample);

  return sqrt(abs((a - b) / raysPerSample));
}


float calcPhiAse ( unsigned &hRaysPerSample,
		   const unsigned maxRaysPerSample,
		   const Mesh& dMesh,
		   const Mesh& hMesh,
		   const std::vector<double>& hSigmaA,
		   const std::vector<double>& hSigmaE,
		   const std::vector<float>& mseThreshold,
		   const bool useReflections,
		   std::vector<float> &hPhiAse,
		   std::vector<double> &mse,
		   std::vector<unsigned> &totalRays,
		   unsigned gpu_i,
		   unsigned minSample_i,
		   unsigned maxSample_i,
		   float &runtime){

  // Optimization to use more L1 cache
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  cudaSetDevice(gpu_i);

  using thrust::device_vector;
  using thrust::raw_pointer_cast;

  // variable Definitions CPU
  time_t starttime                = time(0);
  unsigned hRaysPerSampleSave     = hRaysPerSample;
  unsigned maxReflections         = useReflections ? hMesh.getMaxReflections() : 0;
  unsigned reflectionSlices       = 1 + (2 * maxReflections);
  unsigned numberOfWavelengths    = hSigmaE.size();
  bool distributeRandomly         = true;
  dim3 blockDim(128);             
  dim3 gridDim(200);              //can't be more than 200 due to restrictions from the Mersenne Twister

  // Memory allocation/init and copy for device memory
  device_vector<unsigned> dNumberOfReflections(maxRaysPerSample,  0);
  device_vector<unsigned> dIndicesOfPrisms    (maxRaysPerSample,  0);
  device_vector<float>    dPhiAse             (hPhiAse.size(), 0);
  device_vector<unsigned> dRaysPerPrism       (hMesh.numberOfPrisms * reflectionSlices, 1);
  device_vector<unsigned> dPrefixSum          (hMesh.numberOfPrisms * reflectionSlices, 0);
  device_vector<double>   dImportance         (hMesh.numberOfPrisms * reflectionSlices, 0);
  device_vector<float>    dPhiAseSquare       (hMesh.numberOfSamples * numberOfWavelengths, 0); //OPTIMIZE: use only 1 value

  // OUTPUT DATA
  // thrust::host_vector<unsigned> hNumberOfReflections(maxRaysPerSample,0);
  // thrust::host_vector<unsigned> hIndicesOfPrisms(maxRaysPerSample,0);
  //thrust::host_vector<unsigned> hRaysPerPrism(hMesh.numberOfPrisms * reflectionSlices, 0);
  //unsigned midRaysPerSample=0;
 
  // CUDA Mersenne twister (can not have more than 200 blocks!)
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.x  * sizeof(curandStateMtgp32)));
  CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));
  CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));
  CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED));

  // Calculate Phi Ase for each wavelength
  for(unsigned wave_i = 0; wave_i < numberOfWavelengths; ++wave_i){
    //time_t progressStartTime = time(0);
    //calculation for each sample point
    for(unsigned sample_i = minSample_i; sample_i < maxSample_i; ++sample_i){
      float mseRunZero = 0.0;
      // MSE BUG TEST
      //for(unsigned sample_i = 71; sample_i < 72; ++sample_i){
      //dout(V_DEBUG) << "SAMPLE " << sample_i << std::endl;
      unsigned sampleOffset  = sample_i + hMesh.numberOfSamples * wave_i;
      hRaysPerSample = hRaysPerSampleSave;

      unsigned hRaysPerSampleDump = 0; while(true){
	unsigned run = 0;
	hRaysPerSampleDump = importanceSampling(
						sample_i, reflectionSlices, dMesh, hRaysPerSample, hSigmaA[wave_i], hSigmaE[wave_i],
						raw_pointer_cast(&dImportance[0]), 
						raw_pointer_cast(&dRaysPerPrism[0]),
						distributeRandomly, blockDim, gridDim
						);

        // Prism scheduling for gpu threads
        mapRaysToPrisms(dIndicesOfPrisms, dNumberOfReflections, dRaysPerPrism, dPrefixSum, reflectionSlices, hRaysPerSampleDump, hMesh.numberOfPrisms);

	// OUTPUT DATA
        // if(sample_i == 0){
        //   thrust::copy(dNumberOfReflections.begin(),dNumberOfReflections.end(),hNumberOfReflections.begin());
        //   thrust::copy(dIndicesOfPrisms.begin(),dIndicesOfPrisms.end(),hIndicesOfPrisms.begin());
	// thrust::copy(dRaysPerPrism.begin(), dRaysPerPrism.end(), hRaysPerPrism.begin());
        //   midRaysPerSample=hRaysPerSample;
        // }

        // Start Kernel
        calcSamplePhiAse<<< gridDim, blockDim >>>( devMTGPStates,
						   dMesh, 
						   raw_pointer_cast(&dIndicesOfPrisms[0]), 
						   wave_i, 
						   raw_pointer_cast(&dNumberOfReflections[0]), 
						   raw_pointer_cast(&dImportance[0]),
						   hRaysPerSampleDump, 
						   raw_pointer_cast(&dPhiAse[0]), 
						   raw_pointer_cast(&dPhiAseSquare[0]),
						   sample_i, 
						   hSigmaA[wave_i], 
						   hSigmaE[wave_i] );


	float mseTmp = calcExpectation(dPhiAse[sampleOffset], dPhiAseSquare[sampleOffset], hRaysPerSampleDump);
	if(run == 0){
	  mseRunZero = mseTmp;
	  run++;
	}

	// MSE TESTs
	//dout(V_DEBUG) << "MSE: " << mseTmp << " with " << hRaysPerSampleDump << " rays,[" << dPhiAse[sampleOffset] << " || " << dPhiAseSquare[sampleOffset] << "]"<< std::endl;
	 //if(mseTmp > mse.at(sampleOffset)){
	 //  double a = dPhiAse[sampleOffset];
	 //  double b = dPhiAseSquare[sampleOffset];
	 //  dout(V_DEBUG) << "\nRaysPerSampleDump: "<< hRaysPerSampleDump << std::endl;
	 //  dout(V_DEBUG) << "RaysPerSample: " << hRaysPerSample << std::endl;
	 //  dout(V_DEBUG) << mseTmp << " > " << mse.at(sampleOffset) << "(" << sample_i << ")\n" << std::endl;
	 //}
        mse.at(sampleOffset) = mseTmp;

        if(mse.at(sampleOffset) < mseThreshold.at(wave_i))     break;
        if(hRaysPerSample * 10 > (unsigned long)maxRaysPerSample)break;

        // If the threshold is still too high, increase the number of rays and reset the previously calculated value
        hRaysPerSample             *= 10;
        dPhiAse[sampleOffset]       = 0;
        dPhiAseSquare[sampleOffset] = 0;

      }
      // Update progressbar
      //if((sample_i+1) % 10 == 0) fancyProgressBar(sample_i-minSample_i, maxSample_i / (gpu_i + 1), 60, progressStartTime);
      fancyProgressBar(maxSample_i / (gpu_i + 1));

      // get phiASE
      hPhiAse.at(sampleOffset) = dPhiAse[sampleOffset];
      hPhiAse.at(sampleOffset)   /= hRaysPerSampleDump * 4.0f * M_PI;
      totalRays.at(sampleOffset)  = hRaysPerSampleDump;
      mse.at(sampleOffset) = mseRunZero;

    }
    

  }

  // JUST OUTPUT
  // std::vector<unsigned> reflectionsPerPrism(hMesh.numberOfPrisms, 0);
  // std::vector<unsigned> raysPerPrism(hMesh.numberOfPrisms, 0);
  
  // for(unsigned i=0; i < midRaysPerSample; ++i){
  //   unsigned index = hIndicesOfPrisms[i];
  //   reflectionsPerPrism[index] = max(reflectionsPerPrism[index], (hNumberOfReflections[i] + 1) / 2);
  // }

  // for(unsigned i=0; i < hMesh.numberOfPrisms; ++i){
  //   for(unsigned j=0; j < reflectionSlices; ++j){
  //     unsigned index = i + hMesh.numberOfPrisms * j;
  //     raysPerPrism[i] += hRaysPerPrism[index];
  //   }
  // }

  //writePrismToVtk(hMesh, reflectionsPerPrism, "octrace_0_reflections", hRaysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, 0);
  //writePrismToVtk(hMesh, raysPerPrism, "octrace_0_rays", hRaysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, 0);

  dout(V_INFO | V_NOLABEL) << "\n" << std::endl;
  // Free Memory
  cudaFree(devMTGPStates);
  cudaFree(devKernelParams);

  runtime = difftime(time(0),starttime);
  return runtime;
}
