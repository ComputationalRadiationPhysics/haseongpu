/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

// CLIB
#include <cmath>
#include <cassert>
#include <cstdlib>

// STL
#include <vector>
#include <utility> /* std::forward */
#include <type_traits>

// CUDA
// #include <curand_kernel.h> /*curand_uniform*/
// #include <curand_mtgp32_host.h>
// #include <cuda_runtime_api.h>
// #include <thrust/device_vector.h>
// #include <vector_types.h> /* dim3 */

// HASEonGPU
//#include <write_to_vtk.hpp>
//#include <calc_phi_ase.hpp>
//#include <map_rays_to_prisms.hpp>
//#include <cudachecks.hpp>
//#include <importance_sampling.hpp>
//#include <calc_sample_gain_sum.hpp>
//#include <mesh.hpp>
//#include <progressbar.hpp> /*progressBar */
//#include <logging.hpp>
#include <types.hpp> /* ExperimentParameter, ComputeParameter, Result */

// Alpaka
#include <alpaka/alpaka.hpp>
//#include <alpaka/accs/EnabledAccs.hpp> 

#define SEED 4321

// double calcMSE(const double phiAse, const double phiAseSquare, const unsigned raysPerSample){
//   double a = phiAseSquare / raysPerSample;
//   double b = (phiAse / raysPerSample) * (phiAse / raysPerSample);

//   return sqrt(abs((a - b) / raysPerSample));
// }

// std::vector<int> generateRaysPerSampleExpList(int minRaysPerSample, int maxRaysPerSample, int steps){
//   std::vector<int> raysPerSample;

//   if((minRaysPerSample == maxRaysPerSample) || steps < 2){
//     raysPerSample.push_back(minRaysPerSample);
//     return raysPerSample;
//   }

//   for(int i = 0; i < steps; ++i){
//     int step_val = minRaysPerSample * pow((maxRaysPerSample / minRaysPerSample), (i / (float)(steps - 1)));
//     raysPerSample.push_back(step_val);

//   }
  
//   return raysPerSample;

// }



struct EmptyKernel {
    template <typename T_Acc>
    ALPAKA_FN_ACC void operator()(T_Acc const & acc) const {

	    
    }

};


template< typename T>
struct has_destructor
{   
    /* Has destructor :) */
    template <typename A> 
    static std::true_type test(decltype(std::declval<A>().~A()) *) {
        return std::true_type();
    }

    /* Has no destructor :( */
    template<typename A>
    static std::false_type test(...) {
        return std::false_type(); 
    }

    /* This will be either `std::true_type` or `std::false_type` */
    typedef decltype(test<T>(0)) type;

    static const bool value = type::value; /* Which is it? */
};

namespace alpaka {
    namespace acc {
        template <typename T1, typename T2>
        class AccGpuCudaRt;

        template <typename T1, typename T2>
        class AccCpuSerial;
    }
}

float calcPhiAse ( const ExperimentParameters& experiment,
		   const ComputeParameters& compute,
		   Result& result,
		   const unsigned minSample_i,
		   const unsigned maxSample_i,
		   float &runtime ){

     // Set types 
     using Dim  = alpaka::dim::DimInt<1u>;
     using Size = std::size_t;

     // Cpu Serial
     using Acc  = alpaka::acc::AccCpuSerial<Dim, Size>;
     using Stream = alpaka::stream::StreamCpuSync;

     // Cpu Threads
     // using Acc  = alpaka::acc::AccCpuThreads<Dim, Size>;
     // using Stream = alpaka::stream::StreamCpuSync;


     // Cpu OpenMP2 Blocks
     // using Acc    = alpaka::acc::AccCpuOmp2Blocks<Dim, Size>;
     // using Stream = alpaka::stream::StreamCpuAsync;

     // Cpu OpenMP2 Threads
     // using Acc    = alpaka::acc::AccCpuOmp2Threads<Dim, Size>;
     // using Stream = alpaka::stream::StreamCpuAsync;

     // Cpu OpenMP4
     // using Acc    = alpaka::acc::AccCpuOmp4<Dim, Size>;
     // using Stream = alpaka::stream::StreamCpuAsync;
     
     // CUDA
     // Is not working. I think there need to be special options set
     // in the cmake file. See alpaka examples.
     // using Acc    = alpaka::acc::AccGpuCudaRt<Dim, Size>;
     // using Stream = alpaka::stream::StreamCpuAsync;

     std::cout << has_destructor<alpaka::acc::AccCpuSerial<Dim, Size> >::value << std::endl;
     //std::cout << has_destructor<alpaka::acc::AccGpuCudaRt<Dim, Size> >::value << std::endl;

     // Get host and device
     auto host = alpaka::dev::cpu::getDev();
     auto nDevs  = alpaka::dev::DevMan<Acc>::getDevCount();

     // Print some device Information
     for(unsigned dev_i = 0; dev_i < nDevs; ++dev_i){
       auto dev = alpaka::dev::DevMan<Acc>::getDevByIdx(dev_i);
       std::cout << "Found device: " << alpaka::dev::getName(dev) << std::endl;
       std::cout << "              " << alpaka::dev::getMemBytes(dev) << " bytes total" << std::endl;
       std::cout << "              " << alpaka::dev::getFreeMemBytes(dev) << " bytes free" << std::endl;

       Stream stream(dev);
       EmptyKernel kernel;

       auto const workDiv = alpaka::workdiv::WorkDivMembers<Dim, Size>(1ul,1ul);
       auto const exec    = alpaka::exec::create<Acc>(workDiv, kernel);


       // TODO get type of Exec

       alpaka::stream::enqueue(stream, exec);
      
     }

     //alpaka::core::forEachType<alpaka::examples::accs::EnabledAccs<Dim, Size> >(AccInfo());

    

    
  // // Optimization to use more L1 cache
  // cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  // cudaSetDevice(compute.gpu_i);

  // using thrust::device_vector;
  // using thrust::raw_pointer_cast;

  // // variable Definitions CPU
  // time_t starttime                = time(0);
  // unsigned maxReflections         = experiment.useReflections ? mesh.getMaxReflections() : 0;
  // unsigned reflectionSlices       = 1 + (2 * maxReflections);
  // // In some cases distributeRandomly has to be true !
  // // Otherwise bad or no ray distribution possible.
  // bool distributeRandomly         = true;
  // dim3 blockDim(128);             //can't be more than 256 due to restrictions from the Mersenne Twister
  //                                 // MUST be 128, since in the kernel we use a bitshift << 7
  // dim3 gridDim(200);              //can't be more than 200 due to restrictions from the Mersenne Twister

  // // Divide RaysPerSample range into steps
  // std::vector<int>  raysPerSampleList = generateRaysPerSampleExpList(experiment.minRaysPerSample,
  // 								     experiment.maxRaysPerSample,
  // 								     compute.adaptiveSteps);
  
  // std::vector<int>::iterator raysPerSampleIter = raysPerSampleList.begin();

  // // // Calc max sigmaA / sigmaE
  // // double maxSigmaE = 0;
  // // double maxSigmaA = 0;
  // // for(unsigned i = 0; i < hSigmaE.size(); ++i){
  // //   if(hSigmaE.at(i) > maxSigmaE){
  // //     maxSigmaE = hSigmaE.at(i);
  // //     maxSigmaA = hSigmaA.at(i);
  // //   }
  // // }

  // // Memory allocation/init and copy for device memory
  // device_vector<unsigned> dNumberOfReflectionSlices(experiment.maxRaysPerSample, 0);
  // device_vector<float>    dGainSum            (1, 0);
  // device_vector<float>    dGainSumSquare      (1, 0);
  // device_vector<unsigned> dRaysPerPrism       (mesh.numberOfPrisms * reflectionSlices, 1);
  // device_vector<unsigned> dPrefixSum          (mesh.numberOfPrisms * reflectionSlices, 0);
  // device_vector<double>   dImportance         (mesh.numberOfPrisms * reflectionSlices, 0);
  // device_vector<double>   dPreImportance      (mesh.numberOfPrisms * reflectionSlices, 0);
  // device_vector<unsigned> dIndicesOfPrisms    (experiment.maxRaysPerSample,  0);
  // device_vector<double>   dSigmaA             (experiment.sigmaA.begin(), experiment.sigmaA.end());
  // device_vector<double>   dSigmaE             (experiment.sigmaE.begin(),experiment.sigmaE.end());

  // // CUDA Mersenne twister (can not have more than 200 blocks!)
  // curandStateMtgp32 *devMTGPStates;
  // mtgp32_kernel_params *devKernelParams;
  // CUDA_CALL(cudaMalloc((void **)&devMTGPStates, gridDim.x  * sizeof(curandStateMtgp32)));
  // CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));
  // CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));
  // CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED + minSample_i));

  // // Calculation for each sample point
  // for(unsigned sample_i = minSample_i; sample_i < maxSample_i; ++sample_i){
  //   unsigned hRaysPerSampleDump = 0; 
  //   raysPerSampleIter = raysPerSampleList.begin();
  //   bool mseTooHigh=true;

  //   importanceSamplingPropagation(sample_i,
  // 				  reflectionSlices,
  // 				  mesh,
  // 				  experiment.maxSigmaA,
  // 				  experiment.maxSigmaE,
  // 				  raw_pointer_cast(&dPreImportance[0]), 
  // 				  blockDim,
  // 				  gridDim);

  //   float hSumPhi = thrust::reduce(dPreImportance.begin(), dPreImportance.end(),0.);

  //   while(mseTooHigh){
  //     CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED + sample_i));
  //     unsigned run = 0;
  //     while(run < compute.maxRepetitions && mseTooHigh){
  // 	run++;

  // 	hRaysPerSampleDump = importanceSamplingDistribution(reflectionSlices,
  // 							    mesh,
  // 							    *raysPerSampleIter,
  // 							    raw_pointer_cast(&dPreImportance[0]), 
  // 							    raw_pointer_cast(&dImportance[0]), 
  // 							    raw_pointer_cast(&dRaysPerPrism[0]),
  // 							    hSumPhi,
  // 							    distributeRandomly,
  // 							    blockDim,
  // 							    gridDim);
          
  // 	// Prism scheduling for gpu threads
  // 	mapRaysToPrisms(dIndicesOfPrisms, dNumberOfReflectionSlices, dRaysPerPrism, dPrefixSum, reflectionSlices, hRaysPerSampleDump, mesh.numberOfPrisms);

  // 	// Start Kernel
  // 	dGainSum[0]       = 0;
  // 	dGainSumSquare[0] = 0;

  // 	if(experiment.useReflections){
  // 	  calcSampleGainSumWithReflection<<< gridDim, blockDim >>>(devMTGPStates,
  // 								   mesh, 
  // 								   raw_pointer_cast(&dIndicesOfPrisms[0]), 
  // 								   raw_pointer_cast(&dNumberOfReflectionSlices[0]), 
  // 								   raw_pointer_cast(&dImportance[0]),
  // 								   hRaysPerSampleDump, 
  // 								   raw_pointer_cast(&dGainSum[0]), 
  // 								   raw_pointer_cast(&dGainSumSquare[0]),
  // 								   sample_i, 
  // 								   raw_pointer_cast(&dSigmaA[0]),
  // 								   raw_pointer_cast(&dSigmaE[0]),
  // 								   experiment.sigmaA.size(),
  // 								   raw_pointer_cast(&(device_vector<unsigned> (1,0))[0]));
  // 	}
  // 	else{
  // 	  calcSampleGainSum<<< gridDim, blockDim >>>(devMTGPStates,
  // 						     mesh, 
  // 						     raw_pointer_cast(&dIndicesOfPrisms[0]), 
  // 						     raw_pointer_cast(&dImportance[0]),
  // 						     hRaysPerSampleDump, 
  // 						     raw_pointer_cast(&dGainSum[0]), 
  // 						     raw_pointer_cast(&dGainSumSquare[0]),
  // 						     sample_i, 
  // 						     raw_pointer_cast(&dSigmaA[0]),
  // 						     raw_pointer_cast(&dSigmaE[0]),
  // 						     experiment.sigmaA.size(),
  // 						     raw_pointer_cast(&(device_vector<unsigned> (1,0))[0]));
  // 	}

  // 	float mseTmp = calcMSE(dGainSum[0], dGainSumSquare[0], hRaysPerSampleDump);

  // 	assert(!isnan(dGainSum[0]));
  // 	assert(!isnan(dGainSumSquare[0]));
  // 	assert(!isnan(mseTmp));

  // 	if(result.mse.at(sample_i) > mseTmp){
  // 	  result.mse.at(sample_i) = mseTmp;
  // 	  result.phiAse.at(sample_i) = dGainSum[0]; 
  // 	  result.phiAse.at(sample_i)   /= *raysPerSampleIter * 4.0f * M_PI;
  // 	  result.totalRays.at(sample_i) = *raysPerSampleIter;
  // 	}
  // 	if(result.mse.at(sample_i) < experiment.mseThreshold) mseTooHigh = false;
  //     }

  //     // Increase rays per sample or break, when mseThreshold was not met
  //     raysPerSampleIter++;
  //     if(raysPerSampleIter == raysPerSampleList.end())
  // 	break;
      
	  
  //   }

  //    if(verbosity & V_PROGRESS){
  //      fancyProgressBar(mesh.numberOfSamples);
  //    }

  // }
    
  // // Free Memory
  // cudaFree(devMTGPStates);
  // cudaFree(devKernelParams);

  // runtime = difftime(time(0),starttime);
  // return runtime;
    return .0;
}
