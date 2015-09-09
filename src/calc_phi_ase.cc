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
//#include <calc_sample_gain_sum.hpp>
//#include <mesh.hpp>
//#include <progressbar.hpp> /*progressBar */
//#include <logging.hpp>
#include <importance_sampling.hpp> /* importanceSamplingPropagation */
#include <types.hpp>               /* ExperimentParameter, ComputeParameter, Result */
#include <mesh.hpp>
#include <parser.hpp>              /* parseMesh */

// Alpaka
#include <alpaka/alpaka.hpp>
#include <alpaka/core/EnabledAccs.hpp> 

#define SEED 4321

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






// template< typename T>
// struct has_destructor
// {   
//     /* Has destructor :) */
//     template <typename A> 
//     static std::true_type test(decltype(std::declval<A>().~A()) *) {
//         return std::true_type();
//     }

//     /* Has no destructor :( */
//     template<typename A>
//     static std::false_type test(...) {
//         return std::false_type(); 
//     }

//     /* This will be either `std::true_type` or `std::false_type` */
//     typedef decltype(test<T>(0)) type;

//     static const bool value = type::value; /* Which is it? */
// };

// namespace alpaka {
//     namespace acc {
//         template <typename T1, typename T2>
//         class AccGpuCudaRt;

//         template <typename T1, typename T2>
//         class AccCpuSerial;
//     }
// }

// struct AccInfo{
//     template <typename T_Acc, typename T_Size>
//     void operator()(T_Size const n){

//         using Acc    = T_Acc;
//         using Size   = T_Size;
//         using Dim    = alpaka::dim::DimInt<1u>;
//         using Stream = alpaka::stream::StreamCpuAsync;
        
//         auto nDevs  = alpaka::dev::DevMan<Acc>::getDevCount();

//         // Print device Information
//         for(unsigned dev_i = 0; dev_i < nDevs; ++dev_i){
//             auto dev = alpaka::dev::DevMan<Acc>::getDevByIdx(dev_i);
//             std::cout << "Found device: " << alpaka::dev::getName(dev) << std::endl;
//             std::cout << "              " << alpaka::dev::getMemBytes(dev) << " bytes total" << std::endl;
//             std::cout << "              " << alpaka::dev::getFreeMemBytes(dev) << " bytes free" << std::endl;

//             Stream stream(dev);
//             EmptyKernel kernel;

//             auto const workDiv = alpaka::workdiv::WorkDivMembers<Dim, Size>(1ul,1ul);
//             auto const exec    = alpaka::exec::create<Acc>(workDiv, kernel);


//             alpaka::stream::enqueue(stream, exec);
      
//         }

        
//     }
    
    
// };

template <typename T_Stream, typename T_Buf, typename T_Value, typename T_Extents>
T_Value reduce(T_Stream &stream, T_Buf const &buf, T_Extents const extents, T_Value const init){
    using DevHost = alpaka::dev::DevCpu;
    DevHost devHost (alpaka::dev::cpu::getDev());

    // Create host buffer and copy buffer to host buffer
    auto hostBuf( alpaka::mem::buf::alloc<T_Value, T_Extents, T_Extents, DevHost>(devHost, extents));    
    alpaka::mem::view::copy(stream, hostBuf, buf, extents);

    // Reduce host buffer
    T_Value result = init;
    for(T_Extents i = 0; i < extents; ++i){
	result += alpaka::mem::view::getPtrNative(hostBuf)[i];
	
    }
    return result;

}

template <typename T_Buf, typename T_Extents, typename T_Value>
void initHostBuffer(T_Buf &buf, T_Extents const extents, T_Value const value){
    
    for(T_Extents i = 0; i < extents; ++i){
	alpaka::mem::view::getPtrNative(buf)[static_cast<T_Extents>(i)] = value;
    }

}

template <typename T_Buf, typename T_Extents, typename T_It>
void initHostBuffer(T_Buf &buf, T_Extents const extents, T_It const &begin, T_It const& end){


    T_It it = begin;
    for(unsigned i = 0; i < extents; ++i){
	alpaka::mem::view::getPtrNative(buf)[static_cast<T_Extents>(i)] = *it;
	it++;
    }

}


float calcPhiAse ( const ExperimentParameters& experiment,
		   const ComputeParameters& compute,
		   Result& result,
		   const unsigned minSample_i,
		   const unsigned maxSample_i,
		   float &runtime ){

    /*****************************************************************************
     * Choose device
     ****************************************************************************/
    // Set types 
    using Dim     = alpaka::dim::DimInt<3>;  
    using Size    = std::size_t;
    using Extents = Size;
    using Host    = alpaka::acc::AccCpuSerial<Dim, Size>;
    using Acc     = alpaka::acc::AccCpuOmp2Threads<Dim, Size>;
    using Stream  = alpaka::stream::StreamCpuSync;
    using DevAcc  = alpaka::dev::Dev<Acc>;
    using DevHost = alpaka::dev::DevCpu;


    // Get the first device
    DevAcc  devAcc  (alpaka::dev::DevMan<Acc>::getDevByIdx(0));
    DevHost devHost (alpaka::dev::cpu::getDev());
    Stream  stream  (devAcc);

    /*****************************************************************************
     * Parse mesh
     ****************************************************************************/
    Mesh<Acc,  DevAcc>  dMesh = parseMesh<Acc,  DevAcc> (compute.inputPath, devAcc);
    Mesh<Host, DevHost> hMesh = parseMesh<Host, DevHost>(compute.inputPath, devHost);

    
    //alpAka::core::forEachType<alpaka::core::accs::EnabledAccs<Dim, Size> >(AccInfo(), 5ul);

    
    // // Optimization to use more L1 cache
    // cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    // cudaSetDevice(compute.gpu_i);

    // variable Definitions CPU
    time_t starttime                = time(0);
    unsigned maxReflections         = experiment.useReflections ? hMesh.getMaxReflections() : 0;
    unsigned reflectionSlices       = 1 + (2 * maxReflections);


    alpaka::Vec<Dim, Size> importanceGrid (static_cast<Size>(reflectionSlices), // z
					   static_cast<Size>(1), // y
					   static_cast<Size>(200)); // x
    alpaka::Vec<Dim, Size> grid (static_cast<Size>(1),
				 static_cast<Size>(1),
				 static_cast<Size>(200));
    alpaka::Vec<Dim, Size> blocks (static_cast<Size>(1),
				   static_cast<Size>(1),
				   static_cast<Size>(128));

    auto const importanceWorkdiv(alpaka::workdiv::WorkDivMembers<Dim, Size>(importanceGrid, blocks));
    auto const workdiv(alpaka::workdiv::WorkDivMembers<Dim, Size>(grid, blocks));

    // // In some cases distributeRandomly has to be true !
    // // Otherwise bad or no ray distribution possible.
    bool distributeRandomly         = true;
    // dim3 blockDim(128);             //can't be more than 256 due to restrictions from the Mersenne Twister
    //                                 // MUST be 128, since in the kernel we use a bitshift << 7
    // dim3 gridDim(200);              //can't be more than 200 due to restrictions from the Mersenne Twister



    
    // Divide RaysPerSample range into steps
    std::vector<int>  raysPerSampleList = generateRaysPerSampleExpList(experiment.minRaysPerSample,
								       experiment.maxRaysPerSample,
								       compute.adaptiveSteps);
  
    std::vector<int>::iterator raysPerSampleIter = raysPerSampleList.begin();

    // This lines can be reduced, since maxSigma is stored in experiment
    // // // Calc max sigmaA / sigmaE
    // // double maxSigmaE = 0;
    // // double maxSigmaA = 0;
    // // for(unsigned i = 0; i < hSigmaE.size(); ++i){
    // //   if(hSigmaE.at(i) > maxSigmaE){
    // //     maxSigmaE = hSigmaE.at(i);
    // //     maxSigmaA = hSigmaA.at(i);
    // //   }
    // // }


  

    /*****************************************************************************
     * Memory allocation/init and copy for device memory
     ****************************************************************************/
 

  
    // Define memory buffers
    auto hNumberOfReflectionSlices ( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevHost>(devHost, static_cast<Extents>(experiment.maxRaysPerSample)));
    auto hGainSum                  ( alpaka::mem::buf::alloc<float,    Size, Extents, DevHost>(devHost, static_cast<Extents>(1)));
    auto hGainSumSquare            ( alpaka::mem::buf::alloc<float,    Size, Extents, DevHost>(devHost, static_cast<Extents>(1)));
    auto hRaysPerPrism             ( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevHost>(devHost, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices)));
    auto hPrefixSum                ( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevHost>(devHost, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices)));
    auto hImportance               ( alpaka::mem::buf::alloc<double  , Size, Extents, DevHost>(devHost, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices)));
    auto hPreImportance            ( alpaka::mem::buf::alloc<double  , Size, Extents, DevHost>(devHost, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices)));
    auto hIndicesOfPrisms          ( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevHost>(devHost, static_cast<Extents>(experiment.maxRaysPerSample)));
    auto hSigmaA                   ( alpaka::mem::buf::alloc<double  , Size, Extents, DevHost>(devHost, static_cast<Extents>(experiment.sigmaA.size())));
    auto hSigmaE                   ( alpaka::mem::buf::alloc<double  , Size, Extents, DevHost>(devHost,static_cast<Extents>(experiment.sigmaE.size())));
  
    auto dNumberOfReflectionSlices ( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevAcc>(devAcc, static_cast<Extents>(experiment.maxRaysPerSample)));
    auto dGainSum                  ( alpaka::mem::buf::alloc<float,    Size, Extents, DevAcc>(devAcc, static_cast<Extents>(1)));
    auto dGainSumSquare            ( alpaka::mem::buf::alloc<float,    Size, Extents, DevAcc>(devAcc, static_cast<Extents>(1)));
    auto dRaysPerPrism             ( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevAcc>(devAcc, static_cast<Extents>(hMesh.numberOfPrisms  * reflectionSlices)));
    auto dPrefixSum                ( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevAcc>(devAcc, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices)));
    auto dImportance               ( alpaka::mem::buf::alloc<double  , Size, Extents, DevAcc>(devAcc, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices)));
    auto dPreImportance            ( alpaka::mem::buf::alloc<double  , Size, Extents, DevAcc>(devAcc, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices)));
    auto dIndicesOfPrisms          ( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevAcc>(devAcc, static_cast<Extents>(experiment.maxRaysPerSample)));
    auto dSigmaA                   ( alpaka::mem::buf::alloc<double  , Size, Extents, DevAcc>(devAcc, static_cast<Extents>(experiment.sigmaA.size())));
    auto dSigmaE                   ( alpaka::mem::buf::alloc<double  , Size, Extents, DevAcc>(devAcc, static_cast<Extents>(experiment.sigmaE.size())));

    // std::cout << "numberOfPrisms: " << hMesh.numberOfPrisms << std::endl;
    // std::cout << "reflectionsSlices: " << reflectionSlices << std::endl;
    // std::cout << "maxRaysPerSample: " << experiment.maxRaysPerSample << std::endl;
    // std::cout << "Memory was allocated: " << hMesh.numberOfPrisms * reflectionSlices << std::endl;

    // Init memory buffer
    initHostBuffer(hNumberOfReflectionSlices, experiment.maxRaysPerSample, 0);
    initHostBuffer(hGainSum,                  1, 0);
    initHostBuffer(hGainSumSquare,            1, 0);
    initHostBuffer(hRaysPerPrism,             hMesh.numberOfPrisms * reflectionSlices, 1);
    initHostBuffer(hPrefixSum,                hMesh.numberOfPrisms * reflectionSlices, 0);
    initHostBuffer(hImportance,               hMesh.numberOfPrisms * reflectionSlices, 0);
    initHostBuffer(hPreImportance,            hMesh.numberOfPrisms * reflectionSlices, 0);
    initHostBuffer(hIndicesOfPrisms,          experiment.maxRaysPerSample, 0);
    initHostBuffer(hSigmaA,                   experiment.sigmaA.size(), experiment.sigmaA.begin(), experiment.sigmaA.end());
    initHostBuffer(hSigmaE,                   experiment.sigmaA.size(), experiment.sigmaE.begin(), experiment.sigmaA.end());

    // Copy host to device buffer
    alpaka::mem::view::copy(stream, dNumberOfReflectionSlices, hNumberOfReflectionSlices, static_cast<Extents>(experiment.maxRaysPerSample));
    alpaka::mem::view::copy(stream, dGainSum, hGainSum, static_cast<Extents>(1));
    alpaka::mem::view::copy(stream, dGainSumSquare, hGainSumSquare, static_cast<Extents>(1));
    alpaka::mem::view::copy(stream, dRaysPerPrism, hRaysPerPrism, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices));
    alpaka::mem::view::copy(stream, hPrefixSum, dPrefixSum, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices));
    alpaka::mem::view::copy(stream, hImportance, dImportance, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices));
    alpaka::mem::view::copy(stream, hPreImportance, dPreImportance, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices));
    alpaka::mem::view::copy(stream, hIndicesOfPrisms, dIndicesOfPrisms, static_cast<Extents>(experiment.maxRaysPerSample));
    alpaka::mem::view::copy(stream, hSigmaA, dSigmaA, static_cast<Extents>(experiment.sigmaA.size()));
    alpaka::mem::view::copy(stream, hSigmaE, dSigmaE, static_cast<Extents>(experiment.sigmaE.size()));

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


    /*****************************************************************************
     * Calculation for each sample point
     ****************************************************************************/
    for(unsigned sample_i = minSample_i; sample_i < maxSample_i; ++sample_i){
	unsigned hRaysPerSampleDump = 0; 
	raysPerSampleIter           = raysPerSampleList.begin();
	bool mseTooHigh             = true;

	importanceSamplingPropagation<Acc>(stream,
					   importanceWorkdiv,
					   sample_i,
					   dMesh,
					   experiment.maxSigmaA,
					   experiment.maxSigmaE,
					   alpaka::mem::view::getPtrNative(dPreImportance));
      
	float hSumPhi = reduce(stream, dPreImportance, hMesh.numberOfPrisms * reflectionSlices, 0.);

	std::cout << hSumPhi << std::endl;

	while(mseTooHigh){

	    // Reset random seed ?
	    //     CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED + sample_i));
	  
	    unsigned run = 0;
	    while(run < compute.maxRepetitions && mseTooHigh){
		run++;

		hRaysPerSampleDump = importanceSamplingDistribution<Acc>(stream,
									 importanceWorkdiv,
									 reflectionSlices,
									 dMesh,
									 *raysPerSampleIter,
									 alpaka::mem::view::getPtrNative(dPreImportance), 
									 alpaka::mem::view::getPtrNative(dImportance), 
									 alpaka::mem::view::getPtrNative(dRaysPerPrism),
									 hSumPhi,
									 distributeRandomly);
          
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
	    }

	    //     // Increase rays per sample or break, when mseThreshold was not met
	    //     raysPerSampleIter++;
	    //     if(raysPerSampleIter == raysPerSampleList.end())
	    // 	break;
      
	  
	}

	//    if(verbosity & V_PROGRESS){
	//      fancyProgressBar(mesh.numberOfSamples);
	//    }

    }
    
    // // Free Memory
    // cudaFree(devMTGPStates);
    // cudaFree(devKernelParams);
    //alpaka::mem::free(hNumberOfReflectionSlices ( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevHost>(devHost, static_cast<Extents>(experiment.maxRaysPerSample)));
  

    // runtime = difftime(time(0),starttime);
    // return runtime;
    return .0;
}
