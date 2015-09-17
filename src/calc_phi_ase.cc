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

// Alpaka
#include <alpaka/alpaka.hpp>

// HASEonGPU
//#include <write_to_vtk.hpp>
//#include <progressbar.hpp>         /*progressBar */
//#include <logging.hpp>
#include <importance_sampling.hpp>  /* importanceSamplingPropagation, importanceSamplingDistribution */
#include <types.hpp>                /* ExperimentParameter, ComputeParameter, Result */
#include <mesh.hpp>
#include <parser.hpp>               /* parseMesh */
#include <map_rays_to_prisms.hpp>   /* mapRaysToPrisms */
#include <calc_sample_gain_sum.hpp> /* CalcSampleGainSum*/



#define SEED 4321

/******************************************************************************
 * Alpaka utilities
 ******************************************************************************/
template <typename T_Buf, typename T_Extents, typename T_Value>
void initHostBuffer(T_Buf &buf, T_Extents const extents, T_Value const value){
    
    for(T_Extents i = 0; i < extents; ++i){
	alpaka::mem::view::getPtrNative(buf)[static_cast<T_Extents>(i)] = value;
    }

}

template <typename T_Buf, typename T_Extents, typename T_It>
void initHostBuffer(T_Buf &buf, T_Extents const extents, T_It const &begin, T_It const& end){

    
    T_It it = begin;
    (void) end;
    for(unsigned i = 0; i < extents; ++i){
	alpaka::mem::view::getPtrNative(buf)[static_cast<T_Extents>(i)] = *it;
	it++;
    }

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

template <typename T_Value, typename T_Stream, typename T_Buf, typename T_Extents>
void exclusivePrefixSum(T_Stream &stream, T_Buf const &inBuf, T_Buf &outBuf, const T_Extents extents){

    using DevHost = alpaka::dev::DevCpu;
    DevHost devHost (alpaka::dev::cpu::getDev());

    auto hostBuf( alpaka::mem::buf::alloc<T_Value, T_Extents, T_Extents, DevHost>(devHost, extents));
    initHostBuffer(hostBuf, extents, 0);
    alpaka::mem::view::copy(stream, hostBuf, inBuf, extents);
    
    T_Value value = 0;

    for(T_Extents i = 0; i < extents; ++i){
	T_Value prevValue = value;
	value   = value + alpaka::mem::view::getPtrNative(hostBuf)[i];
	alpaka::mem::view::getPtrNative(hostBuf)[i] = prevValue;
    }

    alpaka::mem::view::copy(stream, outBuf, hostBuf, extents);
  
}




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

/******************************************************************************
 * Phi ASE calculations
 ******************************************************************************/
double calcMSE(const double phiAse, const double phiAseSquare, const unsigned raysPerSample){
  double a = phiAseSquare / raysPerSample;
  double b = (phiAse / raysPerSample) * (phiAse / raysPerSample);

  return sqrt(fabs((a - b) / raysPerSample));
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
    using Dim     = alpaka::dim::DimInt<2>;  
    using Size    = std::size_t;
    using Extents = Size;
    using Host    = alpaka::acc::AccCpuSerial<Dim, Size>;
    using Acc     = alpaka::acc::AccCpuOmp2Blocks<Dim, Size>;
    //using Acc     = alpaka::acc::AccCpuSerial<Dim, Size>;
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
    const Mesh<Acc,  DevAcc>  dMesh(parseMesh<Acc,  DevAcc> (compute.inputPath, devAcc));
    const Mesh<Host, DevHost> hMesh(parseMesh<Host, DevHost>(compute.inputPath, devHost));

    
    //alpAka::core::forEachType<alpaka::core::accs::EnabledAccs<Dim, Size> >(AccInfo(), 5ul);

    
    // // Optimization to use more L1 cache
    // cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    // cudaSetDevice(compute.gpu_i);

    // variable Definitions CPU
    const time_t starttime                = time(0);
    const unsigned maxReflections         = experiment.useReflections ? hMesh.getMaxReflections() : 0;
    const unsigned reflectionSlices       = 1 + (2 * maxReflections);

    const alpaka::Vec<Dim, Size> importanceGrid (static_cast<Size>(200), //can't be more than 200 due to restrictions from the Mersenne Twister 
						 static_cast<Size>(reflectionSlices));
    const alpaka::Vec<Dim, Size> grid (static_cast<Size>(200), //can't be more than 200 due to restrictions from the Mersenne Twister, MapPrefixSumToPrism needs number of prisms threads
				       static_cast<Size>(1));
    const alpaka::Vec<Dim, Size> blocks (static_cast<Size>(1),
					 static_cast<Size>(1)); // can't be more than 256 due to restrictions from the Mersenne Twister
                                                                // MUST be 128, since in the kernel we use a bitshift << 7

    auto const importanceWorkdiv(alpaka::workdiv::WorkDivMembers<Dim, Size>(importanceGrid, blocks));
    auto const workdiv(alpaka::workdiv::WorkDivMembers<Dim, Size>(grid, blocks));


    // In some cases distributeRandomly has to be true !
    // Good comment, but next time describe why !
    // Otherwise bad or no ray distribution possible.
    bool distributeRandomly         = true;
    
    // Divide RaysPerSample range into steps
    std::vector<int>  raysPerSampleList = generateRaysPerSampleExpList(experiment.minRaysPerSample,
									     experiment.maxRaysPerSample,
									     compute.adaptiveSteps);
  
    std::vector<int>::iterator raysPerSampleIter = raysPerSampleList.begin();

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
    auto hSigmaE                   ( alpaka::mem::buf::alloc<double  , Size, Extents, DevHost>(devHost, static_cast<Extents>(experiment.sigmaE.size())));
    auto hGlobalOffsetMultiplicator( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevHost>(devHost, static_cast<Extents>(1)));
    
    
  
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
    auto dGlobalOffsetMultiplicator( alpaka::mem::buf::alloc<unsigned, Size, Extents, DevAcc>(devAcc, static_cast<Extents>(1)));

    // std::cout << "numberOfPrisms: " << hMesh.numberOfPrisms << std::endl;
    // std::cout << "reflectionSlices: " << reflectionSlices << std::endl;
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
    initHostBuffer(hGlobalOffsetMultiplicator,1, 0);    

    // Copy host to device buffer
    alpaka::mem::view::copy(stream, dNumberOfReflectionSlices, hNumberOfReflectionSlices, static_cast<Extents>(experiment.maxRaysPerSample));
    alpaka::mem::view::copy(stream, dGainSum, hGainSum, static_cast<Extents>(1));
    alpaka::mem::view::copy(stream, dGainSumSquare, hGainSumSquare, static_cast<Extents>(1));
    alpaka::mem::view::copy(stream, dRaysPerPrism, hRaysPerPrism, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices));
    alpaka::mem::view::copy(stream, dPrefixSum, hPrefixSum, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices));
    alpaka::mem::view::copy(stream, dImportance, hImportance, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices));
    alpaka::mem::view::copy(stream, dPreImportance, hPreImportance, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices));
    alpaka::mem::view::copy(stream, dIndicesOfPrisms, hIndicesOfPrisms, static_cast<Extents>(experiment.maxRaysPerSample));
    alpaka::mem::view::copy(stream, dSigmaA, hSigmaA, static_cast<Extents>(experiment.sigmaA.size()));
    alpaka::mem::view::copy(stream, dSigmaE, hSigmaE, static_cast<Extents>(experiment.sigmaE.size()));
    alpaka::mem::view::copy(stream, dGlobalOffsetMultiplicator, hGlobalOffsetMultiplicator, static_cast<Extents>(1));    		    

    // Kept this lines to have the correct dimension of the vectors for safety
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

    /*****************************************************************************
     * Calculation for each sample point
     ****************************************************************************/
    for(unsigned sample_i = minSample_i; sample_i < maxSample_i; ++sample_i){
    	unsigned raysPerSampleDump = 0; 
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

    	//std::cout << hSumPhi << std::endl;

    	while(mseTooHigh){

    	//     // FIXIT: Reset random seed ?
    	//     // CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, gridDim.x, SEED + sample_i));
	  
     	    unsigned run = 0;
     	    while(run < compute.maxRepetitions && mseTooHigh){
     	 	run++;

    		raysPerSampleDump = importanceSamplingDistribution<Acc>(stream,
    									importanceWorkdiv,
    									workdiv,
    									dMesh,
    									*raysPerSampleIter,
    									alpaka::mem::view::getPtrNative(dPreImportance), 
    									alpaka::mem::view::getPtrNative(dImportance), 
    									alpaka::mem::view::getPtrNative(dRaysPerPrism),
    									hSumPhi,
    									distributeRandomly);
    	 	//std::cout << run << std::endl;

		//std::cout << raysPerSampleDump << std::endl;
		
    		// alpaka::mem::view::copy(stream, hRaysPerPrism, dRaysPerPrism, static_cast<Extents>(hMesh.numberOfPrisms  * reflectionSlices));
    		// for(unsigned i = 0; i < (hMesh.numberOfPrisms  * reflectionSlices); ++i){
    		//     std::cout << alpaka::mem::view::getPtrNative(hRaysPerPrism)[i] << " ";
    		// }

		
		 exclusivePrefixSum<unsigned>(stream, dRaysPerPrism, dPrefixSum, hMesh.numberOfPrisms * reflectionSlices);


    		// alpaka::mem::view::copy(stream, hPrefixSum, dPrefixSum, static_cast<Extents>(hMesh.numberOfPrisms * reflectionSlices));
    		// for(unsigned i = 0; i < (hMesh.numberOfPrisms  * reflectionSlices); ++i){
    		//     std::cout << alpaka::mem::view::getPtrNative(hPrefixSum)[i] << " ";
    		// }

    		// Prism scheduling for gpu threads
    		mapPrefixSumToPrisms<Acc>(stream,
    					  workdiv,
    					  hMesh.numberOfPrisms,
    					  reflectionSlices,
    					  alpaka::mem::view::getPtrNative(dRaysPerPrism),
    					  alpaka::mem::view::getPtrNative(dPrefixSum),
    					  alpaka::mem::view::getPtrNative(dIndicesOfPrisms),
    					  alpaka::mem::view::getPtrNative(dNumberOfReflectionSlices));


    		// alpaka::mem::view::copy(stream, hIndicesOfPrisms, dIndicesOfPrisms, static_cast<Extents>(experiment.maxRaysPerSample));
    		// for(unsigned i = 0; i < (experiment.maxRaysPerSample); ++i){
    		//     std::cout << alpaka::mem::view::getPtrNative(hIndicesOfPrisms)[i] << " ";
    		// }
		
    		alpaka::mem::view::getPtrNative(hGainSum)[0]                   = 0;
    		alpaka::mem::view::getPtrNative(hGainSumSquare)[0]             = 0;
    		//alpaka::mem::view::getPtrNative(hGlobalOffsetMultiplicator)[0] = 0;

    		alpaka::mem::view::copy(stream, dGainSum, hGainSum, static_cast<Extents>(1));
    		alpaka::mem::view::copy(stream, dGainSumSquare, hGainSumSquare, static_cast<Extents>(1));
    		//alpaka::mem::view::copy(stream, dGlobalOffsetMultiplicator, hGlobalOffsetMultiplicator, static_cast<Extents>(1));    		

    		// FIXIT: following 4 lines just to have some temporary global variable in the kernel


		
    		// Start Kernel
    		if(experiment.useReflections){
    		    // CalcSampleGainSumWithReflection calcSampleGainSumWithReflection;

    		    // auto const exec (alpaka::exec::create<Acc>(workdiv,
    		    // 					       calcSampleGainSumWithReflection,
    		    // 					       dMesh, 
    		    // 					       alpaka::mem::view::getPtrNative(dIndicesOfPrisms),
    		    // 					       alpaka::mem::view::getPtrNative(dNumberOfReflectionSlices),
    		    // 					       alpaka::mem::view::getPtrNative(dImportance),
    		    // 					       raysPerSampleDump, 
    		    // 					       alpaka::mem::view::getPtrNative(dGainSum), 
    		    // 					       alpaka::mem::view::getPtrNative(dGainSumSquare),
    		    // 					       sample_i, 
    		    // 					       alpaka::mem::view::getPtrNative(dSigmaA),
    		    // 					       alpaka::mem::view::getPtrNative(dSigmaE),
    		    // 					       experiment.sigmaA.size(),
    		    // 					       alpaka::mem::view::getPtrNative(dGlobalOffsetMultiplicator)));
    		    // alpaka::stream::enqueue(stream, exec);
		    
    		}
    		else{
    		    CalcSampleGainSum calcSampleGainSum;

    		    auto const exec (alpaka::exec::create<Acc>(workdiv,
    		    					       calcSampleGainSum,
    		    					       dMesh, 
    		    					       alpaka::mem::view::getPtrNative(dIndicesOfPrisms), 
    		    					       alpaka::mem::view::getPtrNative(dImportance),
    		    					       raysPerSampleDump, 
    		    					       alpaka::mem::view::getPtrNative(dGainSum), 
    		    					       alpaka::mem::view::getPtrNative(dGainSumSquare),
    		    					       sample_i, 
    		    					       alpaka::mem::view::getPtrNative(dSigmaA),
    		    					       alpaka::mem::view::getPtrNative(dSigmaE),
    		    					       experiment.sigmaA.size(),
    		    					       alpaka::mem::view::getPtrNative(dGlobalOffsetMultiplicator)));
    		    alpaka::stream::enqueue(stream, exec);
    		}


		alpaka::mem::view::copy(stream, hGlobalOffsetMultiplicator, dGlobalOffsetMultiplicator, static_cast<Extents>(1));
    		alpaka::mem::view::copy(stream, hGainSum, dGainSum, static_cast<Extents>(1));
    		alpaka::mem::view::copy(stream, hGainSumSquare, dGainSumSquare, static_cast<Extents>(1));
		
    		double mseTmp = calcMSE(alpaka::mem::view::getPtrNative(hGainSum)[0],
    				       alpaka::mem::view::getPtrNative(hGainSumSquare)[0],
    				       raysPerSampleDump);

    		//std::cout << mseTmp << " " << raysPerSampleDump << " " << alpaka::mem::view::getPtrNative(hGlobalOffsetMultiplicator)[0] << " " << alpaka::mem::view::getPtrNative(hGainSum)[0] << " " << alpaka::mem::view::getPtrNative(hGainSumSquare)[0] << std::endl;
		
    		assert(!isnan(alpaka::mem::view::getPtrNative(hGainSum)[0]));
    		assert(!isnan(alpaka::mem::view::getPtrNative(hGainSumSquare)[0]));
    		assert(!isnan(mseTmp));

    		if(result.mse.at(sample_i) > mseTmp){
    		    //std::cout << alpaka::mem::view::getPtrNative(hGainSum)[0] << std::endl;		    
    		    result.mse.at(sample_i)       = mseTmp;
    		    result.phiAse.at(sample_i)    = alpaka::mem::view::getPtrNative(hGainSum)[0]; 
    		    result.phiAse.at(sample_i)   /= *raysPerSampleIter * 4.0f * M_PI;
    		    result.totalRays.at(sample_i) = *raysPerSampleIter;
    		    std::cout << result.phiAse.at(sample_i) << std::endl;
    		}


		
    	 	if(result.mse.at(sample_i) < experiment.mseThreshold) mseTooHigh = false;
		
     	    }
    	    // Increase rays per sample or break, when mseThreshold was not met
    	    raysPerSampleIter++;
    	    if(raysPerSampleIter == raysPerSampleList.end()) break;
      
	  
     	}

    	// if(verbosity & V_PROGRESS){
    	//     fancyProgressBar(mesh.numberOfSamples);
    	// }

    }

    runtime = difftime(time(0),starttime);
    return runtime;

	
}
