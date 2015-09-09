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


/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once

// STL
#include <iostream>

// Alpaka
#include <alpaka/alpaka.hpp>

// HASEonGPU
#include <mesh.hpp>
#include <propagate_ray.hpp>


/**
 * @brief uses a given importance distribution to decide how many rays will be launched from each prism
 *
 * @param *raysDump will contain the number of rays which were mapped to a specific prism
 * 
 * for other parameters, see documentation of importanceSampling()
 */
struct DistributeRaysByImportance {
    template <typename T_Acc,
	      typename T_Mesh>
    ALPAKA_FN_ACC void operator()( T_Acc const &acc,
				   T_Mesh const &mesh,
				   unsigned *raysPerPrism,
				   double *importance,
				   float *sumPhi,
				   const unsigned raysPerSample,
				   unsigned *raysDump) const {

	// unsigned reflection_i = blockIdx.z;
	// unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;

	// int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
	// if(startPrism >= mesh.numberOfPrisms) return;

	// raysPerPrism[startPrism + reflectionOffset] = (unsigned) floor(importance[startPrism + reflectionOffset] / (*sumPhi) * raysPerSample);
	// if(raysPerPrism[startPrism + reflectionOffset] > raysPerSample){
	//     printf("importance: %f sumPhi: %f raysPerPrism[%d]: %d (max %d)\n",importance[startPrism+reflectionOffset],*sumPhi,startPrism+reflectionOffset,raysPerPrism[startPrism+reflectionOffset],raysPerSample);
	// }
	// assert(raysPerPrism[startPrism + reflectionOffset] <= raysPerSample);
	// atomicAdd(raysDump, raysPerPrism[startPrism + reflectionOffset]);
    }

};

// __global__ void distributeRaysByImportance(Mesh mesh,
// 					   unsigned *raysPerPrism,
// 					   double *importance,
// 					   float *sumPhi,
// 					   unsigned raysPerSample,
// 					   unsigned *raysDump){

//   unsigned reflection_i = blockIdx.z;
//   unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;

//   int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
//   if(startPrism >= mesh.numberOfPrisms) return;

//   raysPerPrism[startPrism + reflectionOffset] = (unsigned) floor(importance[startPrism + reflectionOffset] / (*sumPhi) * raysPerSample);
//   if(raysPerPrism[startPrism + reflectionOffset] > raysPerSample){
// 	  printf("importance: %f sumPhi: %f raysPerPrism[%d]: %d (max %d)\n",importance[startPrism+reflectionOffset],*sumPhi,startPrism+reflectionOffset,raysPerPrism[startPrism+reflectionOffset],raysPerSample);
//   }
//   assert(raysPerPrism[startPrism + reflectionOffset] <= raysPerSample);
//   atomicAdd(raysDump, raysPerPrism[startPrism + reflectionOffset]);

// }


/**
 * @brief calculates a first estimate on the importance of each prism, based on a single ray started in the center of each prism
 *
 * @param *importance will contain the initial importance for each prism
 * @param *sumPhi will contain the cumulative sum of the importance values
 *
 * For other parameters, see documentation of importanceSampling()
 *
 */
struct PropagateFromTriangleCenter {
    template <typename T_Acc,
	      typename T_Mesh>
    ALPAKA_FN_ACC void operator()( T_Acc const &acc,
				   T_Mesh const &mesh,
				   double *importance,
				   const unsigned sample_i,
				   const double sigmaA,
				   const double sigmaE) const {
	double gain = 0;
	unsigned reflection_i = alpaka::idx::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0];
	unsigned reflections  = (reflection_i + 1) / 2;
	ReflectionPlane reflectionPlane  = (reflection_i % 2 == 0)? BOTTOM_REFLECTION : TOP_REFLECTION;

	unsigned startPrism = alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0];
	
	if(startPrism >= mesh.numberOfPrisms){
	    return;
	    
	}
	
	unsigned startLevel       = startPrism/(mesh.numberOfTriangles);
	unsigned startTriangle    = startPrism - (mesh.numberOfTriangles * startLevel);
	Point startPoint          = mesh.getCenterPoint(startTriangle, startLevel);
	Point samplePoint         = mesh.getSamplePoint(sample_i);
	unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;

	gain = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, mesh, sigmaA, sigmaE); 
	importance[startPrism + reflectionOffset] = mesh.getBetaVolume(startPrism) * gain;
	if(mesh.getBetaVolume(startPrism) < 0 || gain < 0 || importance[startPrism+reflectionOffset] < 0){
	    printf("beta: %f importance: %f gain: %f\n", mesh.getBetaVolume(startPrism), importance[startPrism + reflectionOffset], gain);
	}

    }
    
};


/**
 * @brief Calculates preImportance which needs only to be done once.
 *
 * @param sample_i         Index of sample point for which importance sampling should be done.
 * @param reflectionSlices 1 + (2 * maxReflections) - Coded information about how many
 *                         reflections a ray should do and on which surface to start.
 * @param dMesh            All information about triangles, points, contants. 
 *                         Is located on device memory. See mesh.h for details.
 * @param sigmaA           Absorption value of the ray.
 * @param sigmaE           Emission value of the ray.
 * @param preImportance    Importance based on gain from test rays (will be returned as pointer).
 * @param blockDim         Number of threads per block.
 * @param gridDim          Number of blocks per grid.
 *
 */
template <typename T_Acc, typename T_Workdiv, typename T_Stream, typename T_Mesh>
ALPAKA_FN_HOST void importanceSamplingPropagation( T_Stream &stream,
						   T_Workdiv const &workdiv,
						   const unsigned sample_i,
						   T_Mesh &dMesh,
						   const double maxSigmaA,
						   const double maxSigmaE,
						   double *dPreImportance){

    PropagateFromTriangleCenter propagateFromTriangleCenter;

    auto const exec (alpaka::exec::create<T_Acc>(workdiv, propagateFromTriangleCenter, dMesh, dPreImportance, sample_i, maxSigmaA, maxSigmaE));    
    alpaka::stream::enqueue(stream, exec);
    
}


/**
 * @brief Calculates importance and ray distribution on prisms.
 *        Based on preImportance and triangle surfaces.
 *
 * @param reflectionSlices   1 + (2 * maxReflections) - Coded information about how many
 *                           reflections a ray should do and on which surface to start.
 * @param dMesh              All information about triangles, points, contants. 
 *                           Is located on device memory. See mesh.h for details.
 * @param raysPerSample      Number of rays that should be distributed on gain medium.
 * @param preImportance      Values calculated from importanceSamplingPropagation.
 * @param importance         Final importance values for this number of rays per sample.
 * @param raysPerPrism       Information about how many rays will start in each prism.
 * @param hSumPhi            Global memory where all threads can sum their phiAse on.
 * @param distributeRandomly In case that not all rays can be distributed by importance
 *                           sampling, the rest will be distributed randomly.
 * @param blockDim           Number of threads per block.
 * @param gridDim            Number of blocks per grid.
 *
 */
template <typename T_Acc, typename T_Workdiv, typename T_Stream, typename T_Mesh>
unsigned importanceSamplingDistribution(T_Stream &stream,
					T_Workdiv &workdiv,
					const unsigned reflectionSlices,
					const T_Mesh &dMesh,
					const unsigned raysPerSample,
					double *dPreImportance,
					double *dImportance,
					unsigned *dRaysPerPrism,
					const float sumPhi,
					const bool distributeRandomly){

    using DevAcc  = alpaka::dev::Dev<T_Stream>;
    using DevHost = alpaka::dev::DevCpu;

    DevAcc  devAcc  (alpaka::dev::getDev(stream));
    DevHost devHost (alpaka::dev::cpu::getDev());

    auto hSumPhi   ( alpaka::mem::buf::alloc<float,    std::size_t, std::size_t, DevHost>(devHost, static_cast<std::size_t>(1)));
    auto hRaysDump ( alpaka::mem::buf::alloc<unsigned, std::size_t, std::size_t, DevHost>(devHost, static_cast<std::size_t>(1)));  
    auto dSumPhi   ( alpaka::mem::buf::alloc<float,    std::size_t, std::size_t, DevAcc> (devAcc, static_cast<std::size_t>(1)));
    auto dRaysDump ( alpaka::mem::buf::alloc<unsigned, std::size_t, std::size_t, DevAcc> (devAcc, static_cast<std::size_t>(1)));  

    alpaka::mem::view::getPtrNative(hSumPhi)[0]   = sumPhi;
    alpaka::mem::view::getPtrNative(hRaysDump)[0] = 0;

    alpaka::mem::view::copy(stream, dSumPhi, hSumPhi, static_cast<std::size_t>(1));
    alpaka::mem::view::copy(stream, dRaysDump, hRaysDump, static_cast<std::size_t>(1));    


    DistributeRaysByImportance distributeRaysByImportance;


    
    // dim3 gridDimReflection(gridDim.x, 1, reflectionSlices);

    // CUDA_CHECK_KERNEL_SYNC(distributeRaysByImportance<<< gridDimReflection, blockDim >>>(deviceMesh, raysPerPrism, preImportance, dSumPhi, raysPerSample, dRaysDump));

    //auto const exec (alpaka::exec::create<T_Acc>(workdiv, distributeRaysByImportance, dMesh, dRaysPerPrism, dPreImportance, dSumPhi, raysPerSample, dRaysDump));
    auto const exec (alpaka::exec::create<T_Acc>(workdiv, distributeRaysByImportance, dMesh, dRaysPerPrism, dPreImportance, alpaka::mem::view::getPtrNative(dSumPhi), raysPerSample, alpaka::mem::view::getPtrNative(dRaysDump)));
    alpaka::stream::enqueue(stream, exec);


    // // Distribute remaining rays randomly if wanted
    // if(distributeRandomly){
    //   CUDA_CHECK_KERNEL_SYNC(distributeRemainingRaysRandomly<<< 200,blockDim >>>(deviceMesh ,raysPerPrism, raysPerSample, dRaysDump));
    //   hRaysDump = raysPerSample;
    // }
    // else {
    //   hRaysDump = copyFromDevice(dRaysDump);
    // }

    // CUDA_CHECK_KERNEL_SYNC(recalculateImportance<<< gridDimReflection, blockDim >>>(deviceMesh, raysPerPrism, hRaysDump, importance));

    // cudaFree(dSumPhi);
    // cudaFree(dRaysDump);

    return alpaka::mem::view::getPtrNative(hRaysDump)[0];
}
