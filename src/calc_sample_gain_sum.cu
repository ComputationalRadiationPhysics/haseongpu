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


#include <cassert> /* assert */

#include <mesh.hpp>
#include <geometry.hpp> /* generateRay */
#include <propagate_ray.hpp> /* propagateRay */
#include <reflection.hpp> /* ReflectionPlane */

#include <curand_kernel.h> /*curand_uniform*/

/**
 * @brief get the offset for accessing indicesOfPrisms and numberOfReflectionSlices (slow!).
 *
 * @param blockOffset shared memory location that holds the offset for the whole block (4 warps)
 * @param raysPerSample number of raysPerSample (can be any number higher than raysPerSample/warpsize)
 * @param globalOffsetMultiplicator is incremented by 1 each time a warp asks for a new workload
 * @return an unused offset in the global arrays indicesOfPrisms/numberOfReflectionSlices
 *
 */
__device__ unsigned getRayNumberWarpbased(unsigned* blockOffset,unsigned raysPerSample, unsigned *globalOffsetMultiplicator){
	// if this is warpID 0
	if((threadIdx.x &31) == 0){
		//get a new offset for the warp (threadId % 32)
		blockOffset[(threadIdx.x>>5)] = atomicInc(globalOffsetMultiplicator,raysPerSample);
	}
	__syncthreads();

	// multiply blockoffset by 32 (size of warp)
	return (threadIdx.x &31) + (blockOffset[(threadIdx.x>>5)] <<5) ;

}

/**
 * @brief get the offset for accessing indicesOfPrisms and numberOfReflectionSlices.
 *        Warning: works only for a blocksize of 128 threads!
 *
 * @param blockOffset shared memory location that holds the offset for the whole block
 * @param raysPerSample number of raysPerSample (can be any number higher than raysPerSample/blocksize)
 * @param globalOffsetMultiplicator is incremented by 1 each time a block asks for a new workload
 * @return an unused offset in the global arrays indicesOfPrisms/numberOfReflectionSlices
 *
 */
__device__ unsigned getRayNumberBlockbased(unsigned* blockOffset,unsigned raysPerSample,unsigned *globalOffsetMultiplicator){
	// The first thread in the threadblock increases the globalOffsetMultiplicator (without real limit) 
	if(threadIdx.x == 0){
		//blockOffset is the new value of the globalOffsetMultiplicator
		blockOffset[0] = atomicInc(globalOffsetMultiplicator,raysPerSample);
	}
	__syncthreads();

	//multiply blockOffset by 128 (size of the threadblock) 
	return threadIdx.x + (blockOffset[0] <<7) ;
}

/**
 * @brief get a random number from [0..length)
 *
 * @param length the maximum number to return (exclusive)
 * @param globalState State for random number generation (mersenne twister).
 *                    The state need to be initialized before. See
 *                    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MTGP/
 *                    for more information.
 *
 * @return a random number
 *
 */
__device__ __inline__ unsigned genRndSigmas(unsigned length, curandStateMtgp32* globalState) {
  return unsigned(curand_uniform(&globalState[blockIdx.x])*(length-1));
}

__global__ void calcSampleGainSumWithReflection(curandStateMtgp32* globalState,
						const Mesh mesh, 
						const unsigned* indicesOfPrisms, 
						const unsigned reflection_i,
						const double* importance,
						const unsigned raysPerSample,
						float *gainSum, 
						float *gainSumSquare,
						const unsigned sample_i,
						const double *sigmaA, 
						const double *sigmaE,
						const unsigned maxInterpolation,
						unsigned *globalOffsetMultiplicator
						) {

  int rayNumber = 0;
  double gainSumTemp = 0;
  double gainSumSquareTemp = 0;
  Point samplePoint = mesh.getSamplePoint(sample_i);
  __shared__ unsigned blockOffset[4]; // 4 in case of warp-based raynumber

  // One thread can compute multiple rays
  while (true) {
	// the whole block gets a new offset (==workload)
	rayNumber = getRayNumberBlockbased(blockOffset,raysPerSample,globalOffsetMultiplicator);
	if(rayNumber >= raysPerSample) break;

    // Get triangle/prism to start ray from
    unsigned startPrism             = indicesOfPrisms[rayNumber];
    unsigned reflections            = (reflection_i + 1) / 2;
    ReflectionPlane reflectionPlane = (reflection_i % 2 == 0) ? BOTTOM_REFLECTION : TOP_REFLECTION;
    unsigned startLevel             = startPrism / mesh.numberOfTriangles;
    unsigned startTriangle          = startPrism - (mesh.numberOfTriangles * startLevel);
    Point startPoint                = mesh.genRndPoint(startTriangle, startLevel, globalState);
	
	//get a random index in the wavelength array
    unsigned sigma_i                = genRndSigmas(maxInterpolation, globalState);

    // Calculate reflections as different ray propagations
    double gain    = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, mesh, sigmaA[sigma_i], sigmaE[sigma_i]);

	// include the stimulus from the starting prism and the importance of that ray
    gain          *= mesh.getBetaVolume(startPrism) * importance[startPrism];
    
    assert(!isnan(mesh.getBetaVolume(startPrism)));
    assert(!isnan(importance[startPrism]));
    assert(!isnan(gain));

    gainSumTemp       += gain;
    gainSumSquareTemp += gain * gain;


  }
  atomicAdd(&(gainSum[0]), float(gainSumTemp));
  atomicAdd(&(gainSumSquare[0]), float(gainSumSquareTemp));

}

__global__ void calcSampleGainSum(curandStateMtgp32* globalState,
				  const Mesh mesh, 
				  const unsigned* indicesOfPrisms, 
				  const double* importance,
				  const unsigned raysPerSample,
				  float *gainSum, 
				  float *gainSumSquare,
				  const unsigned sample_i,
				  const double* sigmaA, 
				  const double* sigmaE,
				  const unsigned lambdaResolution,
				  unsigned *globalOffsetMultiplicator
				  ) {

  int rayNumber = 0; 
  double gainSumTemp = 0;
  double gainSumSquareTemp = 0;
  Point samplePoint = mesh.getSamplePoint(sample_i);
  __shared__ unsigned blockOffset[4]; // 4 in case of warp-based raynumber
  
  // One thread can compute multiple rays
  while(true){
	// the whole block gets a new offset (==workload)
    rayNumber = getRayNumberBlockbased(blockOffset,raysPerSample,globalOffsetMultiplicator);
    if(rayNumber>=raysPerSample) break;

    // Get triangle/prism to start ray from
    unsigned startPrism             = indicesOfPrisms[rayNumber];
    unsigned startLevel             = startPrism/mesh.numberOfTriangles;
    unsigned startTriangle          = startPrism - (mesh.numberOfTriangles * startLevel);
    Point startPoint                = mesh.genRndPoint(startTriangle, startLevel, globalState);
    Ray ray                         = generateRay(startPoint, samplePoint);

	// get a random index in the wavelength array
    unsigned sigma_i                = genRndSigmas(lambdaResolution, globalState);
    assert(sigma_i < lambdaResolution);

	// calculate the gain for the whole ray at once
    double gain    = propagateRay(ray, &startLevel, &startTriangle, mesh, sigmaA[sigma_i], sigmaE[sigma_i]);
    gain          /= ray.length * ray.length; // important, since usually done in the reflection device function

	// include the stimulus from the starting prism and the importance of that ray
    gain          *= mesh.getBetaVolume(startPrism) * importance[startPrism];

    gainSumTemp       += gain;
    gainSumSquareTemp += gain * gain;

  }
  atomicAdd(&(gainSum[0]), float(gainSumTemp));
  atomicAdd(&(gainSumSquare[0]), float(gainSumSquareTemp));

}
