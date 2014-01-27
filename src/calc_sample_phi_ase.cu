//#include "calc_sample_phi_ase.h"
#include <stdio.h> /* printf */
#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */
#include <assert.h> /* assert */
#include <reflection.h> /* ReflectionPlane */



__device__ unsigned getRayNumberWarpbased(unsigned* blockOffset,unsigned raysPerSample, unsigned *globalOffsetMultiplicator){
	if((threadIdx.x &31) == 0){
		blockOffset[(threadIdx.x>>5)] = atomicInc(globalOffsetMultiplicator,raysPerSample);
	}
	__syncthreads();

	return (threadIdx.x &31) + (blockOffset[(threadIdx.x>>5)] <<5) ;

}
__device__ unsigned getRayNumberBlockbased(unsigned* blockOffset,unsigned raysPerSample,unsigned *globalOffsetMultiplicator){
	if(threadIdx.x == 0){
		blockOffset[0] = atomicInc(globalOffsetMultiplicator,raysPerSample);
	}
	__syncthreads();

	return threadIdx.x + (blockOffset[0] <<7) ;

}

__global__ void calcSampleGainSumWithReflection(curandStateMtgp32* globalState,
				 const Mesh mesh, 
				 const unsigned* indicesOfPrisms, 
				 const unsigned wave_i, 
				 const unsigned* numberOfReflectionSlices,
				 const double* importance,
				 const unsigned raysPerSample,
				 float *gainSum, 
				 float *gainSumSquare,
				 const unsigned sample_i,
				 const double sigmaA, 
				 const double sigmaE,
				 unsigned *globalOffsetMultiplicator
				 ) {

  int rayNumber = 0;
  double gainSumTemp = 0;
  double gainSumSquareTemp = 0;
  Point samplePoint = mesh.getSamplePoint(sample_i);
  __shared__ unsigned blockOffset[4];

  // One thread can compute multiple rays
  while (true) {
	rayNumber = getRayNumberBlockbased(blockOffset,raysPerSample,globalOffsetMultiplicator);
	if(rayNumber >= raysPerSample) break;

    // Get triangle/prism to start ray from
    unsigned startPrism             = indicesOfPrisms[rayNumber];
    unsigned reflection_i           = numberOfReflectionSlices[rayNumber]; //numberOfReflectio == ReflectionSlice
    unsigned reflections            = (reflection_i + 1) / 2;
    ReflectionPlane reflectionPlane = (reflection_i % 2 == 0) ? BOTTOM_REFLECTION : TOP_REFLECTION;
    unsigned startLevel             = startPrism / mesh.numberOfTriangles;
    unsigned startTriangle          = startPrism - (mesh.numberOfTriangles * startLevel);
    unsigned reflectionOffset       = reflection_i * mesh.numberOfPrisms;

    //Point startPoint = mesh.getCenterPoint(startTriangle, startLevel);
    Point startPoint = mesh.genRndPoint(startTriangle, startLevel, globalState);

    // Calculate reflections as different ray propagations
    double gain    = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, mesh, sigmaA, sigmaE);
    gain          *= mesh.getBetaValue(startPrism) * importance[startPrism + reflectionOffset];
    
    assert(!isnan(mesh.getBetaValue(startPrism)));
    assert(!isnan(importance[startPrism + reflectionOffset]));
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
				 const unsigned wave_i, 
				 const double* importance,
				 const unsigned raysPerSample,
				 float *gainSum, 
				 float *gainSumSquare,
				 const unsigned sample_i,
				 const double sigmaA, 
				 const double sigmaE,
				 unsigned *globalOffsetMultiplicator
				 ) {

  int rayNumber = 0; 
  double gainSumTemp = 0;
  double gainSumSquareTemp = 0;
  Point samplePoint = mesh.getSamplePoint(sample_i);
  __shared__ unsigned blockOffset[4];

  // One thread can compute multiple rays
  // The current ray which we compute is based on the gid and an offset (number of threads*blocks)
  while(true){
	  rayNumber = getRayNumberBlockbased(blockOffset,raysPerSample,globalOffsetMultiplicator);
	  if(rayNumber>=raysPerSample) break;

	  // Get triangle/prism to start ray from
	  unsigned startPrism             = indicesOfPrisms[rayNumber];
	  unsigned startLevel             = startPrism/mesh.numberOfTriangles;
	  unsigned startTriangle          = startPrism - (mesh.numberOfTriangles * startLevel);

	  Point startPoint = mesh.genRndPoint(startTriangle, startLevel, globalState);
	  Ray ray   = generateRay(startPoint, samplePoint);

	  double gain    = propagateRay(ray, &startLevel, &startTriangle, mesh, sigmaA, sigmaE);

	  gain          /= ray.length * ray.length; // important, since usually done in the reflection device function!
	  gain          *= mesh.getBetaValue(startPrism) * importance[startPrism];

	  gainSumTemp       += gain;
	  gainSumSquareTemp += gain * gain;

  }
  atomicAdd(&(gainSum[0]), float(gainSumTemp));
  atomicAdd(&(gainSumSquare[0]), float(gainSumSquareTemp));

}
