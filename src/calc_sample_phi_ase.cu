#include <curand_kernel.h> /* curand_uniform */
#include <stdio.h> /* printf */
#include "propagate_ray.h"

/**
 * Does the raytracing for a single Sample point (in a defined level).
 * This Kernel has to be started for each sample point with the same value for iterations
 * and the same number of blocks/threads.
 *
 * \var globalState the state of the mersenneTwister PRNG
 * 		(has a maximum of 200 positions!)
 * \var phi points to a memory region which is initialized with 0
 * 		(can hold one value for each sample point)
 * \var point2D the index of the current sample point (points to p_in)
 * \var level the level of the current sample point (how deep we are through the material)
 * 		(always for the same combination of startprism+samplepoint
 */
__global__ void calcSamplePhiAse(
		curandStateMtgp32* globalState,
		float* phiASE,
		int point2D,
		int level,
		double *points,
		double *xOfNormals,
		double *yOfNormals,
		int *positionsOfNormalVectors,
		int *neighbors,
		int *forbidden,
		int* triangleIndices,
		double* betaValues,
		double* importance,
		unsigned* indicesOfPrisms,
		unsigned raysPerSample,
		double nTot,
		double sigmaE,
		double sigmaA,
		double thicknessOfPrism,
		int numberOfLevels,
		int numberOfPoints,
		int numberOfTriangles
		){

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  double endPointX = points[point2D];
  double endPointY = points[numberOfPoints + point2D];
  double endPointZ = level * thicknessOfPrism;
  __shared__ double threadGain[256]; //MUST be the same as number of threads
  threadGain[threadIdx.x] = 0.;

  // on thread can compute multiple rays
  for (int i=0; ; ++i){

	  // the current ray which we compute is based on the id and an offset (number of threads*blocks)
	  int rayNumber = id + (blockDim.x*gridDim.x * i);
	  if(rayNumber >= raysPerSample){
		  break;
	  }

	  // get a new prism to start from
	  int startPrism = indicesOfPrisms[rayNumber];
	  int startLevel = startPrism/numberOfTriangles;
	  int startTriangle = startPrism - (numberOfTriangles * startLevel);

	  // Get triangle vertex indicies
	  int t1 = triangleIndices[startTriangle];
	  int t2 = triangleIndices[startTriangle + numberOfTriangles];
	  int t3 = triangleIndices[startTriangle + 2 * numberOfTriangles];

	  // random startpoint generation
	  double u = curand_uniform(&globalState[blockIdx.x]);
	  double v = curand_uniform(&globalState[blockIdx.x]);

	  if((u+v)>1)
	  {
		  u = 1-u;
		  v = 1-v;
	  }
	  double w = 1-u-v;

	  // convert the random startpoint into coordinates
	  double xRand = (points[t1] * u) + (points[t2] * v) + (points[t3] * w);
	  double yRand = (points[numberOfPoints + t1] * u) + (points[numberOfPoints + t2] * v) + (points[numberOfPoints + t3] * w);
	  double zRand = (startLevel + curand_uniform(&globalState[blockIdx.x])) * thicknessOfPrism;

	  // propagate the ray
	  double gain = propagateRay(xRand, yRand, zRand, endPointX, endPointY, endPointZ, 
				   startTriangle, startLevel, points, xOfNormals, yOfNormals, 
				   positionsOfNormalVectors, neighbors, forbidden,  betaValues,
				   nTot, sigmaE, sigmaA, thicknessOfPrism, numberOfLevels, numberOfPoints, numberOfTriangles);

	  threadGain[threadIdx.x] += gain * betaValues[startPrism] * importance[startPrism];
  }

  // reduce the shared memory to one element (CUDA by Example, Chapter 5.3)
  __syncthreads();
  unsigned i = blockDim.x/2;
  while(i != 0){
	  if(threadIdx.x < i){
		  threadGain[threadIdx.x] += threadGain[threadIdx.x + i];
	  }
	  __syncthreads();
	  i /= 2;
  }

  // thread 0 writes it to the global memory
  if(threadIdx.x == 0){
	  atomicAdd(&(phiASE[point2D + (level * numberOfPoints)]), float(threadGain[threadIdx.x]));
  }
}
