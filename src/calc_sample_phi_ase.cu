#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */

#define BLOCKDIM 256

__global__ void calcSamplePhiAse(
				 curandStateMtgp32* globalState,
				 Mesh mesh, 
				 const unsigned* indicesOfPrisms, 
				 const double* importance,
				 const unsigned raysPerSample, 
				 float *phiAse, 
				 float *phiAseSquare,
				 const unsigned sample_i,
				 double *sigmaA, 
				 double *sigmaE
				 ) {

  // Get global ID
  int gid = threadIdx.x + blockIdx.x * blockDim.x;
  int rayNumber = 0;
  unsigned stride = 0;
  unsigned wave_i = blockIdx.y;
  double gainSum = 0;
  double gainSumSquare = 0;

  __shared__ double threadGain[BLOCKDIM];
  __shared__ double threadGainSquare[BLOCKDIM];
  threadGain[threadIdx.x] = 0.;
  threadGainSquare[threadIdx.x] = 0.;
  Point samplePoint = mesh.getSamplePoint(sample_i);

  // One thread can compute multiple rays
  // The current ray which we compute is based on the gid and an offset (number of threads*blocks)
  while ((rayNumber = gid + stride) < raysPerSample) {
          stride += blockDim.x * gridDim.x;
  	  // Get triangle prism to start from
  	  int startPrism = indicesOfPrisms[rayNumber + wave_i * raysPerSample];
  	  int startLevel = startPrism/mesh.numberOfTriangles;
  	  int startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);

	  Point startPoint = mesh.genRndPoint(startTriangle, startLevel, &(globalState[wave_i * gridDim.x]));
	  Ray ray          = generateRay(startPoint, samplePoint);
	  double gain      = propagateRay(ray, startLevel, startTriangle, &mesh, sigmaA[wave_i], sigmaE[wave_i]);

	  gain *= mesh.getBetaValue(startPrism);
	  gain *= importance[startPrism + wave_i * mesh.numberOfPrisms];

	  gainSum += gain;
	  gainSumSquare += gain * gain;

  }

  threadGain[threadIdx.x] = gainSum;
  threadGainSquare[threadIdx.x] = gainSumSquare;
  
  // Reduce the threadGain array (CUDA by Example, Chapter 5.3)  
  __syncthreads();
  unsigned i = blockDim.x/2;
  while(i != 0){
	  if(threadIdx.x < i){
		  threadGain[threadIdx.x] += threadGain[threadIdx.x + i];
		  threadGainSquare[threadIdx.x] += threadGainSquare[threadIdx.x + i];
	  }
	  __syncthreads();
	  i /= 2;
  }
  // thread 0 writes it to the global memory
  if(threadIdx.x == 0){
    atomicAdd(&(phiAse[sample_i  + wave_i * mesh.numberOfSamples]), float(threadGain[threadIdx.x]));
    atomicAdd(&(phiAseSquare[sample_i  + wave_i * mesh.numberOfSamples]), float(threadGain[threadIdx.x]));
  }
}
