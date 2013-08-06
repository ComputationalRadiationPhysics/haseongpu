
#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */


__global__ void calcSamplePhiAse(
		curandStateMtgp32* globalState,
		Mesh mesh, 
		const unsigned* indicesOfPrisms, 
		const double* importance,
		const unsigned raysPerSample, 
		float *phiAse, 
		const unsigned sample_i,
		double *sigmaA, 
		double *sigmaE
		) {

  // Get global ID
  int gid = threadIdx.x + blockIdx.x * blockDim.x;
  int rayNumber = 0;
  unsigned stride = 0;
  unsigned wavelengthOffset = blockIdx.y;

  extern __shared__ double threadGain[]; // Size is set by Kernelparameter
  threadGain[threadIdx.x] = 0.;
  Point samplePoint = mesh.getSamplePoint(sample_i);

  // One thread can compute multiple rays
  // The current ray which we compute is based on the gid and an offset (number of threads*blocks)
  while ((rayNumber = gid + stride) < raysPerSample) {
          stride += blockDim.x * gridDim.x;
  	  // Get triangle prism to start from
  	  int startPrism = indicesOfPrisms[rayNumber + wavelengthOffset * raysPerSample];
  	  int startLevel = startPrism/mesh.numberOfTriangles;
  	  int startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);

	  Point startPoint = mesh.genRndPoint(startTriangle, startLevel, &(globalState[wavelengthOffset * gridDim.x]));
	  Ray ray          = generateRay(startPoint, samplePoint);
	  double gain      = propagateRay(ray, startLevel, startTriangle, &mesh, sigmaA[wavelengthOffset], sigmaE[wavelengthOffset]);

	  gain *= mesh.getBetaValue(startPrism);
	  gain *= importance[startPrism + wavelengthOffset * mesh.numberOfPrisms];

	  threadGain[threadIdx.x] += gain;
  }

  // Reduce the threadGain array (CUDA by Example, Chapter 5.3)  
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
	  atomicAdd(&(phiAse[sample_i  + wavelengthOffset * mesh.numberOfSamples]), float(threadGain[threadIdx.x]));
  }
}
