#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */


__global__ void calcSamplePhiAse(curandStateMtgp32* globalState,
				 Mesh mesh, 
				 const unsigned* indicesOfPrisms, 
				 const int* indicesOfWavelengths, 
				 const double* importance,
				 const unsigned raysPerSample, 
				 float *phiAse, 
				 float *phiAseSquare,
				 const unsigned sample_i,
				 double *sigmaA, 
				 double *sigmaE
				 ) {

  int wave_i = indicesOfWavelengths[blockIdx.y];
  int gid = threadIdx.x + blockIdx.x * blockDim.x;
  int rayNumber = 0;
  unsigned stride = 0;
  double gainSum = 0;
  double gainSumSquare = 0;
  Point samplePoint = mesh.getSamplePoint(sample_i);

  if(wave_i == -1) return;

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
  atomicAdd(&(phiAse[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSum));
  atomicAdd(&(phiAseSquare[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSumSquare));


}
