//#include "calc_sample_phi_ase.h"
#include <stdio.h> /* printf */
#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */
#include <assert.h> /* assert */

__global__ void calcSamplePhiAse(curandStateMtgp32* globalState,
				 Mesh mesh, 
				 const unsigned* indicesOfPrisms, 
				 const int* indicesOfWavelengths, 
				 const unsigned* numberOfReflections,
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
  unsigned wavelengthOffset = wave_i * mesh.numberOfPrisms;

  // Return if this wavelength block should be ignored
  if(wave_i == -1) return;

  // One thread can compute multiple rays
  // The current ray which we compute is based on the gid and an offset (number of threads*blocks)
  while ((rayNumber = gid + stride) < raysPerSample) {
    stride += blockDim.x * gridDim.x;
    // Get triangle/prism to start ray from
    unsigned startPrism       = indicesOfPrisms[rayNumber + wave_i * raysPerSample];
    unsigned reflection_i     = numberOfReflections[rayNumber + wave_i * raysPerSample];
    unsigned reflections      = (reflection_i + 1) / 2;
    int reflectionPlane       = (reflection_i % 2 == 0) ? -1 : 1;
    unsigned startLevel       = startPrism/mesh.numberOfTriangles;
    unsigned startTriangle    = startPrism - (mesh.numberOfTriangles * startLevel);
    unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms * blockDim.y;

    Point startPoint = mesh.genRndPoint(startTriangle, startLevel, &(globalState[wave_i * gridDim.x]));

    // Calculate reflections as different ray propagations
    double gain    = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, &mesh, sigmaA[wave_i], sigmaE[wave_i]);
    gain          *= mesh.getBetaValue(startPrism) * importance[startPrism + wavelengthOffset + reflectionOffset];

    gainSum       += gain;
    gainSumSquare += gain * gain;

  }
  atomicAdd(&(phiAse[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSum));
  atomicAdd(&(phiAseSquare[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSumSquare));

}
