//#include "calc_sample_phi_ase.h"
#include <stdio.h> /* printf */
#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */
#include <assert.h> /* assert */
#include <reflection.h> /* ReflectionPlane */

__global__ void calcSamplePhiAse(curandStateMtgp32* globalState,
				 const Mesh mesh, 
				 const unsigned* indicesOfPrisms, 
				 const unsigned wave_i, 
				 const unsigned* numberOfReflections,
				 const double* importance,
				 const unsigned raysPerSample,
				 float *phiAse, 
				 float *phiAseSquare,
				 const unsigned sample_i,
				 double sigmaA, 
				 double sigmaE
				 ) {

  int gid = threadIdx.x + blockIdx.x * blockDim.x;
  int rayNumber = 0;
  unsigned stride = 0;
  double gainSum = 0;
  double gainSumSquare = 0;
  Point samplePoint = mesh.getSamplePoint(sample_i);

  // One thread can compute multiple rays
  // The current ray which we compute is based on the gid and an offset (number of threads*blocks)
  while ((rayNumber = gid + stride) < raysPerSample) {
    stride += blockDim.x * gridDim.x;
    // Get triangle/prism to start ray from
    unsigned startPrism             = indicesOfPrisms[rayNumber];
    unsigned reflection_i           = numberOfReflections[rayNumber]; 
    unsigned reflections            = (reflection_i + 1) / 2;
    ReflectionPlane reflectionPlane = (reflection_i % 2 == 0) ? BOTTOM_REFLECTION : TOP_REFLECTION;
    unsigned startLevel             = startPrism/mesh.numberOfTriangles;
    unsigned startTriangle          = startPrism - (mesh.numberOfTriangles * startLevel);
    unsigned reflectionOffset       = reflection_i * mesh.numberOfPrisms;

    Point startPoint = mesh.genRndPoint(startTriangle, startLevel, globalState);

    // Calculate reflections as different ray propagations
    double gain    = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, mesh, sigmaA, sigmaE);
    gain          *= mesh.getBetaValue(startPrism) * importance[startPrism + reflectionOffset];

    gainSum       += gain;
    gainSumSquare += gain * gain;

  }
  atomicAdd(&(phiAse[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSum));
  atomicAdd(&(phiAseSquare[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSumSquare));

}

__global__ void calcSamplePhiAseWithoutReflections(curandStateMtgp32* globalState,
				 const Mesh mesh, 
				 const unsigned* indicesOfPrisms, 
				 const unsigned wave_i, 
				 const double* importance,
				 const unsigned raysPerSample,
				 float *phiAse, 
				 float *phiAseSquare,
				 const unsigned sample_i,
				 double sigmaA, 
				 double sigmaE
				 ) {

  int gid = threadIdx.x + blockIdx.x * blockDim.x;
  int rayNumber = 0;
  unsigned stride = 0;
  double gainSum = 0;
  double gainSumSquare = 0;
  Point samplePoint = mesh.getSamplePoint(sample_i);

  // One thread can compute multiple rays
  // The current ray which we compute is based on the gid and an offset (number of threads*blocks)
  while ((rayNumber = gid + stride) < raysPerSample) {
    stride += blockDim.x * gridDim.x;
    // Get triangle/prism to start ray from
    unsigned startPrism             = indicesOfPrisms[rayNumber];
    unsigned startLevel             = startPrism/mesh.numberOfTriangles;
    unsigned startTriangle          = startPrism - (mesh.numberOfTriangles * startLevel);

    Point startPoint = mesh.genRndPoint(startTriangle, startLevel, globalState);
    Ray ray   = generateRay(startPoint, samplePoint);

    double gain    = propagateRay(ray, &startLevel, &startTriangle, mesh, sigmaA, sigmaE);

    gain          /= ray.length * ray.length; // important, since usually done in the reflection device function!
    gain          *= mesh.getBetaValue(startPrism) * importance[startPrism];

    gainSum       += gain;
    gainSumSquare += gain * gain;

  }
  atomicAdd(&(phiAse[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSum));
  atomicAdd(&(phiAseSquare[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSumSquare));
  // if(phiAseSquare[sample_i  + wave_i * mesh.numberOfSamples] > 100000000){
  //   printf("phiAse² > 100000000: ",phiAseSquare[sample_i  + wave_i * mesh.numberOfSamples]);
  // }

}
