#include <stdio.h> /* printf */
#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */
#include <assert.h> /* assert */
#include <reflection.h> /* calcNextReflection */

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
  // Const data for reflection test
  const unsigned reflections = 1;         // Should be random generated (0 to X)
  const float reflectivity = 0.5;
  const float totalReflectionAngle = 45;
  int reflectionPlane = 1;                // Should be random generated (-1 or 1) // -1 = BOTTOM, 1 = TOP

  // Return if this wavelength block should be ignored
  if(wave_i == -1) return;

  // One thread can compute multiple rays
  // The current ray which we compute is based on the gid and an offset (number of threads*blocks)
  while ((rayNumber = gid + stride) < raysPerSample) {
    stride += blockDim.x * gridDim.x;
    // Get triangle/prism to start ray from
    unsigned startPrism = indicesOfPrisms[rayNumber + wave_i * raysPerSample];
    unsigned startLevel = startPrism/mesh.numberOfTriangles;
    unsigned startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
    Point startPoint = mesh.genRndPoint(startTriangle, startLevel, &(globalState[wave_i * gridDim.x]));

    // Calculate reflections as different ray propagations
    double gain = 1.0;
    for(unsigned reflection_i = 0; reflection_i < reflections; ++reflection_i){
      Point reflectionPoint = {0,0,0};
      float reflectionAngle = 0;
      Ray   reflectionRay   = {{0,0,0},{0,0,0}};

      // Calc reflectionPoint and reflectionAngle
      calcNextReflection(startPoint, samplePoint, (reflections-reflection_i), reflectionPlane, &reflectionPoint, &reflectionAngle, &mesh);

      // Propagate this part of the ray
      reflectionRay       = generateRay(startPoint, reflectionPoint);
      gain               *= propagateRay(reflectionRay, &startLevel, &startTriangle, &mesh, sigmaA[wave_i], sigmaE[wave_i]);
      if(reflectionAngle < totalReflectionAngle) 
	gain             *= reflectivity;
      startPoint          = reflectionPoint;
      reflectionPlane     = (reflectionPlane * -1);

    }
    // Calculate last part of ray without reflection
    Ray ray = generateRay(startPoint, samplePoint);
    gain   *= propagateRay(ray, &startLevel, &startTriangle, &mesh, sigmaA[wave_i], sigmaE[wave_i]);
    gain   *= mesh.getBetaValue(startPrism) * importance[startPrism + wave_i * mesh.numberOfPrisms];
    
    gainSum += gain;
    gainSumSquare += gain * gain;

  }
  atomicAdd(&(phiAse[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSum));
  atomicAdd(&(phiAseSquare[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSumSquare));

}
