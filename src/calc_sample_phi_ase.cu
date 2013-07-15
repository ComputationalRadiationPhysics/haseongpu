#include <curand_kernel.h> /* curand_uniform */
#include <stdio.h> /* printf */
#include <mesh.h>
#include <geometry.h>
#include <propagate_ray.h> /* propagateRay */

#define ID 2560

// ##############################################################
// # Reconstruction                                             #
// ##############################################################
__device__ Point calcRndStartPoint(Triangle triangle, unsigned level, double thickness, curandStateMtgp32* globalState){
  Point startPoint = {0,0,0};
  double u = curand_uniform(&globalState[blockIdx.x]);
  double v = curand_uniform(&globalState[blockIdx.x]);

  if((u+v)>1)
    {
      u = 1-u;
      v = 1-v;
    }
  double w = 1-u-v;

  // convert the random startpoint into coordinates
  startPoint.x = (triangle.A.x * u) + (triangle.B.x * v) + (triangle.C.x * w);
  startPoint.y = (triangle.A.y * u) + (triangle.B.y * v) + (triangle.C.y * w);
  startPoint.z = (level + curand_uniform(&globalState[blockIdx.x])) * thickness;

  return startPoint;
}


__global__ void calcSamplePhiAseNew(curandStateMtgp32* globalState, Point samplePoint, Mesh mesh, unsigned* indicesOfPrisms, double* importance,unsigned raysPerSample, float *phiAse, const unsigned sample_i, const double sigmaA, const double sigmaE, const double nTot) {
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  Triangle *triangles = mesh.triangles;
  unsigned numberOfTriangles = mesh.numberOfTriangles;

  // One thread can compute multiple rays
  int rayNumber;
  unsigned i=0;
  // the current ray which we compute is based on the id and an offset (number of threads*blocks)
  while ((rayNumber = id + (blockDim.x*gridDim.x * i++)) < raysPerSample) {
  	  // TODO check indices on new structs
  	  // Get triangle prism to start from
  	  int startPrism = indicesOfPrisms[rayNumber];
  	  int startLevel = startPrism/numberOfTriangles;
  	  int startTriangle_i = startPrism - (numberOfTriangles * startLevel);
  	  Triangle startTriangle = triangles[startTriangle_i];

  	  // Random startpoint generation
	  Point startPoint = calcRndStartPoint(startTriangle, startLevel, mesh.thickness, globalState);

	  // Ray generation
	  Ray ray = generateRay(startPoint, samplePoint);

  	  // // propagate the ray
	  double gain = propagateRay(ray, startLevel, startTriangle, sigmaA, sigmaE, nTot, mesh.thickness );

	  gain *= startTriangle.betaValues[startLevel];
	  gain *= importance[startPrism];
	  atomicAdd(&phiAse[sample_i], float(gain));

  }

}


