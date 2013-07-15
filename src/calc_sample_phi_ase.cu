#include <curand_kernel.h> /* curand_uniform */
#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */


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


__global__ void calcSamplePhiAse(
		curandStateMtgp32* globalState,
		Point samplePoint,
		Mesh mesh, 
		unsigned* indicesOfPrisms, 
		double* importance,
		unsigned raysPerSample, 
		float *phiAse, 
		const unsigned sample_i,
		const double sigmaA, 
		const double sigmaE, 
		const double nTot) {

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  Triangle *triangles = mesh.triangles;
  unsigned numberOfTriangles = mesh.numberOfTriangles;

  __shared__ double threadGain[256]; //MUST be the same as number of threads
  threadGain[threadIdx.x] = 0.;

  // One thread can compute multiple rays
  int rayNumber;
  unsigned i=0;
  // the current ray which we compute is based on the id and an offset (number of threads*blocks)
  while ((rayNumber = id + (blockDim.x*gridDim.x * i++)) < raysPerSample) {

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
	  double gain = propagateRayNew(ray, startLevel, startTriangle, sigmaA, sigmaE, nTot, mesh.thickness );

	  gain *= startTriangle.betaValues[startLevel];
	  gain *= importance[startPrism];

	  threadGain[threadIdx.x] += gain;
  }

  // reduce the shared memory to one element (CUDA by Example, Chapter 5.3)
  __syncthreads();
  
  i = blockDim.x/2;
  while(i != 0){
	  if(threadIdx.x < i){
		  threadGain[threadIdx.x] += threadGain[threadIdx.x + i];
	  }
	  __syncthreads();
	  i /= 2;
  }
  // thread 0 writes it to the global memory
  if(threadIdx.x == 0){
	  atomicAdd(&(phiAse[sample_i]), float(threadGain[threadIdx.x]));
  }
}
