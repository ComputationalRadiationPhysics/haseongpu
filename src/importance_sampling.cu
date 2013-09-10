#include "importance_sampling.h"
#include <mesh.h>
#include <stdio.h>
#include <propagate_ray.h>
#include <geometry.h>
#include <assert.h>
#include <curand_kernel.h>
#include <cudachecks.h>

/**
 * @brief calculates a first estimate on the importance of each prism, based on a single ray started in the center of each prism
 *
 * @param *importance will contain the initial importance for each prism
 *
 * @param *sumPhi will contain the cumulative sum of the importance values
 *
 * For other parameters, see documentation of importanceSampling()
 *
 */
__global__ void propagateFromTriangleCenter(
    Mesh mesh,
    double *importance,
    float *sumPhi,
    const unsigned sample_i,
    const double sigmaA,
    const double sigmaE
    ){

  double gain = 0;
  unsigned reflection_i = blockIdx.z;
  unsigned reflections = (reflection_i + 1) / 2;
  int reflectionPlane  = (reflection_i % 2 == 0)? -1 : 1;

  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= mesh.numberOfPrisms){
    return;
  }
  unsigned startLevel = startPrism/(mesh.numberOfTriangles);
  unsigned startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
  Point startPoint = mesh.getCenterPoint(startTriangle, startLevel);
  Point samplePoint = mesh.getSamplePoint(sample_i);
  unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;

  gain = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, &mesh, sigmaA, sigmaE);

  // DEBUG
  //printf("C x:%f y:%f z:%f reflections: %u reflectionPlane %d gain: %f\n", samplePoint.x, samplePoint.y, samplePoint.z, reflections, reflectionPlane, gain);

  importance[startPrism + reflectionOffset] = mesh.getBetaValue(startPrism) * gain;
  
  atomicAdd(sumPhi, float(importance[startPrism + reflectionOffset]));

}

/**
 * @brief uses a given importance distribution to decide how many rays will be launched from each prism
 *
 * @param *raysDump will contain the number of rays which were mapped to a specific prism
 * 
 * for other parameters, see documentation of importanceSampling()
 */
__global__ void distributeRaysByImportance(Mesh mesh,
					   unsigned *raysPerPrism,
					   double *importance,
					   float *sumPhi,
					   unsigned raysPerSample,
					   unsigned *raysDump){

  unsigned reflection_i = blockIdx.z;
  unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;

  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= mesh.numberOfPrisms) return;
  raysPerPrism[startPrism + reflectionOffset] = (unsigned) floor(importance[startPrism + reflectionOffset] / (*sumPhi) * raysPerSample);

  atomicAdd(raysDump, raysPerPrism[startPrism + reflectionOffset]);

}

/**
 * @brief takes a number of rays and distributes them randomly over the available prisms
 *
 * @param *raysPerPrism the number of rays for each prism (will be changed)
 *
 * @param *raysDump the number of rays which were already distributed
 *
 * for other parameters, see documentation of importanceSampling()
 *
 */
__global__ void distributeRemainingRaysRandomly(Mesh mesh,
						unsigned *raysPerPrism,
						unsigned raysPerSample,
						unsigned *raysDump){
  
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  int raysLeft = raysPerSample - (*raysDump);

  if(id < raysLeft){
    curandState randomState;
    curand_init(id,0,0,&randomState);
    int rand_t = (int ) ceil(curand_uniform(&randomState) * mesh.numberOfTriangles) - 1;
    int rand_z = (int ) ceil(curand_uniform(&randomState) * (mesh.numberOfLevels-1)) - 1;
    unsigned randomPrism = rand_t + rand_z * mesh.numberOfTriangles;
    atomicAdd(&(raysPerPrism[randomPrism]),1);
  } 

}


/**
 * @brief corrects the importance to match with the randomly distributed rays
 *
 * @param *raysPerPrism the number of rays to be launced for each prism
 *
 * @param *importance the importance for each prism (will be changed)
 *
 * for other parameters, see documentation of importanceSampling()
 */
__global__ void recalculateImportance(Mesh mesh,
				      unsigned *raysPerPrism,
				      unsigned raysPerSample,
				      double *importance){


  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned reflection_i = blockIdx.z;
  unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;


  if(startPrism >= mesh.numberOfPrisms){
    return;
  }
  int startLevel = startPrism/(mesh.numberOfTriangles);
  int startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
  if(raysPerPrism[startPrism + reflectionOffset] > 0){
    importance[startPrism + reflectionOffset] = raysPerSample * mesh.surfaces[startTriangle] / (mesh.surfaceTotal * raysPerPrism[startPrism + reflectionOffset]);
  }
  else{
    importance[startPrism + reflectionOffset] = 0;
  }
}


// unused, because we didn't find a good way to parallelize it...
// OPTIMIZE
/**
 * @brief maps every ray to a specific prism
 *
 * @param *raysPerPrism the number of rays to launch in each prism
 *
 * @param raysPerSample the total number of rays to launch 
 *
 * @param *indicesOfPrisms a mapping for each ray to a specific prism
 *
 */
__global__ void createCumulativeSum1(
    Mesh mesh,
    unsigned *raysPerPrism,
    unsigned *cumulativeSum){

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  if(id==0){
    cumulativeSum[0] = 0;
  }
  if(id < mesh.numberOfPrisms-1){
    cumulativeSum[id+1] = raysPerPrism[id];
  }
}

__global__ void createCumulativeSum2(
    Mesh mesh,
    unsigned *cumulativeSum){

  for(int i=0;i<mesh.numberOfPrisms;i++){
    cumulativeSum[i+1] += cumulativeSum[i];
  }
  //printf("PartialSum sum: %d\n",partialSums[0]);
}

__global__ void mapRaysToPrism(
    Mesh mesh,
    unsigned *raysPerPrism,
    unsigned raysPerSample,
    unsigned *indicesOfPrisms,
    unsigned *cumulativeSum){

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  if(id >= mesh.numberOfPrisms) return;

  unsigned absoluteRay = cumulativeSum[id];
  for(unsigned prism_i=cumulativeSum[id]; prism_i < indicesOfPrisms[id]; ++prism_i){
    for(unsigned ray_i=0; ray_i < raysPerPrism[prism_i]; ++ray_i){
      indicesOfPrisms[absoluteRay++] = prism_i;
    }
  }
}


unsigned importanceSampling(unsigned sample_i,
			    const unsigned reflectionSlices,
			    Mesh deviceMesh,
			    const unsigned raysPerSample,
			    const double sigmaA,
			    const double sigmaE,
			    double *importance,
			    unsigned *raysPerPrism,
			    const bool distributeRandomly,
			    dim3 blockDim,
			    dim3 gridDim){

  float *hostSumPhi = (float*) malloc(sizeof(float));

  float *sumPhi;
  unsigned *raysDump;
  unsigned hostRaysDump;

  CUDA_CHECK_RETURN(cudaMalloc(&sumPhi, sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&raysDump, sizeof(unsigned)));

  *hostSumPhi = 0.f;
  hostRaysDump = 0;

  CUDA_CHECK_RETURN(cudaMemcpy(sumPhi,hostSumPhi, sizeof(float),cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(raysDump,&hostRaysDump, sizeof(unsigned),cudaMemcpyHostToDevice));

  dim3 gridDimReflection(gridDim.x, 1, reflectionSlices);
  CUDA_CHECK_KERNEL_SYNC(propagateFromTriangleCenter<<< gridDimReflection, blockDim >>>(deviceMesh, importance, sumPhi, sample_i, sigmaA, sigmaE));
  CUDA_CHECK_KERNEL_SYNC(distributeRaysByImportance<<< gridDimReflection, blockDim >>>(deviceMesh, raysPerPrism,importance, sumPhi, raysPerSample, raysDump));

  // Distribute remaining rays randomly if wanted
  if(distributeRandomly){
    CUDA_CHECK_KERNEL_SYNC(distributeRemainingRaysRandomly<<< 200,blockDim >>>(deviceMesh ,raysPerPrism, raysPerSample, raysDump));
    hostRaysDump = raysPerSample;
  }
  else {
    CUDA_CHECK_RETURN(cudaMemcpy(&hostRaysDump, raysDump,  sizeof(unsigned),cudaMemcpyDeviceToHost));
  }

  CUDA_CHECK_KERNEL_SYNC(recalculateImportance<<< gridDimReflection, blockDim >>>(deviceMesh, raysPerPrism, raysPerSample, importance));

  free(hostSumPhi);
  cudaFree(sumPhi);
  cudaFree(raysDump);

  return raysPerSample;
}
