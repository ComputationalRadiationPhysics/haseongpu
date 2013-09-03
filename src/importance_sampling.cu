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
    unsigned sample_i,
    unsigned reflection_i,
    unsigned reflections,
    int reflectionPlane,
    double *sigmaA,
    double *sigmaE
    ){

  double gain = 0;

  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= mesh.numberOfPrisms){
    return;
  }
  unsigned startLevel = startPrism/(mesh.numberOfTriangles);
  unsigned startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
  Point startPoint = mesh.getCenterPoint(startTriangle, startLevel);
  Point samplePoint = mesh.getSamplePoint(sample_i);
  unsigned wavelengthOffset = blockIdx.y * mesh.numberOfPrisms;
  unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms * blockDim.y;

  gain = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, &mesh, sigmaA[blockIdx.y], sigmaE[blockIdx.y]);

  if(reflections > 0){
    assert(gain == 0);
  }

  // DEBUG
  //printf("C x:%f y:%f z:%f reflections: %u reflectionPlane %d gain: %f\n", samplePoint.x, samplePoint.y, samplePoint.z, reflections, reflectionPlane, gain);

  importance[startPrism + wavelengthOffset + reflectionOffset] = mesh.getBetaValue(startPrism) * gain;
  
  atomicAdd(&(sumPhi[blockIdx.y]), float(importance[startPrism + wavelengthOffset + reflectionOffset]));

}

/**
 * @brief uses a given importance distribution to decide how many rays will be launched from each prism
 *
 * @param *raysDump will contain the number of rays which were mapped to a specific prism
 * 
 * for other parameters, see documentation of importanceSampling()
 */
__global__ void distributeRaysByImportance(
    Mesh mesh,
    unsigned reflection_i,
    unsigned *raysPerPrism,
    double *importance,
    float *sumPhi,
    unsigned raysPerSample,
    unsigned *raysDump){

  unsigned wavelengthOffset = blockIdx.y * mesh.numberOfPrisms;
  unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms * blockDim.y;



  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= mesh.numberOfPrisms) return;
  raysPerPrism[startPrism + wavelengthOffset + reflectionOffset] = (unsigned) floor(importance[startPrism + wavelengthOffset + reflectionOffset] / (sumPhi[blockIdx.y]) * raysPerSample);

  if(reflection_i > 0){
    assert(raysPerPrism[startPrism + wavelengthOffset + reflectionOffset] == 0);
  }

  atomicAdd(&(raysDump[blockIdx.y]), raysPerPrism[startPrism + wavelengthOffset + reflectionOffset]);

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
__global__ void distributeRemainingRaysRandomly(
    Mesh mesh,
    unsigned reflection_i,
    unsigned *raysPerPrism,
    unsigned raysPerSample,
    unsigned *raysDump){

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  int raysLeft = raysPerSample-raysDump[blockIdx.y];
  unsigned wavelengthOffset = blockIdx.y * mesh.numberOfPrisms;
  //unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms * blockDim.y;

  if(id < raysLeft){
    curandState randomState;
    curand_init(id,0,0,&randomState);
    int rand_t = (int ) ceil(curand_uniform(&randomState) * mesh.numberOfTriangles) - 1;
    int rand_z = (int ) ceil(curand_uniform(&randomState) * (mesh.numberOfLevels-1)) - 1;
    unsigned randomPrism = rand_t + rand_z * mesh.numberOfTriangles;
    //atomicAdd(&raysPerPrism[randomPrism + wavelengthOffset + reflectionOffset],1);
    atomicAdd(&raysPerPrism[randomPrism + wavelengthOffset],1);
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
__global__ void recalculateImportance(
    Mesh mesh,
    unsigned reflection_i,
    unsigned *raysPerPrism,
    unsigned raysPerSample,
    double *importance){
  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned wavelengthOffset = blockIdx.y * mesh.numberOfPrisms;
  unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms * blockDim.y;

  if(startPrism >= mesh.numberOfPrisms){
    return;
  }
  int startLevel = startPrism/(mesh.numberOfTriangles);
  int startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
  if(raysPerPrism[startPrism + wavelengthOffset + reflectionOffset] > 0){
    importance[startPrism + wavelengthOffset + reflectionOffset] = raysPerSample * mesh.surfaces[startTriangle] / (mesh.surfaceTotal * raysPerPrism[startPrism + wavelengthOffset + reflectionOffset]);
  }
  else{
    importance[startPrism + wavelengthOffset + reflectionOffset] = 0;
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
    unsigned *cumulativeSums){

  unsigned wavelengthOffset = blockIdx.y * mesh.numberOfPrisms;
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  if(id==0){
    cumulativeSums[0 + wavelengthOffset] = 0;
  }
  if(id < mesh.numberOfPrisms-1){
    cumulativeSums[id+1 + wavelengthOffset] = raysPerPrism[id + wavelengthOffset];
  }
}

__global__ void createCumulativeSum2(
    Mesh mesh,
    unsigned *cumulativeSums){

  unsigned wavelengthOffset = blockIdx.y * mesh.numberOfPrisms;
  for(int i=0;i<mesh.numberOfPrisms;i++){
    cumulativeSums[i+1 + wavelengthOffset] += cumulativeSums[i + wavelengthOffset];
  }
  //printf("PartialSum sum: %d\n",partialSums[0]);
}

__global__ void mapRaysToPrism(
    Mesh mesh,
    unsigned *raysPerPrism,
    unsigned raysPerSample,
    unsigned *indicesOfPrisms,
    unsigned *cumulativeSums){

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  if(id >= mesh.numberOfPrisms) return;

  unsigned wavelengthOffset = blockIdx.y;
  unsigned absoluteRay = cumulativeSums[id + wavelengthOffset * mesh.numberOfPrisms];
  for(unsigned prism_i=cumulativeSums[id + wavelengthOffset * mesh.numberOfPrisms]; prism_i < indicesOfPrisms[id + wavelengthOffset * raysPerSample]; ++prism_i){
    for(unsigned ray_i=0; ray_i < raysPerPrism[prism_i + wavelengthOffset * mesh.numberOfPrisms]; ++ray_i){
      indicesOfPrisms[absoluteRay++] = prism_i;
    }
  }
}


unsigned importanceSampling(
    unsigned sample_i,
    unsigned reflectionSlices,
    Mesh deviceMesh,
    unsigned raysPerSample,
    double *sigmaA,
    double *sigmaE,
    double *importance,
    unsigned *raysPerPrism,
    dim3 blockDim,
    dim3 gridDim){

  float *hostSumPhi = (float*) malloc(gridDim.y * sizeof(float));
  unsigned *hostDumpHost = (unsigned*) malloc(gridDim.y * sizeof(unsigned));
  float *sumPhi;
  unsigned *raysDump;

  CUDA_CHECK_RETURN(cudaMalloc(&sumPhi, gridDim.y * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&raysDump, gridDim.y * sizeof(unsigned)));

  for(unsigned i=0; i < gridDim.y; ++i){
    hostSumPhi[i] = 0.f;
    hostDumpHost[i] = 0;
  }

  CUDA_CHECK_RETURN(cudaMemcpy(sumPhi,hostSumPhi, gridDim.y * sizeof(float),cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(raysDump,hostDumpHost, gridDim.y * sizeof(unsigned),cudaMemcpyHostToDevice));

  for(unsigned reflection_i = 0; reflection_i < reflectionSlices; ++reflection_i){
    int reflectionPlane  = (reflection_i % 2 == 0)? -1 : 1;
    unsigned reflections = (reflection_i + 1) / 2;
    CUDA_CHECK_KERNEL_SYNC(propagateFromTriangleCenter<<< gridDim, blockDim >>>(deviceMesh,importance,sumPhi,sample_i, reflection_i, reflections, reflectionPlane,sigmaA, sigmaE));
  }

  for(unsigned reflection_i = 0; reflection_i < reflectionSlices; ++reflection_i){
    CUDA_CHECK_KERNEL_SYNC(distributeRaysByImportance<<< gridDim, blockDim >>>(deviceMesh, reflection_i, raysPerPrism,importance,sumPhi,raysPerSample,raysDump));
  }

  //for(unsigned reflection_i = 0; reflection_i < reflectionSlices; ++reflection_i){
    CUDA_CHECK_KERNEL_SYNC(distributeRemainingRaysRandomly<<< gridDim,blockDim >>>(deviceMesh, 0, raysPerPrism,raysPerSample,raysDump));
    //}

  for(unsigned reflection_i = 0; reflection_i < reflectionSlices; ++reflection_i){
    CUDA_CHECK_KERNEL_SYNC(recalculateImportance<<< gridDim, blockDim >>>(deviceMesh,reflection_i, raysPerPrism,raysPerSample,importance));
  }

  free(hostSumPhi);
  free(hostDumpHost);
  cudaFree(sumPhi);
  cudaFree(raysDump);

  return raysPerSample;
}
