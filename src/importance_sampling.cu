#include <mesh.h>
#include <stdio.h>
#include <propagate_ray.h>
#include <geometry.h>
#include <assert.h>
#include <curand_kernel.h>
#include <cudachecks.h>

__global__ void importanceSamplingKernel1(
    Mesh mesh,
    double *importance,
    float *sumPhi,
    Point samplePoint,
    double sigmaA, 
    double sigmaE, 
    double nTot){

  Triangle *triangles = mesh.triangles;
  __shared__ double threadPhi[256];
  double gain = 0;
  Ray ray;
  Point startPoint;
  Triangle startTriangle;

  threadPhi[threadIdx.x] = 0;

  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= mesh.numberOfPrisms){
    return;
  }
  int level_i = startPrism/(mesh.numberOfTriangles);
  int triangle_i = startPrism - (mesh.numberOfTriangles * level_i);

  startTriangle = triangles[triangle_i];
  startPoint.x = startTriangle.center.x;
  startPoint.y = startTriangle.center.y;
  startPoint.z = (level_i + 0.5) * mesh.thickness;

  ray = generateRay(startPoint, samplePoint);
  gain = propagateRayNew(ray, level_i, startTriangle, triangles, sigmaA, sigmaE, nTot, mesh.thickness);
  importance[startPrism] = startTriangle.betaValues[level_i] * gain;

  threadPhi[threadIdx.x] = importance[triangle_i + level_i * mesh.numberOfTriangles];
  __syncthreads();

  unsigned i = blockDim.x/2;
  while(i != 0){
    if(threadIdx.x < i){
      threadPhi[threadIdx.x] += threadPhi[threadIdx.x + i];
    }
    __syncthreads();
    i /= 2;
  }
  if(threadIdx.x == 0){
    atomicAdd(sumPhi, float(threadPhi[threadIdx.x]));
  }
}

__global__ void importanceSamplingKernel2(
    Mesh mesh,
    unsigned *raysPerPrism,
    double *importance,
    float *sumPhi,
    unsigned raysPerSample,
    unsigned *raysDump){
  __shared__ unsigned raySum[256];
  raySum[threadIdx.x] = 0;
  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= mesh.numberOfPrisms) return;
  raysPerPrism[startPrism] = (unsigned) floor(importance[startPrism] / (*sumPhi) * raysPerSample);
  raySum[threadIdx.x] = raysPerPrism[startPrism];
  __syncthreads();

  unsigned i = blockDim.x/2;
  while(i != 0){
    if(threadIdx.x < i){
      raySum[threadIdx.x] += raySum[threadIdx.x + i];
    }
    __syncthreads();
    i /= 2;
  }
  if(threadIdx.x == 0){
    atomicAdd(raysDump, raySum[threadIdx.x]);
  }
}

__global__ void importanceSamplingKernel3(
    Mesh mesh,
    unsigned *raysPerPrism,
    unsigned raysPerSample,
    unsigned *raysDump){

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  int raysLeft = raysPerSample-(*raysDump);

  if(id < raysLeft){
    curandState randomState;
    curand_init(id,0,0,&randomState);
    int rand_t = (int ) ceil(curand_uniform(&randomState) * mesh.numberOfTriangles) - 1;
    int rand_z = (int ) ceil(curand_uniform(&randomState) * (mesh.numberOfLevels-1)) - 1;
    atomicAdd(&raysPerPrism[rand_t + rand_z * mesh.numberOfTriangles],1);
  }
}

__global__ void importanceSamplingKernel4(
    Mesh mesh,
    unsigned *raysPerPrism,
    unsigned raysPerSample,
    double *importance){ 
  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= mesh.numberOfPrisms){
    return;
  }
  int startLevel = startPrism/(mesh.numberOfTriangles);
  int startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
  if(raysPerPrism[startPrism] > 0){
    importance[startPrism] = raysPerSample * mesh.triangles[startTriangle].surface / (mesh.surface * raysPerPrism[startPrism]);
  }else{
    importance[startPrism] = 0;
  }
}

// unused, because we didn't find a good way to parallelize it...
// OPTIMIZE
// TODO
__global__ void importanceSamplingKernel5(
    Mesh mesh,
    unsigned *raysPerPrism,
    unsigned raysPerSample,
    unsigned *indicesOfPrisms){

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  if(id==0){
    // Prism scheduling for gpu threads
    unsigned absoluteRay = 0;
    for(unsigned prism_i=0; prism_i < mesh.numberOfPrisms; ++prism_i){
      for(unsigned ray_i=0; ray_i < raysPerPrism[prism_i]; ++ray_i){
        indicesOfPrisms[absoluteRay++] = prism_i;
#if TEST_VALUES==true
        assert(absoluteRay <= raysPerSample);
#endif
      }
    }
  }
}

unsigned importanceSamplingGPU(
    Point samplePoint, 
    Mesh deviceMesh,
    unsigned raysPerSample, 
    double sigmaA, 
    double sigmaE, 
    double nTot,  
    double *importance, 
    float *sumPhi,
    unsigned *raysPerPrism,
    unsigned *indicesOfPrisms,
    unsigned *raysDump,
    int threads,
    int blocks){

  float *sumPhiHost = (float*) malloc(sizeof(float));
  unsigned *raysDumpHost = (unsigned*) malloc(sizeof(unsigned));

  *sumPhiHost = 0.f;
  *raysDumpHost = 0;

  CUDA_CHECK_RETURN(cudaMemcpy(sumPhi,sumPhiHost,sizeof(float),cudaMemcpyHostToDevice));
  CUDA_CHECK_RETURN(cudaMemcpy(raysDump,raysDumpHost,sizeof(unsigned),cudaMemcpyHostToDevice));

  importanceSamplingKernel1<<< blocks,threads >>>(deviceMesh,importance,sumPhi,samplePoint,sigmaA,sigmaE,nTot);
  importanceSamplingKernel2<<< blocks,threads >>>(deviceMesh,raysPerPrism,importance,sumPhi,raysPerSample,raysDump);
  importanceSamplingKernel3<<< blocks,threads >>>(deviceMesh,raysPerPrism,raysPerSample,raysDump);
  importanceSamplingKernel4<<< blocks,threads >>>(deviceMesh,raysPerPrism,raysPerSample,importance);

//  CUDA_CHECK_RETURN(cudaMemcpy(hostRaysPerPrism,raysPerPrism, hostMesh.numberOfPrisms*sizeof(unsigned),cudaMemcpyDeviceToHost));
//
//    // Prism scheduling for gpu threads
//  for(unsigned prism_i=0, absoluteRay = 0; prism_i < hostMesh.numberOfPrisms; ++prism_i){
//    for(unsigned ray_i=0; ray_i < hostRaysPerPrism[prism_i]; ++ray_i){
//      hostIndicesOfPrisms[absoluteRay++] = prism_i;
//      assert(absoluteRay <= hostRaysPerSample);
//    }
//  }
//  // Copy dynamic sample data to device
//  CUDA_CHECK_RETURN(cudaMemcpy(indicesOfPrisms, hostIndicesOfPrisms, hostRaysPerSample * sizeof(unsigned), cudaMemcpyHostToDevice));

  //return raysDump;  //use this, if we don't use randomly distributed rays (Kernel3)
  return raysPerSample;
}
