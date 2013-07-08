#include <stdio.h>
#include "propagate_ray.h"
#include <curand_kernel.h>


__global__ void importanceKernel1(
    int point,
    int level,
    double *importance,
    double *points,
    double *xOfNormals,
    double *yOfNormals,
    int *positionsOfNormalVectors,
    int *neighbors,
    int *forbidden,
    double *betaValues,
    double *xOfTriangleCenter,
    double *yOfTriangleCenter,
    unsigned numberOfPoints,
    unsigned numberOfLevels,
    unsigned numberOfTriangles,
    float thicknessOfPrism,
    float sigmaA,
    float sigmaE,
    float nTot)
{
  double xPos, yPos, zPos;
  double prop;

  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  int startLevel = startPrism/numberOfTriangles;
  int startTriangle = startPrism - (numberOfTriangles * startLevel);

  if(startPrism >= numberOfTriangles * (numberOfLevels-1)){
    return;
  }

  xPos = points[point];
  yPos = points[point + numberOfPoints];
  zPos = level * thicknessOfPrism;

  // Calculate importance by propagation from trianglecenter to every other center
  // first FOR-loop
  prop = propagateRay(xOfTriangleCenter[startTriangle], yOfTriangleCenter[startTriangle], 
      thicknessOfPrism * (startLevel+0.5),  xPos, yPos, zPos, startTriangle, startLevel, 
      points, xOfNormals, yOfNormals, positionsOfNormalVectors, 
      neighbors, forbidden, betaValues,
      nTot, sigmaE, sigmaA, thicknessOfPrism, numberOfLevels, numberOfPoints, numberOfTriangles);

  importance[startPrism] = betaValues[startPrism] * prop;
}

__global__ void importanceKernel2(
    unsigned *numberOfImportantRays,
    double *importance,
    float sumPhi,
    int raysPerSample,
    unsigned numberOfPrisms){

  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= numberOfPrisms) return;
  // second FOR-loop
  // Calculate number of rays/prism
  numberOfImportantRays[startPrism] = (unsigned)(floor(importance[startPrism] / sumPhi * raysPerSample));
}


__global__ void importanceKernel3(
    int raysPerSample,
    int raysDump,
    unsigned *numberOfImportantRays,
    unsigned numberOfLevels,
    unsigned numberOfTriangles){
  // TODO What happens with random failure ?
  // Distribute the remaining rays randomly
  // OPTIMIZE: parallelize
  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= numberOfTriangles * (numberOfLevels-1)) return;
  if(startPrism == 0){
    int raysLeft = raysPerSample - raysDump;
    curandState randomState;
    curand_init(1234,0,0,&randomState);
    for (int i_r=0; i_r < raysLeft; i_r++){
      int rand_t = (int) ceil((curand_uniform(&randomState) * numberOfTriangles)) - 1;
      int rand_z = (int) ceil((curand_uniform(&randomState) * (numberOfLevels - 1))) - 1;
      numberOfImportantRays[rand_t + rand_z * numberOfTriangles]++;
    }
  }
}

__global__ void importanceKernel4(
    unsigned *numberOfImportantRays,
    double *importance,
    float *surface,
    int surfaceTotal,
    unsigned raysPerSample,
    unsigned numberOfPrisms,
    unsigned numberOfTriangles){

  int startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  int startLevel = startPrism/numberOfTriangles;
  int startTriangle = startPrism - (numberOfTriangles * startLevel);
  if(startPrism >= numberOfPrisms) return;
  //  Now think about the mount of rays which would come out of this volume(surface)
  //  dividing this number with the new amount of rays gives the final importance weight for this area!
  if (numberOfImportantRays[startPrism] > 0){
    importance[startPrism] = raysPerSample * surface[startTriangle] / (surfaceTotal * numberOfImportantRays[startPrism]);
  }
  else{
    importance[startPrism] = 0; 
  }
  return;
}
