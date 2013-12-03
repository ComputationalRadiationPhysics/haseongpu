#include <importance_sampling.h>
#include <mesh.h>
#include <stdio.h>
#include <propagate_ray.h>
#include <geometry.h>
#include <assert.h>
#include <curand_kernel.h>
#include <cudachecks.h>
#include <cuda_utils.h>
#include <reflection.h> /* ReflectionPlane */
#include <thrust/device_vector.h>

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
__global__ void propagateFromTriangleCenter(const Mesh mesh,
					    double *importance,
					    float *sumPhi,
					    const unsigned sample_i,
					    const double sigmaA,
					    const double sigmaE
					    ){

  double gain = 0;
  unsigned reflection_i = blockIdx.z;
  unsigned reflections = (reflection_i + 1) / 2;
  ReflectionPlane reflectionPlane  = (reflection_i % 2 == 0)? BOTTOM_REFLECTION : TOP_REFLECTION;

  unsigned startPrism = threadIdx.x + blockIdx.x * blockDim.x;
  if(startPrism >= mesh.numberOfPrisms){
    return;
  }
  unsigned startLevel = startPrism/(mesh.numberOfTriangles);
  unsigned startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
  Point startPoint = mesh.getCenterPoint(startTriangle, startLevel);
  Point samplePoint = mesh.getSamplePoint(sample_i);
  unsigned reflectionOffset = reflection_i * mesh.numberOfPrisms;

  gain = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, mesh, sigmaA, sigmaE); 
  importance[startPrism + reflectionOffset] = mesh.getBetaValue(startPrism) * gain;
  if(mesh.getBetaValue(startPrism) < 0 || gain < 0 || importance[startPrism+reflectionOffset] < 0){
    printf("beta: %f importance: %f gain: %f\n", mesh.getBetaValue(startPrism), importance[startPrism + reflectionOffset], gain);
  }


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
  if(raysPerPrism[startPrism + reflectionOffset] > raysPerSample){
	  printf("importance: %f sumPhi: %f raysPerPrism[%d]: %d (max %d)\n",importance[startPrism+reflectionOffset],*sumPhi,startPrism+reflectionOffset,raysPerPrism[startPrism+reflectionOffset],raysPerSample);
  }
  assert(raysPerPrism[startPrism + reflectionOffset] <= raysPerSample);
  atomicAdd(raysDump, raysPerPrism[startPrism + reflectionOffset]);

}

/**
 * @brief takes a number of rays and distributes them randomly over the available prisms
 *        Warning: Does not distribute to reflection slices !!!
 *
 * @param *raysPerPrism the number of rays for each prism (will be changed)
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

float importanceSamplingPropagation(unsigned sample_i,
			    const unsigned reflectionSlices,
			    Mesh deviceMesh,
				const unsigned numberOfPrisms,
			    const double sigmaA,
			    const double sigmaE,
			    thrust::device_vector<double>& importance
				){

  thrust::device_vector<float> dSumPhi(1,0);
  int pBlock=128;
  dim3 gridDimReflection(ceil(float(numberOfPrisms)/pBlock), 1, reflectionSlices);
  CUDA_CHECK_KERNEL_SYNC(propagateFromTriangleCenter<<< gridDimReflection, pBlock >>>(
			  deviceMesh, 
			  thrust::raw_pointer_cast(&importance[0]), 
			  thrust::raw_pointer_cast(&dSumPhi[0]), 
			  sample_i, 
			  sigmaA, 
			  sigmaE
			  ));

  return dSumPhi[0];
}

unsigned importanceSamplingDistribution(
			    const unsigned reflectionSlices,
			    Mesh deviceMesh,
				const unsigned numberOfPrisms,
			    const unsigned raysPerSample,
			    thrust::device_vector<double> &importance,
			    thrust::device_vector<unsigned> &raysPerPrism,
				const float hSumPhi,
			    const bool distributeRandomly
			    ){

  thrust::device_vector<unsigned> dRaysDump(1,0);
  thrust::device_vector<float> dSumPhi(1,hSumPhi);

  int dBlock=288;
  CUDA_CHECK_KERNEL_SYNC(distributeRaysByImportance<<< dim3(ceil(float(numberOfPrisms)/dBlock)), dBlock >>>(
			  deviceMesh, 
			  thrust::raw_pointer_cast(&raysPerPrism[0]),
			  thrust::raw_pointer_cast(&importance[0]), 
			  thrust::raw_pointer_cast(&dSumPhi[0]),
			  raysPerSample, 
			  thrust::raw_pointer_cast(&dRaysDump[0])
			  ));

  // Distribute remaining rays randomly if wanted
  if(distributeRandomly){
	int rBlock = 256;
	int rGrid = ceil(float(raysPerSample-dRaysDump[0])/rBlock);
    CUDA_CHECK_KERNEL_SYNC(distributeRemainingRaysRandomly<<< rGrid,rBlock >>>(
				deviceMesh,
				thrust::raw_pointer_cast(&raysPerPrism[0]),
				raysPerSample,
				thrust::raw_pointer_cast(&dRaysDump[0])
				));
	dRaysDump[0] = raysPerSample;
  }

  int iBlock=256;
  CUDA_CHECK_KERNEL_SYNC(recalculateImportance<<< dim3(ceil(float(numberOfPrisms/iBlock)),1,reflectionSlices), iBlock >>>(
			  deviceMesh, 
			  thrust::raw_pointer_cast(&raysPerPrism[0]), 
			  dRaysDump[0],
			  thrust::raw_pointer_cast(&importance[0])
			  ));

  return dRaysDump[0];
}
