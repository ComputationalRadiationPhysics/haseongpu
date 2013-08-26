#include <stdio.h> /* printf */
#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */
#include <assert.h> /* assert */



__device__ double calcIntersectionAngle(const Ray ray, float *reflectionAngle){
  // Calc intesection angle with z-plane
  float nominator = abs(ray.dir.z);
  float denominator = sqrt((ray.dir.x * ray.dir.x) + (ray.dir.y * ray.dir.y) + (ray.dir.z * ray.dir.z));
  if(denominator != 0.0){
    double radian = asin(nominator / denominator);
    *reflectionAngle = (180 / 3.1415926) * radian;
    return 0;
  }
  return 1;
}

__device__ int calcPlaneIntersectionPoint(const Ray reflectionRay, const int reflectionPlane, const Mesh *mesh, Point *intersectionPoint){
  // Assume that mesh is on x/y axis and parallel to x/y axis
  double planeZ = 0.0;
  if(reflectionPlane == 1){
    planeZ = mesh->thickness * mesh->numberOfLevels;
  }
  double denominator = reflectionRay.dir.z;  
  if (denominator != 0.0){
    double nominator = planeZ - reflectionRay.p.z;
    double length = nominator/denominator;
    if(length > 0){
      intersectionPoint->x = reflectionRay.p.x + length * reflectionRay.dir.x;
      intersectionPoint->y = reflectionRay.p.y + length * reflectionRay.dir.y;
      intersectionPoint->z = reflectionRay.p.z + length * reflectionRay.dir.z;
      return 0;
    }
  }
  return 1;
}


// BUG its not that easy, crystal has to be mirrored
__device__ Ray generateReflectionRay(const Point startPoint, Point endPoint, const int reflectionCount, const int reflectionPlane, const Mesh *mesh){

  if(reflectionPlane % 2 == 0){
    // Even reflectionCounts are postponements
    endPoint.z = endPoint.z + reflectionPlane * (reflectionCount * mesh->thickness * mesh->numberOfLevels); 
  }
  else {
    // Odd reflectionsCounts are reflects
    int mirrorPlaneZ = reflectionPlane * (reflectionCount - 1) * mesh->thickness * mesh->numberOfLevels;
    endPoint.z = mirrorPlaneZ + reflectionPlane * (mirrorPlaneZ - reflectionPlane * endPoint.z);
    
  }

  return generateRay(startPoint, endPoint);

}

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
  const unsigned reflections = 2;
  const float reflectivity = 0.5;
  const float totalReflectionAngle = 45;
  int reflectionPlane = 1; // -1 = BOTTOM, 1 = TOP

  if(wave_i == -1) return;

  // One thread can compute multiple rays
  // The current ray which we compute is based on the gid and an offset (number of threads*blocks)
  while ((rayNumber = gid + stride) < raysPerSample) {
    stride += blockDim.x * gridDim.x;
    // Get triangle prism to start from
    unsigned startPrism = indicesOfPrisms[rayNumber + wave_i * raysPerSample];
    unsigned startLevel = startPrism/mesh.numberOfTriangles;
    unsigned startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);
    Point startPoint = mesh.genRndPoint(startTriangle, startLevel, &(globalState[wave_i * gridDim.x]));

    // Calculate reflections as different ray propagations
    double gain = 1.0;
    for(unsigned reflection_i = 0; reflection_i < reflections; ++reflection_i){
      Point intersectionPoint = {0,0,0};
      float reflectionAngel = 0;

      Ray reflectionRay = generateReflectionRay(startPoint, samplePoint, reflections - reflection_i, reflectionPlane, &mesh);
      calcPlaneIntersectionPoint(reflectionRay, reflectionPlane, &mesh, &intersectionPoint);
      calcIntersectionAngle(reflectionRay, &reflectionAngel);
    
      //Debug output
      // BUG should a ray always have the same angle ???
      // --> wrong generateReflectionRay function
      if(gid == 0) printf("[%d][%d] StartPoint %f, %f, %f EndPoint %f, %f, %f Angle %f IntersectionPoint %f, %f, %f reflectionPlane %d\n", sample_i, reflection_i, startPoint.x, startPoint.y, startPoint.z, samplePoint.x, samplePoint.y, samplePoint.z, reflectionAngel, intersectionPoint.x, intersectionPoint.y, intersectionPoint.z, reflectionPlane);
      
      reflectionRay   = generateRay(startPoint, intersectionPoint);
      gain           *= propagateRay(reflectionRay, &startLevel, &startTriangle, &mesh, sigmaA[wave_i], sigmaE[wave_i]);
      startPoint      = intersectionPoint;
      reflectionPlane = (reflectionPlane * -1);

    }
    gain *= mesh.getBetaValue(startPrism);
    gain *= importance[startPrism + wave_i * mesh.numberOfPrisms];
    
    gainSum += gain;
    gainSumSquare += gain * gain;

  }
  atomicAdd(&(phiAse[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSum));
  atomicAdd(&(phiAseSquare[sample_i  + wave_i * mesh.numberOfSamples]), float(gainSumSquare));


}
