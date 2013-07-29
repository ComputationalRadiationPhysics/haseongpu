
#include <mesh.h>
#include <geometry.h> /* generateRay */
#include <propagate_ray.h> /* propagateRay */

/**
 * @brief Generates a random point inside the given prism (triangle and level).
 *        For random number generation, the mersenne twister (high periodicity)
 *        is used.
 * 
 * @param globalState State for random number generation (mersenne twister).
 *                    The state need to be initialized before. See
 *                    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MTGP/
 *                    for more information.
 *
 * @return point is a random point inside a prism
 **/
// __device__ Point calcRndStartPoint(Triangle triangle, unsigned level, double thickness, curandStateMtgp32* globalState){
//   Point startPoint = {0,0,0};
//   double u = curand_uniform(&globalState[blockIdx.x]);
//   double v = curand_uniform(&globalState[blockIdx.x]);

//   if((u+v)>1)
//     {
//       u = 1-u;
//       v = 1-v;
//     }
//   double w = 1-u-v;

//   // convert the random startpoint into coordinates
//   startPoint.x = (triangle.A.x * u) + (triangle.B.x * v) + (triangle.C.x * w);
//   startPoint.y = (triangle.A.y * u) + (triangle.B.y * v) + (triangle.C.y * w);
//   startPoint.z = (level + curand_uniform(&globalState[blockIdx.x])) * thickness;

//   return startPoint;
// }


__global__ void calcSamplePhiAse(
		curandStateMtgp32* globalState,
		Mesh mesh, 
		const unsigned* indicesOfPrisms, 
		const double* importance,
		const unsigned raysPerSample, 
		float *phiAse, 
		const unsigned sample_i,
		const double sigmaA, 
		const double sigmaE, 
		const double nTot) {

  // Get global ID
  int gid = threadIdx.x + blockIdx.x * blockDim.x;
  int rayNumber = 0;
  unsigned stride = 0;

  __shared__ double threadGain[256]; //MUST be the same as number of threads
  threadGain[threadIdx.x] = 0.;
  const Point samplePoint = mesh.getSamplePoint(sample_i);

  // One thread can compute multiple rays
  // The current ray which we compute is based on the gid and an offset (number of threads*blocks)
  while ((rayNumber = gid + stride) < raysPerSample) {
          stride += blockDim.x * gridDim.x;
  	  // Get triangle prism to start from
  	  int startPrism = indicesOfPrisms[rayNumber];
  	  int startLevel = startPrism/mesh.numberOfTriangles;
  	  int startTriangle = startPrism - (mesh.numberOfTriangles * startLevel);

	  Point startPoint = mesh.genRndPoint(startTriangle, startLevel, globalState);
	  Ray ray          = generateRay(startPoint, samplePoint);
	  double gain      = propagateRay(ray, startLevel, startTriangle, &mesh, sigmaA, sigmaE, nTot, mesh.thickness );

	  gain *= mesh.getBetaValue(startPrism);
	  gain *= importance[startPrism];

	  threadGain[threadIdx.x] += gain;
  }

  // Reduce the threadGain array (CUDA by Example, Chapter 5.3)  
  __syncthreads();
  unsigned i = blockDim.x/2;
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
