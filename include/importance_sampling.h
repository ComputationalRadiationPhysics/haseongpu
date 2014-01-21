#ifndef importance_sampling_H
#define importance_sampling_H

#include <mesh.h>

/**
 * @brief calculates an importance sampling for the given mesh and ray parameters
 *
 * @param samplePoint the sample point for which the importance is destined
 *
 * @param deviceMesh the device mesh
 *
 * @param raysPerSample the number of rays which will be used for the sample point
 *
 * @param sigmaA
 *
 * @param sigmaE
 *
 * @param *importance the importance values on the GPU (pointer to the memory location)
 *
 * @param *sumPhi cumulative values of all the gains from propagation (pointer to the memory location)
 *
 * @param *raysPerPrism for each prism, how many rays will be launced
 *
 * @param *indicesOfPrisms for each ray, from which prism will it start (currently unused)
 *
 * @param raysDump the rays which were mapped to a specific prism (pointer to memory location)
 *
 * @param threads the number of threads for the kernels
 *
 * @param blocks the number of blocks for the kernels
 *
 * @return the number of rays which are used for one sample point
 */
//unsigned importanceSampling(
//			    unsigned sample_i,
//			    unsigned maxReflections,
//			    Mesh deviceMesh,
//			    unsigned raysPerSample, 
//			    double sigmaA, 
//			    double sigmaE, 
//			    double *importance, 
//			    unsigned *raysPerPrism,
//			    bool distributeRandomly,
//			    dim3 threads,
//			    dim3 blocks);
float importanceSamplingPropagation(unsigned sample_i,
			    const unsigned reflectionSlices,
			    Mesh deviceMesh,
			    const double sigmaA,
			    const double sigmaE,
			    double *importance,
			    dim3 blockDim,
			    dim3 gridDim);
unsigned importanceSamplingDistribution(
			    const unsigned reflectionSlices,
			    Mesh deviceMesh,
			    const unsigned raysPerSample,
			    double *importance,
			    unsigned *raysPerPrism,
				float hSumPhi,
			    const bool distributeRandomly,
			    dim3 blockDim,
			    dim3 gridDim);

#endif /* importance_sampling_H */
