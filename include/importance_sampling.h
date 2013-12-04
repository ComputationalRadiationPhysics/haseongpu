#ifndef importance_sampling_H
#define importance_sampling_H

#include <mesh.h>
#include <thrust/device_vector.h>

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
 * @param raysDump the rays which were mapped to a specific prism (pointer to memory location)
 *
 * @return the number of rays which are used for one sample point
 */
float importanceSamplingPropagation(
    unsigned sample_i,
    const unsigned reflectionSlices,
    Mesh deviceMesh,
    const unsigned numberOfPrisms,
    const double sigmaA,
    const double sigmaE,
    thrust::device_vector<double> &importance
    );

unsigned importanceSamplingDistribution(
    const unsigned reflectionSlices,
    Mesh deviceMesh,
    const unsigned numberOfPrisms,
    const unsigned raysPerSample,
    thrust::device_vector<double> &importance,
    thrust::device_vector<unsigned> &raysPerPrism,
    const float hSumPhi,
    const bool distributeRandomly
    );

#endif /* importance_sampling_H */
