/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @licence GPLv3
 *
 */

#ifndef THREAD_H
#define THREAD_H

#include <mesh.h>

#include <pthread.h> /* pthread_t */
#include <vector> /* std::vector */

/**
 * @brief Wrapper for calcPhiAse on pthread base.
 *        This function will spawn a thread for
 *        each function call and start calcPhiAse.
 *
 * @param minRaysPerSample Lower bound for raysPerSample
 *                         in case of adaptive sampling.
 * @param maxRaysPerSample Uppper boud for raysPerSample
 *                         in case of adaptive sampling.
 * @param maxRepetitions   Number of Repetitions will
 *                         be done, when not reaching mse threshold
 * @param dMesh            All information about triangles, points, contants. 
 *                         Is located in device memory. See mesh.h for details.
 * @param hMesh            Same as dMesh, but locatet in host memory.
 * @param sigmaA           Vector with Absorption values
 * @param sigmaE           Vector with Emission values
 * @param mseThreshold     Threshold for adaptive and repetitive sampling.
 *                         Not reaching this threshold leads to recomputations.
 * @param useReflections   Rays can reflect on upper and lower surface of gain medium
 * @param phiAse           Reference to phiAse result (one value for every sample point).
 * @param mse              Reference to mse result (one value for every sample point).
 * @param totalRays        Reference to numberOfRays simulated per sample point.
 * @param gpu_i            Number of device that should be used.
 * @param minSample_i      Smallest Index of sample point to calculate.
 * @param maxSample_i      Biggest Index of sample point to calculate.
 * @param runtime          Reference to the needed runtime.
 *
 * @deprecated will be completly replaced by mpi
 *             or should be replaced by c++11 threads
 * @return     threadId
 */
pthread_t calcPhiAseThreaded( unsigned &minRaysPerSample,
			      const unsigned maxRaysPerSample,
			      const unsigned maxRepetitions,
			      const Mesh& mesh,
			      const std::vector<double>& sigmaA,
			      const std::vector<double>& sigmaE,
			      const double mseThreshold,
			      const bool useReflections,
			      std::vector<float> &phiAse,
			      std::vector<double> &mse,
			      std::vector<unsigned> &totalRays,
			      unsigned gpu_i,
			      unsigned minSample_i,
			      unsigned maxSample_i,
			      float &runtime);
/**
 * @brief Wait for all threads to finish
 *
 */
void joinAll(std::vector<pthread_t> threadIds);

#endif /* THREAD_H */
