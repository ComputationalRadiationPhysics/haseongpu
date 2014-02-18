/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @licence GPLv3
 *
 */

#ifndef calcPhiAseMPI_H
#define calcPhiAseMPI_H

#include <vector>
#include <mesh.h>

/**
 * @brief A wrapper for calcPhiAse, that distributes sample points
 *        to the available MPI nodes.The Nodes will split 
 *        up in one head node and the others as compute nodes. 
 *        The head node distributes the available sample
 *        points by demand.
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
 * @param hsigmaA           Vector with Absorption values
 * @param hsigmaE           Vector with Emission values
 * @param mseThreshold     Threshold for adaptive and repetitive sampling.
 *                         Not reaching this threshold leads to recomputations.
 * @param useReflections   Rays can reflect on upper and lower surface of gain medium
 * @param phiAse           Reference to phiAse result (one value for every sample point).
 * @param mse              Reference to mse result (one value for every sample point).
 * @param totalRays        Reference to numberOfRays simulated per sample point.
 * @param gpu_i            Number of device that should be used.
 *
 * @return number of used compute nodes
 */
float calcPhiAseMPI ( unsigned &minRaysPerSample,
		      const unsigned maxRaysPerSample,
		      const unsigned maxRepetitions,
		      const Mesh& dMesh,
		      const Mesh& hMesh,
		      const std::vector<double>& hSigmaA,
		      const std::vector<double>& hSigmaE,
		      const double mseThreshold,
		      const bool useReflections,
		      std::vector<float> &hPhiAse,
		      std::vector<double> &hMse,
		      std::vector<unsigned> &hTotalRays,
		      unsigned gpu_i);

#endif /* calcPhiAseMPI_H */
