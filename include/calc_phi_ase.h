#ifndef calc_phi_ase_H
#define calc_phi_ase_H
#include <vector>
#include <mesh.h>

/**
 * @brief Calculates the integral over time of ASE for each 
 *        point in the given mesh. The calculation is
 *        accelerated by CUDA (compute capability 2.0).
 *
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 *
 * @licence GPLv3
 *
 * @param &threads Number of threads used by kernel.
 *                 By reference because it can be used by calling
 *                 function for statistics etc.
 *
 * @param &blocks Number of blocks used by kernel is 
 *                maximal 200 because of mersenne twister.
 *                By reference because it can be used by calling
 *                function for statistics etc.
 *
 * @param &hostRaysPerSample  Number of rays propagated for each sample point.
 *                            By reference because it can be used by calling
 *                            function for statistics etc.
 *
 * @param mesh All information about triangles, points, contants. 
 *             See mesh.h for details.
 *
 * @param sigmaA
 *
 * @param sigmaE
 *
 * @param dndtAse Solution vector with dndt ASE values
 *
 **/
float calcPhiAse ( const unsigned hostRaysPerSample,
		   const unsigned maxRaysPerSample,
		   const unsigned maxRepetitions,
		   const Mesh& mesh,
		   const Mesh& hostMesh,
		   const std::vector<double>& sigmaA,
		   const std::vector<double>& sigmaE,
		   const double mseThreshold,
		   const bool useReflections,
		   std::vector<float> &phiAse,
		   std::vector<double> &mse,
		   std::vector<unsigned> &totalRays,
		   const unsigned gpu_i,
		   const unsigned minSample_i,
		   const unsigned maxSample_i,
		   float &runtime);


#endif
