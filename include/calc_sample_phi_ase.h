#ifndef calc_sample_phi_ase_H
#define calc_sample_phi_ase_H

#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>
#include <mesh.h>

/**
 * @brief Calculates the phi ASE value for the given 
 *        samplepoint. This is done by a monte carlo
 *        simulation with random generated rays with
 *        the direction to the samplepoint.
 *
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 * @param globalState State for random number generation (mersenne twister).
 *                    The state need to be initialized before. See
 *                    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MTGP/
 *                    for more information.
 *
 * @param samplePoint the point for which the phi ASE value is calculated
 *
 * @param mesh all information about prism structure. See mesh.h for details.
 *
 * @param indecesOfPrisms tells each thread(ray) on which prism it should start
 *
 * @param importance constant values, calculated by importance sampling.
 *                   Reduces calculation time for this simulation.
 *
 * @param raysPerSample 
 *
 * @param phiAse Solution vector
 *
 * @param sample_i index of samplepoint in the mesh->samples array
 *
 * @param sigmaA
 *
 * @param sigmaE
 *
 * @param nTot
 *
 **/
__global__ void calcSamplePhiAse(curandStateMtgp32* globalState, 
				 Mesh mesh,
				 const unsigned* indicesOfPrisms, 
				 const double* importance,
				 const unsigned raysPerSample, 
				 float *phiAse,
				 const unsigned sample_i,
				 double *sigmaA,
				 double *sigmaE,
				 const double nTot);

#endif /* calc_sample_phi_ase_H */
