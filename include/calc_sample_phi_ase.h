#ifndef calc_sample_phi_ase_H
#define calc_sample_phi_ase_H

#include <curand_mtgp32dc_p_11213.h>
#include <mesh.h>

/**
 * @brief Calculates the gain sum for the given 
 *        sample point with or without reflections. This is done by a monte carlo
 *        simulation with randomly generated rays and
 *        wavelenghts.
 *
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 **/

__global__ void calcSampleGainSumWithReflection(curandStateMtgp32* globalState, 
				 const Mesh mesh,
				 const unsigned* indicesOfPrisms, 
				 const unsigned* numberOfReflections,
				 const double* importance,
				 const unsigned raysPerSample,
				 float *gainSum,
				 float *gainSumSquare,
				 const unsigned sample_i,
				 const double *sigmaA,
				 const double *sigmaE,
         const unsigned maxInterpolation,
				 unsigned *globalOffsetMultiplicator
				 );

__global__ void calcSampleGainSum(curandStateMtgp32* globalState, 
				 const Mesh mesh,
				 const unsigned* indicesOfPrisms, 
				 const double* importance,
				 const unsigned raysPerSample,
				 float *gainSum,
				 float *gainSumSquare,
				 const unsigned sample_i,
				 const double *sigmaA,
				 const double *sigmaE,
         const unsigned maxInterpolation,
				 unsigned *globalOffsetMultiplicator
				 );
#endif /* calc_sample_phi_ase_H */
