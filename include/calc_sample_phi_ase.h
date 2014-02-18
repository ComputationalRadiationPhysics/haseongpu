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

#ifndef calc_sample_phi_ase_H
#define calc_sample_phi_ase_H

#include <curand_mtgp32dc_p_11213.h>
#include <mesh.h>

/**
 * @brief: calculate the sum of all gains for a single sample point s_i, accounting
 *		   for reflections on the top and bottom surface. Different cladding, 
 *		   coating and refractive indices are stored in the mesh.
 *
 * @param globalState State for random number generation (mersenne twister).
 *                    The state need to be initialized before. See
 *                    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MTGP/
 *                    for more information.
 * @param mesh the gain medium mesh 
 * @param indicesOfPrisms a mapping from rays to prisms. The value x on the
 *	      n-th position denotes that ray n is supposed to start propagating 
 *	      from the x-th prism
 * @param numberOfReflectionSlices similar to indicesOfPrisms, but denotes the
 *        reflection slice from which each ray starts (see importance sampling)
 * @param importance the importance for each prism, obtained through the
 *        importance sampling kernel
 * @param raysPerSample the number of rays to use for the sample point s_i
 * @param gainSum the sum of all gain values for sample point s_i
 * @param gainSumSquare the sum of all (gain^2) values for sample point s_i
 *        useful for finding the mean square error
 * @param sample_i the index of the sample point to sample (s_i)
 * @param sigmaA the values for absorption cross section in the selected wavelength interval
 * @param sigmaE the values for emission cross section in the selected wavelength interval
 * @param maxInterpolation the number of elements contained in sigmaA or sigmaE
 * @param globalOffsetMultiplicator holds the number of work-requests done by the blocks
 *        (see getRayNumberBlockbased())
 *
 */
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

/**
 * @brief: calculate the sum of all gains for a single sample point s_i, ignoring
 *         reflections. Different cladding, is stored in the mesh.
 *
 * @param globalState State for random number generation (mersenne twister).
 *                    The state need to be initialized before. See
 *                    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MTGP/
 *                    for more information.
 * @param mesh the gain medium mesh 
 * @param indicesOfPrisms a mapping from rays to prisms. The value x on the
 *	      n-th position denotes that ray n is supposed to start propagating 
 *	      from the x-th prism
 * @param importance the importance for each prism, obtained through the
 *        importance sampling kernel
 * @param raysPerSample the number of rays to use for the sample point s_i
 * @param gainSum the sum of all gain values for sample point s_i
 * @param gainSumSquare the sum of all (gain^2) values for sample point s_i
 *        useful for finding the mean square error
 * @param sample_i the index of the sample point to sample (s_i)
 * @param sigmaA the values for absorption cross section in the selected wavelength interval
 * @param sigmaE the values for emission cross section in the selected wavelength interval
 * @param maxInterpolation the number of elements contained in sigmaA or sigmaE
 * @param globalOffsetMultiplicator holds the number of work-requests done by the blocks
 *        (see getRayNumberBlockbased())
 *
 */
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
