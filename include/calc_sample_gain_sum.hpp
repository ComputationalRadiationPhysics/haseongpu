/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


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

#pragma once

// Alpaka
#include <alpaka/alpaka.hpp>

// /**
//  * @brief: calculate the sum of all gains for a single sample point s_i, accounting
//  *		   for reflections on the top and bottom surface. Different cladding, 
//  *		   coating and refractive indices are stored in the mesh.
//  *
//  * @param globalState State for random number generation (mersenne twister).
//  *                    The state need to be initialized before. See
//  *                    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MTGP/
//  *                    for more information.
//  * @param mesh the gain medium mesh 
//  * @param indicesOfPrisms a mapping from rays to prisms. The value x on the
//  *	      n-th position denotes that ray n is supposed to start propagating 
//  *	      from the x-th prism
//  * @param numberOfReflectionSlices similar to indicesOfPrisms, but denotes the
//  *        reflection slice from which each ray starts (see importance sampling)
//  * @param importance the importance for each prism, obtained through the
//  *        importance sampling kernel
//  * @param raysPerSample the number of rays to use for the sample point s_i
//  * @param gainSum the sum of all gain values for sample point s_i
//  * @param gainSumSquare the sum of all (gain^2) values for sample point s_i
//  *        useful for finding the mean square error
//  * @param sample_i the index of the sample point to sample (s_i)
//  * @param sigmaA the values for absorption cross section in the selected wavelength interval
//  * @param sigmaE the values for emission cross section in the selected wavelength interval
//  * @param maxInterpolation the number of elements contained in sigmaA or sigmaE
//  * @param globalOffsetMultiplicator holds the number of work-requests done by the blocks
//  *        (see getRayNumberBlockbased())
//  *
//  */
// __global__ void calcSampleGainSumWithReflection(curandStateMtgp32* globalState, 
// 						const Mesh mesh,
// 						const unsigned* indicesOfPrisms, 
// 						const unsigned* numberOfReflections,
// 						const double* importance,
// 						const unsigned raysPerSample,
// 						float *gainSum,
// 						float *gainSumSquare,
// 						const unsigned sample_i,
// 						const double *sigmaA,
// 						const double *sigmaE,
// 						const unsigned maxInterpolation,
// 						unsigned *globalOffsetMultiplicator
// 						);

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
struct CalcSampleGainSum {
    template <typename T_Acc,
	      typename T_Mesh>
    ALPAKA_FN_ACC void operator()(const T_Mesh mesh,
				  const unsigned* indicesOfPrisms, 
				  const double* importance,
				  const unsigned raysPerSample,
				  float *gainSum,
				  float *gainSumSquare,
				  const unsigned sample_i,
				  const double *sigmaA,
				  const double *sigmaE,
				  const unsigned maxInterpolation,
				  unsigned *globalOffsetMultiplicator) const {

	int rayNumber = 0; 
	double gainSumTemp = 0;
	double gainSumSquareTemp = 0;
	Point samplePoint = mesh.getSamplePoint(sample_i);
	// __shared__ unsigned blockOffset[4]; // 4 in case of warp-based raynumber
  
	// // One thread can compute multiple rays
	// while(true){
	//     // the whole block gets a new offset (==workload)
	//     rayNumber = getRayNumberBlockbased(blockOffset,raysPerSample,globalOffsetMultiplicator);
	//     if(rayNumber>=raysPerSample) break;

	//     // Get triangle/prism to start ray from
	//     unsigned startPrism             = indicesOfPrisms[rayNumber];
	//     unsigned startLevel             = startPrism/mesh.numberOfTriangles;
	//     unsigned startTriangle          = startPrism - (mesh.numberOfTriangles * startLevel);
	//     Point startPoint                = mesh.genRndPoint(startTriangle, startLevel, globalState);
	//     Ray ray                         = generateRay(startPoint, samplePoint);

	//     // get a random index in the wavelength array
	//     unsigned sigma_i                = genRndSigmas(lambdaResolution, globalState);
	//     assert(sigma_i < lambdaResolution);

	//     // calculate the gain for the whole ray at once
	//     double gain    = propagateRay(ray, &startLevel, &startTriangle, mesh, sigmaA[sigma_i], sigmaE[sigma_i]);
	//     gain          /= ray.length * ray.length; // important, since usually done in the reflection device function

	//     // include the stimulus from the starting prism and the importance of that ray
	//     gain          *= mesh.getBetaVolume(startPrism) * importance[startPrism];

	//     gainSumTemp       += gain;
	//     gainSumSquareTemp += gain * gain;

	// }
	// atomicAdd(&(gainSum[0]), float(gainSumTemp));
	// atomicAdd(&(gainSumSquare[0]), float(gainSumSquareTemp));


	
    }

};
