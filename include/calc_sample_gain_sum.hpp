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

// HASEonGPU
#include <geometry.hpp>      /* generateRay */
#include <propagate_ray.hpp> /* propagateRay */

// FIXIT: We should have some random number generator object and not initilize it again and again !
template <typename T_Acc>
ALPAKA_FN_ACC unsigned genRndSigmas(T_Acc const &acc, unsigned length) {
    using Gen =   decltype(alpaka::rand::generator::createDefault(std::declval<T_Acc const &>(),
								  std::declval<uint32_t &>(),
								  std::declval<uint32_t &>()));
    using Dist =  decltype(alpaka::rand::distribution::createUniformReal<float>(std::declval<T_Acc const &>()));
    
    // FIXIT: the orginal used blockid.x to 
    Gen gen(alpaka::rand::generator::createDefault(acc, alpaka::idx::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0], 0));
    Dist dist(alpaka::rand::distribution::createUniformReal<float>(acc));
    
    return static_cast<unsigned>(dist(gen) * (length-1));
}


/**
 * @brief get the offset for accessing indicesOfPrisms and numberOfReflectionSlices.
 *        Warning: works only for a blocksize of 128 threads!
 *
 * @param blockOffset shared memory location that holds the offset for the whole block
 * @param raysPerSample number of raysPerSample (can be any number higher than raysPerSample/blocksize)
 * @param globalOffsetMultiplicator is incremented by 1 each time a block asks for a new workload
 * @return an unused offset in the global arrays indicesOfPrisms/numberOfReflectionSlices
 *
 */
template <typename T_Acc>
ALPAKA_FN_ACC unsigned getRayNumberBlockbased(T_Acc const & acc, unsigned* const blockOffset, unsigned const &raysPerSample, unsigned * globalOffsetMultiplicator){

    // FIXIT: Extra sync since alpaka
    alpaka::block::sync::syncBlockThreads(acc);

    // The first thread in the threadblock increases the globalOffsetMultiplicator (without real limit)
    if(alpaka::idx::getIdx<alpaka::Block, alpaka::Threads>(acc)[0] == 0){
	// blockOffset is the new value of the globalOffsetMultiplicator
	blockOffset[0] = alpaka::atomic::atomicOp<alpaka::atomic::op::Add>(acc, globalOffsetMultiplicator, 1u);
	//std::cout << globalOffsetMultiplicator[0] << std::endl;
    }

    alpaka::block::sync::syncBlockThreads(acc);

    return alpaka::idx::getIdx<alpaka::Block, alpaka::Threads>(acc)[0] + (blockOffset[0] * alpaka::workdiv::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0]) ;
}


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
struct CalcSampleGainSumWithReflection {
    template <typename T_Acc,
	      typename T_Mesh>
    ALPAKA_FN_ACC void operator()(const T_Acc &acc,
				  const T_Mesh &mesh, 
				  const unsigned* indicesOfPrisms, 
				  const unsigned* numberOfReflectionSlices,
				  const double* importance,
				  const unsigned raysPerSample,
				  float *gainSum, 
				  float *gainSumSquare,
				  const unsigned sample_i,
				  const double *sigmaA, 
				  const double *sigmaE,
				  const unsigned maxInterpolation,
				  unsigned * const globalOffsetMultiplicator) const {

	unsigned rayNumber = 0;
	double gainSumTemp = 0;
	double gainSumSquareTemp = 0;
	Point samplePoint = mesh.getSamplePoint(sample_i);

	auto * blockOffset(alpaka::block::shared::allocArr<unsigned, 4>(acc)); // 4 in case of warp-based raynumber
	blockOffset[0] = 0;

	const unsigned nElementsPerThread = 1;

	auto localTId = alpaka::workdiv::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0];
	
	// One thread can compute multiple rays
	while(blockOffset[0] * localTId * nElementsPerThread < raysPerSample){

	    // the whole block gets a new offset (==workload)
	    //rayNumber = getRayNumberBlockbased(acc, blockOffset, raysPerSample, globalOffsetMultiplicator);
	    // HACK: ONLY WITH BLOCKSIZE 1, 1, 1
	    blockOffset[0] = alpaka::atomic::atomicOp<alpaka::atomic::op::Add>(acc, globalOffsetMultiplicator, 1u);

	    // rayNumber = alpaka::idx::getIdx<alpaka::Block, alpaka::Threads>(acc)[0] + (blockOffset[0] * alpaka::workdiv::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0]);

	    auto threadNumberOffset =
		(blockOffset[0] * alpaka::workdiv::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0] * nElementsPerThread) +
		(alpaka::idx::getIdx<alpaka::Block, alpaka::Threads>(acc)[0] * nElementsPerThread);
	    
	    for(unsigned rayNumber = threadNumberOffset; rayNumber < (threadNumberOffset + nElementsPerThread) && rayNumber < raysPerSample; ++rayNumber){

		// rayNumber =
		//     (blockOffset[0] * alpaka::workdiv::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0] * nElementsPerThread) +
		//     (alpaka::idx::getIdx<alpaka::Block, alpaka::Threads>(acc)[0] * nElementsPerThread) +
		//     nthRay;
		
		//rayNumber = alpaka::idx::getIdx<alpaka::Block, alpaka::Threads>(acc)[0] + (blockOffset[0] * alpaka::workdiv::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0])
		
		//if(rayNumber < raysPerSample) {

		// Get triangle/prism to start ray from
		unsigned startPrism             = indicesOfPrisms[rayNumber];
		unsigned reflection_i           = numberOfReflectionSlices[rayNumber];
		unsigned reflections            = (reflection_i + 1) / 2;
		ReflectionPlane reflectionPlane = (reflection_i % 2 == 0) ? BOTTOM_REFLECTION : TOP_REFLECTION;
		unsigned startLevel             = startPrism / mesh.numberOfTriangles;
		unsigned startTriangle          = startPrism - (mesh.numberOfTriangles * startLevel);
		unsigned reflectionOffset       = reflection_i * mesh.numberOfPrisms;
		Point startPoint                = mesh.genRndPoint(acc, startTriangle, startLevel);
	
		//get a random index in the wavelength array
		unsigned sigma_i                = genRndSigmas(acc, maxInterpolation);

		// Calculate reflections as different ray propagations
		double gain    = propagateRayWithReflection(startPoint, samplePoint, reflections, reflectionPlane, startLevel, startTriangle, mesh, sigmaA[sigma_i], sigmaE[sigma_i]);
		//std::cout << gain << std::endl;

		// include the stimulus from the starting prism and the importance of that ray
		gain          *= mesh.getBetaVolume(startPrism) * importance[startPrism + reflectionOffset];
    
		assert(!isnan(mesh.getBetaVolume(startPrism)));
		assert(!isnan(importance[startPrism + reflectionOffset]));
		assert(!isnan(gain));

		gainSumTemp       += gain;
		gainSumSquareTemp += gain * gain;
 
		//}
	    }


	}

	alpaka::atomic::atomicOp<alpaka::atomic::op::Add>(acc, gainSum,       static_cast<float>(gainSumTemp));
	alpaka::atomic::atomicOp<alpaka::atomic::op::Add>(acc, gainSumSquare, static_cast<float>(gainSumSquareTemp));

    }

};


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
    ALPAKA_FN_ACC void operator()(const T_Acc &acc,
				  const T_Mesh &mesh,
				  const unsigned* indicesOfPrisms, 
				  const double* importance,
				  const unsigned raysPerSample,
				  float * const gainSum,
				  float * const gainSumSquare,
				  const unsigned sample_i,
				  const double *sigmaA,
				  const double *sigmaE,
				  const unsigned maxInterpolation,
				  unsigned * globalOffsetMultiplicator) const {

	unsigned rayNumber       = 0; 
	double gainSumTemp       = 0;
	double gainSumSquareTemp = 0;
	Point samplePoint        = mesh.getSamplePoint(sample_i);
	
	auto * blockOffset(alpaka::block::shared::allocArr<unsigned, 4>(acc)); // 4 in case of warp-based raynumber
	blockOffset[0] = 0;
	// sync threads here
	
	// One thread can compute multiple rays
	// Need to replace this while true
	size_t textent = alpaka::workdiv::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0];
	while((blockOffset[0] * textent) < raysPerSample){
	    // the whole block gets a new offset (==workload)
	    rayNumber = getRayNumberBlockbased(acc, blockOffset, raysPerSample, globalOffsetMultiplicator);
	    if(rayNumber < raysPerSample) {

		// Get triangle/prism to start ray from
		unsigned startPrism             = indicesOfPrisms[rayNumber]; 
		unsigned startLevel             = startPrism/mesh.numberOfTriangles;
		unsigned startTriangle          = startPrism - (mesh.numberOfTriangles * startLevel);
		Point startPoint                = mesh.genRndPoint(acc, startTriangle, startLevel);
		Ray ray                         = generateRay(startPoint, samplePoint);

		// get a random index in the wavelength array
		unsigned sigma_i                = genRndSigmas(acc, maxInterpolation);
		assert(sigma_i < maxInterpolation);

		// calculate the gain for the whole ray at once
		double gain    = propagateRay(ray, &startLevel, &startTriangle, mesh, sigmaA[sigma_i], sigmaE[sigma_i]);
		//std::cout << gain << std::endl;

		//std::cout << rayNumber << " " << startPrism << std::endl;
		
		gain          /= ray.length * ray.length; // important, since usually done in the reflection device function

		// include the stimulus from the starting prism and the importance of that ray
		gain          *= mesh.getBetaVolume(startPrism) * importance[startPrism];

		//std::cout << gain << " " << importance[startPrism] << std::endl;
		
		gainSumTemp       += gain;
		gainSumSquareTemp += gain * gain;
		//std::cout << rayNumber << " " << globalOffsetMultiplicator[0] << " "<< omp_get_thread_num()<<std::endl;	
	    }
	    

	}

	

	
	alpaka::atomic::atomicOp<alpaka::atomic::op::Add>(acc, gainSum,       static_cast<float>(gainSumTemp));
	alpaka::atomic::atomicOp<alpaka::atomic::op::Add>(acc, gainSumSquare, static_cast<float>(gainSumSquareTemp));
	
    }

};
