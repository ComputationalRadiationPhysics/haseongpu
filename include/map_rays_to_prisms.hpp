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
 * @author Carlchristian Eckert
 * @license GPLv3
 */
#pragma once

#include <alpaka/alpaka.hpp>

// /**
//  * @brief takes an array of values and writes it into a unary representation. The value
//  *        from the input at index i describes the number of elements to write,
//  *        whereas i itself is the new value to be stored in the output array.
//  *
//  *        example:
//  *        raysPerPrism [3,0,2,1] 
//  *
//  *        3 elements in output should be 0 
//  *        0 elements in output should be 1 
//  *        2 elements in output should be 2 
//  *        1 element  in output should be 3
//  *        
//  *        resulting output arrays:
//  *        [0 0 0 2 2 3] (indicesOfPrisms)
//  *
//  *        output numberOfReflections is handled in a similar way
//  *
//  *
//  * @param indicesOfPrisms a pointer to the OUTPUT generated like described in the example
//  * @param numberOfReflections a pointer to the OUTPUT similar to indicesOfPrisms
//  * @param raysPerPrism the input array
//  *                       must be equal to the the length of the prefixSum.
//  * @param reflectionSlices the number of reflectionSlices. see numberOfPrisms
//  * @param raysPerSample the size of indicesOfPrisms/numberOfReflections. Actually 
//  *                      identical to the sum of all values in raysPerPrism
//  * @param numberOfPrisms the number of prisms. numberOfPrisms * reflectionSlices
//  */
// void mapRaysToPrisms(unsigned *indicesOfPrisms, 
// 		     unsigned const *numberOfReflections,
// 		     unsigned const *raysPerPrism, 
// 		     unsigned const *prefixSum, 
// 		     const unsigned reflectionSlices,
// 		     const unsigned raysPerSample,
// 		     const unsigned numberOfPrisms
// 		     ){

//   // blocksize chosen by occupancyCalculator
//   const unsigned blocksize = 256;
//   const unsigned gridsize  = (raysPerPrism.size()+blocksize-1)/blocksize;

//   thrust::exclusive_scan(raysPerPrism.begin(), raysPerPrism.end(), prefixSum.begin());

//   exclusivePrefixSum();

  

//   mapPrefixSumToPrisms <<<gridsize,blocksize>>> (
// 						 numberOfPrisms, 
// 						 raysPerSample, 
// 						 reflectionSlices,
// 						 raw_pointer_cast( &raysPerPrism[0] ),
// 						 raw_pointer_cast( &prefixSum[0] ), 
// 						 raw_pointer_cast( &indicesOfPrisms[0] ),
// 						 raw_pointer_cast( &numberOfReflections[0] )
// 						 );
// }


/**
 * @brief takes a prefix sum (obtained by an exclusive scan over raysPerPrism) 
 *        and writes it into a unary representation. The value from the 
 *        prefixSum at index i describes the offset where to start writing,
 *        whereas i itself is the new value to be stored in the output array.
 *
 *        example:
 *        raysPerPrism [3,0,2,1] 
 *
 *        -> exclusive prefixSum [0,3,3,5] 
 *
 *        beginning from place 0 in the output should be 0 (length 3 according to raysPerPrism[0])
 *        beginning from place 3 in the output should be 1 (EMPTY range at raysPerPrism[1])
 *        beginning from place 3 in the output should be 2 (length 2 according to raysPerPrism[2])
 *        beginning from place 5 in the output should be 3 (length 1 according to raysPerPrism[3])
 *        
 *        resulting output arrays:
 *        [0 0 0 2 2 3] (indicesOfPrisms)
 *
 *        output numberOfReflections is handled in a similar way
 *
 *
 * @param numberOfPrisms the number of prisms. numberOfPrisms * reflectionSlices
 *                       must be equal to the the length of the prefixSum.
 * @param raysPerSample the size of indicesOfPrisms/numberOfReflections. Actually 
 *                      identical to the sum of all values in raysPerPrism
 * @param reflectionSlices the number of reflectionSlices. see numberOfPrisms
 * @param raysPerPrism the input array from which prefixSum was generated
 * @param prefixSum the prefixSum generated from raysPerPrism
 * @param indicesOfPrisms a pointer to the OUTPUT generated like described in the example
 * @param numberOfReflections a pointer to the OUTPUT similar to indicesOfPrisms
 */

struct MapPrefixSumToPrisms {
    template <typename T_Acc>
    ALPAKA_FN_ACC void operator()(T_Acc const &acc,
				  const unsigned numberOfPrisms,
				  const unsigned reflectionSlices,
				  const unsigned* raysPerPrism,
				  const unsigned* prefixSum,
				  unsigned *indicesOfPrisms,
				  unsigned *numberOfReflections) const {

	auto threadsVec = alpaka::workdiv::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc) *  alpaka::workdiv::getWorkDiv<alpaka::Block, alpaka::Threads>(acc);
	auto nThreads = threadsVec[0] * threadsVec[1];
	
	for(unsigned id = alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0]; id < numberOfPrisms * reflectionSlices; id += nThreads){
	
	    const unsigned count            = raysPerPrism[id];
	    const unsigned startingPosition = prefixSum[id];
	    const unsigned reflection_i     = id / numberOfPrisms;
	    const unsigned prism_i          = id % numberOfPrisms;

	    for(unsigned i=0; i < count ; ++i){
		indicesOfPrisms[startingPosition + i] = prism_i;
		numberOfReflections[startingPosition + i] = reflection_i; 
	    }
	
	}
	
    }

};
template <typename T_Acc, typename T_Workdiv, typename T_Stream>
void mapPrefixSumToPrisms( T_Stream &stream,
			   T_Workdiv &workdiv,
			   const unsigned numberOfPrisms,
			   const unsigned reflectionSlices,
			   const unsigned* raysPerPrism,
			   const unsigned* prefixSum,
			   unsigned *indicesOfPrisms,
			   unsigned *numberOfReflections){

    MapPrefixSumToPrisms mapPrefixSumToPrisms;

    auto const exec (alpaka::exec::create<T_Acc>(workdiv,
						 mapPrefixSumToPrisms,
						 numberOfPrisms,
						 reflectionSlices,
						 raysPerPrism,
						 prefixSum,
						 indicesOfPrisms,
						 numberOfReflections));
    alpaka::stream::enqueue(stream, exec);
    

}
