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


#include <thrust/scan.h>

#include <thrust_device_vector_nowarn.hpp> 
#include <map_rays_to_prisms.hpp>


using thrust::device_vector;
using thrust::host_vector;
using thrust::raw_pointer_cast;

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
__global__ void mapPrefixSumToPrisms(
    const unsigned numberOfPrisms,
    const unsigned raysPerSample,
    const unsigned reflectionSlices,
    const unsigned* raysPerPrism,
    const unsigned* prefixSum,
    unsigned *indicesOfPrisms,
    unsigned *numberOfReflections
    ){

  int id = threadIdx.x + (blockIdx.x * blockDim.x);
  // break if we have too many threads (this is likely)
  if(id >= numberOfPrisms*reflectionSlices) return;

  const unsigned count            = raysPerPrism[id];
  const unsigned startingPosition = prefixSum[id];
  const unsigned reflection_i     = id / numberOfPrisms;
  const unsigned prism_i          = id % numberOfPrisms;

  for(unsigned i=0; i < count ; ++i){
    indicesOfPrisms[startingPosition + i] = prism_i;     
    numberOfReflections[startingPosition + i] = reflection_i; 
  }
}

void mapRaysToPrisms(
    device_vector<unsigned> &indicesOfPrisms, 
    device_vector<unsigned> &numberOfReflections,
    const device_vector<unsigned> &raysPerPrism, 
    device_vector<unsigned> &prefixSum, 
    const unsigned reflectionSlices,
    const unsigned raysPerSample,
    const unsigned numberOfPrisms
    ){

  // blocksize chosen by occupancyCalculator
  const unsigned blocksize = 256;
  const unsigned gridsize  = (raysPerPrism.size()+blocksize-1)/blocksize;

  thrust::exclusive_scan(raysPerPrism.begin(), raysPerPrism.end(),prefixSum.begin());

  mapPrefixSumToPrisms <<<gridsize,blocksize>>> (
      numberOfPrisms, 
      raysPerSample, 
      reflectionSlices,
      raw_pointer_cast( &raysPerPrism[0] ),
      raw_pointer_cast( &prefixSum[0] ), 
      raw_pointer_cast( &indicesOfPrisms[0] ),
      raw_pointer_cast( &numberOfReflections[0] )
      );
}
