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
#include <thrust/device_vector.h>

/**
 * @brief takes an array of values and writes it into a unary representation. The value
 *        from the input at index i describes the number of elements to write,
 *        whereas i itself is the new value to be stored in the output array.
 *
 *        example:
 *        raysPerPrism [3,0,2,1] 
 *
 *        3 elements in output should be 0 
 *        0 elements in output should be 1 
 *        2 elements in output should be 2 
 *        1 element  in output should be 3
 *        
 *        resulting output arrays:
 *        [0 0 0 2 2 3] (indicesOfPrisms)
 *
 *        output numberOfReflections is handled in a similar way
 *
 *
 * @param indicesOfPrisms a reference to the OUTPUT generated like described in the example
 * @param raysPerPrismStart iterator to the input array
 * @param raysPerPrismEnd iterator to the input array
 */
void mapRaysToPrisms(
        thrust::device_vector<unsigned> &indicesOfPrisms,
        const thrust::device_vector<unsigned>::iterator raysPerPrismBegin,
        const thrust::device_vector<unsigned>::iterator raysPerPrismEnd,
        const thrust::device_vector<unsigned>::iterator prefixSumBegin,
        const thrust::device_vector<unsigned>::iterator prefixSumEnd,
        const unsigned offset
        );
