/**
 * @author Carlchristian Eckert
 * @license GPLv3
 */
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
 * @param indicesOfPrisms a pointer to the OUTPUT generated like described in the example
 * @param numberOfReflections a pointer to the OUTPUT similar to indicesOfPrisms
 * @param raysPerPrism the input array
 *                       must be equal to the the length of the prefixSum.
 * @param reflectionSlices the number of reflectionSlices. see numberOfPrisms
 * @param raysPerSample the size of indicesOfPrisms/numberOfReflections. Actually 
 *                      identical to the sum of all values in raysPerPrism
 * @param numberOfPrisms the number of prisms. numberOfPrisms * reflectionSlices
 */
void mapRaysToPrisms(
    thrust::device_vector<unsigned> &indicesOfPrisms,
    thrust::device_vector<unsigned> &numberOfReflections,
    const thrust::device_vector<unsigned> &raysPerPrism,
    thrust::device_vector<unsigned> &prefixSum,
    const unsigned reflectionSlices,
    const unsigned raysPerSample,
    const unsigned numberOfPrisms
    );
