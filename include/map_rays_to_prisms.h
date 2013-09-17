#include <thrust/device_vector.h>

void mapRaysToPrisms(
    thrust::device_vector<unsigned> &indicesOfPrisms,
    thrust::device_vector<unsigned> &numberOfReflections,
    const thrust::device_vector<unsigned> &raysPerPrism,
    thrust::device_vector<unsigned> &prefixSum,
    const unsigned reflectionSlices,
    const unsigned raysPerSample,
    const unsigned numberOfPrisms
    );
