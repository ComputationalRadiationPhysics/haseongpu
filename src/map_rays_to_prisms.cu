#include "map_rays_to_prisms.h"
#include "cudachecks.h"
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <iterator>
#include <thrust/scan.h>
#include <thrust/device_vector.h>


using thrust::device_vector;
using thrust::host_vector;
using thrust::raw_pointer_cast;

__global__ void mapPrefixSumToPrisms(
    const unsigned numberOfPrisms,
    const unsigned raysPerSample,
    const unsigned reflectionSlices,
    const unsigned* raysPerPrism,
    const unsigned* prefixSum,
    unsigned *indicesOfPrisms,
    unsigned *numberOfReflections
    ){

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  if(id >= numberOfPrisms*reflectionSlices) return;

  const unsigned count = raysPerPrism[id];
  const unsigned startingPosition = prefixSum[id];
  const unsigned reflection_i = id/numberOfPrisms;
  const unsigned prism_i      = id%numberOfPrisms;

  for(unsigned i=0; i<count ; ++i){
    indicesOfPrisms[startingPosition + i] = prism_i;     
    numberOfReflections[startingPosition + i] = reflection_i; 
  }
}



void GPU_algorithm(
    const unsigned numberOfPrisms, 
    const unsigned raysPerSample,
    const unsigned reflectionSlices,
    const thrust::device_vector<unsigned>& raysPerPrism, 
    thrust::device_vector<unsigned> &prefixSum, 
    thrust::device_vector<unsigned> &indicesOfPrisms, 
    thrust::device_vector<unsigned> &numberOfReflections
    )
{
  const unsigned blocksize = 256;
  const unsigned gridsize  = (raysPerPrism.size()+blocksize-1)/blocksize;

  thrust::exclusive_scan(raysPerPrism.begin(), raysPerPrism.end(),prefixSum.begin());

  CUDA_CHECK_KERNEL_SYNC(mapPrefixSumToPrisms <<<gridsize,blocksize>>> (
      numberOfPrisms, 
      raysPerSample, 
      reflectionSlices,
      raw_pointer_cast( &raysPerPrism[0] ),
      raw_pointer_cast( &prefixSum[0] ), 
      raw_pointer_cast( &indicesOfPrisms[0] ),
      raw_pointer_cast( &numberOfReflections[0] )
      ));

}


void CPU_algorithm(
    thrust::host_vector<unsigned>& indicesOfPrisms, 
    const thrust::host_vector<unsigned>& raysPerPrism, 
    const unsigned numberOfPrisms, 
    const unsigned raysPerSample,
    thrust::host_vector<unsigned>& numberOfReflections,
    const unsigned reflectionSlices
    ){
  unsigned absoluteRay = 0;
  for(unsigned reflection_i=0; reflection_i < reflectionSlices ; ++reflection_i){
    unsigned reflectionOffset = reflection_i * numberOfPrisms;

    for(unsigned prism_i=0 ; prism_i < numberOfPrisms; ++prism_i){
      for(unsigned ray_i=0; ray_i < raysPerPrism[prism_i + reflectionOffset]; ++ray_i){
        indicesOfPrisms[absoluteRay] = prism_i;
        numberOfReflections[absoluteRay] = reflection_i;
        absoluteRay++;
        assert(absoluteRay <= raysPerSample);
      }
    }
  }
}

void mapRaysToPrisms(
    thrust::device_vector<unsigned> &indicesOfPrisms,
    thrust::device_vector<unsigned> &numberOfReflections,
    const thrust::device_vector<unsigned> &raysPerPrism,
    thrust::device_vector<unsigned> &prefixSum,
    const unsigned reflectionSlices,
    const unsigned raysPerSample,
    const unsigned numberOfPrisms
    ){


  //fill(indicesOfPrisms.begin(),indicesOfPrisms.end(),0);
  //fill(numberOfReflections.begin(),numberOfReflections.end(),0);

  //time_t before_GPU = clock();
  GPU_algorithm(
      numberOfPrisms,
      raysPerSample,
      reflectionSlices,
      raysPerPrism,
      prefixSum,
      indicesOfPrisms,
      numberOfReflections
      );
  //time_t after_GPU = clock();

  // only for error-checking!
  //time_t before_CPU = clock();
  //host_vector<unsigned> indicesOfPrisms2(indicesOfPrisms);
  //host_vector<unsigned> numberOfReflections2(numberOfReflections);
  //CPU_algorithm(indicesOfPrisms2, host_vector<unsigned>(raysPerPrism), numberOfPrisms, raysPerSample,numberOfReflections2,reflectionSlices);
  //indicesOfPrisms = indicesOfPrisms2;
  //numberOfReflections = numberOfReflections2; 
  //time_t after_CPU = clock();

  //some timing
  //int timeGPU = after_GPU - before_GPU;
  //int timeCPU = after_CPU - before_CPU;
  //std::cout << "time GPU including malloc: " << timeGPU/1000 << "k Cycles" << std::endl;
  //std::cout << "time CPU: " << timeCPU/1000 << "k Cycles" << std::endl;

  // some errorchecking
  //for(unsigned i=0; i<indicesOfPrisms2.size(); ++i){
  //  assert(indicesOfPrisms2[i] == indicesOfPrisms[i]);
  //}
}
