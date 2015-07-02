#pragma once

// STL
#include <vector> /* std::vector */

#include <cuda_utils.hpp> /* copyToDevice */


/**
 * @brief Vector on host and array on device
 *        with transparent access. Access is
 *        triggered based on a compiler macro.
 *
 **/
template<class T>
class ConstHybridVector {
public:

  ConstHybridVector(std::vector<T> &srcV) :
    hostV(srcV),
    deviceV(copyToDevice(srcV)){
      
  }

  __forceinline__ __host__ __device__ T at(int i) const{
#ifdef __CUDA_ARCH__
    return deviceV[i];
#else
    return hostV.at(i);
#endif
  }

  __forceinline__ __host__ __device__ const T operator[] (int i) const {
#ifdef __CUDA_ARCH__
    return deviceV[i];
#else
    return hostV[i];
#endif
  }

  __host__  operator T*(){
    return &(hostV.at(0));
  }

  __host__ T* toArray() const{
    std::vector<T> copyV = hostV;
    return &(copyV.at(0));
  }

  __host__ std::vector<T> toVector() const{
    return hostV;
  }

    __host__ void free() {
        cudaFree(deviceV);
    }
  
private:
  T *deviceV;
  std::vector<T> hostV;

};
