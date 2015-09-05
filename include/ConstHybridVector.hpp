#pragma once

// STL
#include <vector> /* std::vector */

// ALPAKA
#include <alpaka/alpaka.hpp> /* ALPAKA_FN_HOST_ACC */


//#include <cuda_utils.hpp> /* copyToDevice */
//#include <cuda_runtime_api.h>
//#include <host_defines.h> /* __host__ __device__ */


/**
 * @brief Vector on host and array on device
 *        with transparent access. Access is
 *        triggered based on a compiler macro.
 *
 **/
template<typename  T_Data, typename T_Dev>
class ConstHybridVector {
public:

    ConstHybridVector(std::vector<T_Data> &srcV, T_Dev &dev) :
	hostV(srcV),
	deviceV(alpaka::mem::buf::alloc<T_Data>( dev, srcV.size())) {
      
    }
    
  ALPAKA_FN_HOST_ACC T_Data at(int i) const{
      /*
#ifdef __CUDA_ARCH__
    return deviceV[i];
#else
    return hostV.at(i);
#endif
      */
  }

    ALPAKA_FN_HOST_ACC const T_Data operator[] (int i) const {
      /*
#ifdef __CUDA_ARCH__
    return deviceV[i];
#else
    return hostV[i];
#endif
      */
  }

  ALPAKA_FN_HOST  operator T_Data*(){
    return &(hostV.at(0));
  }

  ALPAKA_FN_HOST T_Data* toArray() const{
    std::vector<T_Data> copyV = hostV;
    return &(copyV.at(0));
  }

  ALPAKA_FN_HOST std::vector<T_Data> toVector() const{
    return hostV;
  }

    ALPAKA_FN_HOST void free() {
	/*
        cudaFree(deviceV);
	*/
    }
  
private:
    alpaka::mem::buf::Buf<T_Dev, T_Data, std::size_t, std::size_t> deviceV;
    std::vector<T_Data> hostV;

};
