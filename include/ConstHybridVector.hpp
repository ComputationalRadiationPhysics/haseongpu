#pragma once

// STL
#include <cuda_runtime_api.h>
#include <cuda_utils.hpp> /* copyToDevice */

#include <span>
#include <vector> /* std::vector */

/**
 * @brief Vector on host and array on device
 *        with transparent access. Access is
 *        triggered based on a compiler macro.
 *
 **/
template<class T>
class ConstHybridVector
{
public:
    explicit ConstHybridVector(std::vector<T>& srcV) : hostV(srcV), deviceV(copyToDevice(srcV))
    {
    }

    __forceinline__ __host__ __device__ T at(int i) const
    {
#ifdef __CUDA_ARCH__
        return deviceV[i];
#else
        return hostV.at(i);
#endif
    }

    __forceinline__ __host__ __device__ T operator[](int i) const
    {
#ifdef __CUDA_ARCH__
        return deviceV[i];
#else
        return hostV[i];
#endif
    }

    __host__ operator T*()
    {
        return &(hostV.at(0));
    }

    __host__ auto toArray() const
    {
        return std::span<T const>(hostV.data(), hostV.size());
    }

    __host__ std::vector<T> toVector() const
    {
        return hostV;
    }

    __host__ void free()
    {
        CUDA_CHECK_RETURN(cudaFree(deviceV));
    }

private:
    T* deviceV;
    std::vector<T> hostV;
};
