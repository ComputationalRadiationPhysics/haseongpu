#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <vector>
#include <iostream>
#include <cudachecks.h>

/**
 * @brief Copy data from host to device
 *
 **/
template <class T>
T* copyToDevice(const std::vector<T> &v){
  T* deviceV;
  CUDA_CHECK_RETURN(cudaMalloc((void**)&deviceV,  v.size()* sizeof(T)));
  CUDA_CHECK_RETURN(cudaMemcpy(deviceV, &(v[0]), v.size() * sizeof(T), cudaMemcpyHostToDevice));

  return deviceV;
}

template <class T>
void copyToDevice(const std::vector<T> &v,  T* deviceV){
  CUDA_CHECK_RETURN(cudaMemcpy(deviceV, &(v[0]), v.size() * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
void copyToDevice(const std::vector<T> &v,  T* deviceV, unsigned size){
  assert(size <= v.size());
  CUDA_CHECK_RETURN(cudaMemcpy(deviceV, &(v[0]), size * sizeof(T), cudaMemcpyHostToDevice));
}

template <class T>
T* copyToDevice(T a){
  T* deviceV;
  CUDA_CHECK_RETURN(cudaMalloc((void**)&deviceV, sizeof(T)));
  CUDA_CHECK_RETURN(cudaMemcpy(deviceV, &a, sizeof(T), cudaMemcpyHostToDevice));

  return deviceV;
}

template <class T>
void copyToDevice(const T a, T* deviceV){
  CUDA_CHECK_RETURN(cudaMemcpy(deviceV, &a, sizeof(T), cudaMemcpyHostToDevice));
}


/**
 * @brief Copy data from device to host
 *
 **/
template <class T>
void copyFromDevice(std::vector<T> &v, const T* deviceV){
  CUDA_CHECK_RETURN(cudaMemcpy(&(v[0]), deviceV, v.size() * sizeof(T), cudaMemcpyDeviceToHost));

}

template <class T>
T copyFromDevice(const T* deviceV){
  T a;
  CUDA_CHECK_RETURN(cudaMemcpy(&a, deviceV, sizeof(T), cudaMemcpyDeviceToHost));
  
  return a;
}



#endif /* CUDA_UTILS_H */
