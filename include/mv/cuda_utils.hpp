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


#pragma once
#include <vector>

#include <cudachecks.hpp>
#include <cuda_runtime_api.h>

static const unsigned MIN_COMPUTE_CAPABILITY_MAJOR = 2;
static const unsigned MIN_COMPUTE_CAPABILITY_MINOR = 0;


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




/** 
 * @brief Queries for devices on the running mashine and collects
 *        them on the devices array. Set the first device in this 
 *        array as computation-device. On Errors the programm will
 *        be stoped by exit(). 
 * 
 * @param maxGpus max. devices which should be allocated
 * @return vector of possible devices
 */
std::vector<unsigned> getFreeDevices(unsigned maxGpus);

