/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASENonGPU
 *
 * HASENonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASENonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASENonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


#pragma once
#include <stdio.h> /* fprintf, printf, stderr */
#include <cuda_runtime_api.h> /* cuda constants etc. */
#include <stdlib.h> /* exit */

#define CUDA_CHECK(cmd) {cudaError_t error = cmd; if(error!=cudaSuccess){printf("<%s>:%i ",__FILE__,__LINE__); printf("[CUDA] Error: %s\n", cudaGetErrorString(error));}}
/*start kernel, wait for finish and check errors*/
#define CUDA_CHECK_KERNEL_SYNC(...) __VA_ARGS__;CUDA_CHECK(cudaDeviceSynchronize())
/*only check if kernel start is valid*/
#define CUDA_CHECK_KERNEL(...) __VA_ARGS__;CUDA_CHECK(cudaGetLastError())

#define CUDA_CHECK_RETURN(value) {					\
    cudaError_t _m_cudaStat = value;					\
    if (_m_cudaStat != cudaSuccess) {					\
      fprintf(stderr, "Error %s at line %d in file %s\n",		\
	      cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);	\
      exit(1);								\
    }									\
  }
#define CUDA_CALL(x) do { if((x) != cudaSuccess) {	\
      printf("Error at %s:%d\n",__FILE__,__LINE__);	\
      return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x) != CURAND_STATUS_SUCCESS) {	\
      printf("Error at %s:%d\n",__FILE__,__LINE__);		\
      return EXIT_FAILURE;}} while(0)

