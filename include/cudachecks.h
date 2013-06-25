#ifndef CUDACHECKS_H
#define CUDACHECKS_H 

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

#endif /* CUDACHECKS_H */
