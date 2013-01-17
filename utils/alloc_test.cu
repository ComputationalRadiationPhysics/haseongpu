#include <stdio.h>
#include <assert.h>

#define CUDA_CHECK(cmd) {cudaError_t error = cmd; if(error!=cudaSuccess){printf("<%s>:%i ",__FILE__,__LINE__); printf("[CUDA] Error: %s\n", cudaGetErrorString(error));}}
/*start kernel, wait for finish and check errors*/
#define CUDA_CHECK_KERNEL_SYNC(...) __VA_ARGS__;CUDA_CHECK(cudaDeviceSynchronize())
/*only check if kernel start is valid*/
#define CUDA_CHECK_KERNEL(...) __VA_ARGS__;CUDA_CHECK(cudaGetLastError())

struct octNode {
  int test;
  int *prism_list;
};

__global__ void test(struct octNode *octtree) {
  int idx = threadIdx.x;
  octtree[idx].test = octtree[idx].prism_list[idx%2];
}

int main (int argc, char **argv){
  struct octNode octtree[100];
  struct octNode *d_octtree;

  for(int i=0; i<100; ++i) {
    int prisms[2] = {1,2};
    octtree[i].test = 0;
    CUDA_CHECK(cudaMalloc((void**) &(octtree[i].prism_list), 2*sizeof(int)));
    CUDA_CHECK(cudaMemcpy(octtree[i].prism_list, prisms, 2*sizeof(int), cudaMemcpyHostToDevice));
  }

  CUDA_CHECK(cudaMalloc((void**)&d_octtree, 100*sizeof(struct octNode)));
  CUDA_CHECK(cudaMemcpy(d_octtree, octtree, 100*sizeof(struct octNode), cudaMemcpyHostToDevice));

  CUDA_CHECK_KERNEL_SYNC(test<<<1,100>>>(d_octtree));

  CUDA_CHECK(cudaMemcpy(octtree, d_octtree, 100*sizeof(struct octNode), cudaMemcpyDeviceToHost));

  for(int i=0; i<100; ++i) {
    printf("%i\n", octtree[i].test);
  }

  return 0;
}
