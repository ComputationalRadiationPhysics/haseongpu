#include <iostream>
#include <cstdio>
#include <cassert>

#include "grid.h"

using std::endl;
using std::cout;

#define CUDA_CHECK(cmd) {cudaError_t error = cmd; if(error!=cudaSuccess){printf("<%s>:%i ",__FILE__,__LINE__); printf("[CUDA] Error: %s\n", cudaGetErrorString(error));}}
/*start kernel, wait for finish and check errors*/
#define CUDA_CHECK_KERNEL_SYNC(...) __VA_ARGS__;CUDA_CHECK(cudaDeviceSynchronize())
/*only check if kernel start is valid*/
#define CUDA_CHECK_KERNEL(...) __VA_ARGS__;CUDA_CHECK(cudaGetLastError())

//important!!!
//declare variable for the gird on GPU
//-> declare d_grid not as pointer
//commit d_grid  -> &d_grid to the function
//use it as call-by-value-Variable in the Kernel
void copyGrid2GPU(const Grid *grid, Grid *d_grid)
{
	//number of Cells
	int numCells = grid->dim.x * grid->dim.y * grid->dim.z;

	//temporary cell list
	GridCell* cellListtmp=new GridCell[numCells];
	for(int i=0;i<numCells;++i){
		//copy length of prismIdxList for each cell
		cellListtmp[i].length=grid->cellList[i].length;
		//allocate mem for prismIdxList for each cell on GPU
		CUDA_CHECK(cudaMalloc((void**) &(cellListtmp[i].prismIdxList), grid->cellList[i].length * sizeof(int)));
		//copy prismIdxList to GPU for each cell
		CUDA_CHECK(cudaMemcpy(cellListtmp[i].prismIdxList, grid->cellList[i].prismIdxList, grid->cellList[i].length * sizeof(int), cudaMemcpyHostToDevice));
	}
	
	//allocate memory for cellList on GPU
	CUDA_CHECK(cudaMalloc((void**) &(d_grid->cellList), numCells * sizeof(GridCell)));
	//copy temporary list to GPU
	CUDA_CHECK(cudaMemcpy(d_grid->cellList, cellListtmp, numCells * sizeof(GridCell), cudaMemcpyHostToDevice));
	
	//copy dim and aabb
	d_grid->dim=grid->dim;
	d_grid->aabb=grid->aabb;
}
