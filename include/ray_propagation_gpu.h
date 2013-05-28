#ifndef RAY_PROPAGATION_GPU_KERNEL_H
#define RAY_PROPAGATION_GPU_KERNEL_H
#include <vector>
#include "curand_kernel.h"
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>


__device__ double cladAbsorption;
__device__ double nTot;
__device__ double sigmaE;
__device__ double sigmaA;
__device__ double thicknessOfPrism;
__device__ int numberOfLevels;
__device__ int cladNumber;
__device__ int numberOfPoints;
__device__ int numberOfTriangles;

__device__ int selectTriangle(int id, int trianglesInOneLevel);
__device__ int selectLevel(int id, int totalNumberOfLevels);
__device__ float rayPropagationGpu(
		double xPos,
		double yPos,
		double zPos,
		double xDestination,
		double yDestination,
		double zDestination,
		int firstTriangle,
		int firstLevel,
		double *points,
		double *xOfNormals,
		double *yOfNormals,
		int *positionsOfNormalVectors,
		int *neighbors,
		int *forbidden,
		int* cellTypes,
		double* betaValues,
		int id);


__global__ void setupGlobalVariablesKernel (
		double host_sigma_e,
		double host_sigma_a,
		int host_clad_num,
		int host_clad_abs,
		double host_N_tot, 
		int host_N_cells, 
		double host_z_mesh, 
		int host_mesh_z, 
		int hostNumberOfPoints );

__global__ void testKernel();

__global__ void raytraceStep( 
		curandStateMtgp32* globalState, 
		float* phi,
		int point2D, 
		int level, 
		int iterations, 
		double *p_in,
		double *n_x,
		double *n_y, 
		int *n_p, 
		int *neighbors, 
		int *forbidden, 
		int* t_in, 
		int* cell_type, 
		int host_size_t, 
		double* beta_v);
