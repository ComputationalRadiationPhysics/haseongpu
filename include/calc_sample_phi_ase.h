#ifndef calc_sample_phi_ase_H
#define calc_sample_phi_ase_H

#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>
#include <mesh.h>

__global__ void calcSamplePhiAseNew(curandStateMtgp32* globalState, 
				    Point samplePoint, 
				    Mesh mesh,
				    unsigned* indicesOfPrisms, 
				    double* importance,
				    unsigned raysPerSample, 
				    float *phiAse);

__global__ void calcSamplePhiAse(curandStateMtgp32* globalState,
				 float* phiASE,
				 const int point2D,
				 const int level,
				 const int raysPerThread,
				 double *points,
				 double *xOfNormals,
				 double *yOfNormals,
				 int *positionsOfNormalVectors,
				 int *neighbors,
				 int *forbidden,
				 int* triangleIndices,
				 double* betaValues,
				 double* importance,
				 unsigned* indicesOfPrisms,
				 unsigned raysPerSample);

__global__ void setupGlobalVariablesKernel ( 
					    double hostSigmaE,
					    double hostSigmaA, 
					    double hostNTot, 
					    int hostNumberOfTriangles, 
					    double hostThicknessOfPrism, 
					    int hostNumberOfLevels, 
					    int hostNumberOfPoints );


#endif /* calc_sample_phi_ase_H */
