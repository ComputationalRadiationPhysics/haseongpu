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

__global__ void calcSamplePhiAse(
    curandStateMtgp32* globalState,
    float* phiASE,
    int point2D,
    int level,
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
    unsigned raysPerSample,
    double nTot,
    double sigmaE,
    double sigmaA,
    double thicknessOfPrism,
    int numberOfLevels,
    int numberOfPoints,
    int numberOfTriangles);



#endif /* calc_sample_phi_ase_H */
