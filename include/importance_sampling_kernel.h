#ifndef importance_sampling_kernel_H
#define importance_sampling_kernel_H

__global__ void importanceKernel1(
    int point,
    int level,
    double *importance,
    double *points,
    double *xOfNormals,
    double *yOfNormals,
    int *positionsOfNormalVectors,
    int *neighbors,
    int *forbidden,
    double *betaValues,
    double *xOfTriangleCenter,
    double *yOfTriangleCenter,
    unsigned numberOfPoints,
    unsigned numberOfLevels,
    unsigned numberOfTriangles,
    float thicknessOfPrism,
    float sigmaA,
    float sigmaE,
    float nTot);



__global__ void importanceKernel2(
    unsigned *numberOfImportantRays,
    double *importance,
    float sumPhi,
    int raysPerSample,
    unsigned numberOfPrisms);



__global__ void importanceKernel3(
    int raysPerSample,
    int raysDump,
    unsigned *numberOfImportantRays,
    unsigned numberOfLevels,
    unsigned numberOfTriangles);


__global__ void importanceKernel4(
    unsigned *numberOfImportantRays,
    double *importance,
    float *surface,
    int surfaceTotal,
    unsigned raysPerSample,
    unsigned numberOfPrisms,
    unsigned numberOfTriangles);



#endif /* importance_sampling_H */
