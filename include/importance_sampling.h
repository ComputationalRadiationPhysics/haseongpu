#ifndef importance_sampling_H
#define importance_sampling_H

#include <mesh.h>

unsigned importanceSampling(
    Point samplePoint, 
    Mesh deviceMesh,
    unsigned raysPerSample, 
    double sigmaA, 
    double sigmaE, 
    double nTot,  
    double *importance, 
    float *sumPhi,
    unsigned *raysPerPrism,
    unsigned *indicesOfPrisms,
    unsigned *raysDump,
    int threads,
    int blocks);

#endif /* importance_sampling_H */
