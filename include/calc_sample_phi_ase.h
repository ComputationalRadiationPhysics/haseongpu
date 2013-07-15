#ifndef calc_sample_phi_ase_H
#define calc_sample_phi_ase_H

#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>
#include <mesh.h>

__global__ void calcSamplePhiAse(curandStateMtgp32* globalState, 
				    Point samplePoint, 
				    Mesh mesh,
				    unsigned* indicesOfPrisms, 
				    double* importance,
				    unsigned raysPerSample, 
				    float *phiAse,
				    unsigned sample_i,
				    const double sigmaA,
				    const double sigmaE,
				    const double nTot);



#endif /* calc_sample_phi_ase_H */
