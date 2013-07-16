#ifndef calc_sample_phi_ase_H
#define calc_sample_phi_ase_H

#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>
#include <mesh.h>

__global__ void calcSamplePhiAse(curandStateMtgp32* globalState, 
				    const Point samplePoint, 
				    const Mesh mesh,
				    const unsigned* indicesOfPrisms, 
				    const double* importance,
				    const unsigned raysPerSample, 
				    float *phiAse,
				    const unsigned sample_i,
				    const double sigmaA,
				    const double sigmaE,
				    const double nTot);



#endif /* calc_sample_phi_ase_H */
