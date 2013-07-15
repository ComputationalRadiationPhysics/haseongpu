#ifndef RAY_PROPAGATION_GPU_KERNEL_H
#define RAY_PROPAGATION_GPU_KERNEL_H
#include <vector>
#include <mesh.h>

float calcDndtAse (unsigned &threads, 
		      unsigned &blocks, 
		      unsigned &hostRaysPerSample,
		      Mesh mesh,
		      Mesh hostMesh,
		      std::vector<double> *betaCellsVector,
		      float nTot,
		      float sigmaA,
		      float sigmaE,
		      float crystalFluorescence,
		      std::vector<double> *dndtAse);


#endif
