#ifndef RAY_PROPAGATION_GPU_KERNEL_H
#define RAY_PROPAGATION_GPU_KERNEL_H
#include <vector>
#include <mesh.h>

float calcDndtAseNew (unsigned &threads, 
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


float calcDndtAse(
		std::vector<double> *dndtAse, 
		unsigned &threads, 
		unsigned &blocks, 
		unsigned &hostRaysPerSample,
		std::vector<double> *betaValuesVector,
		std::vector<double> *xOfNormalsVector,
		std::vector<double> *yOfNormalsVector,
		std::vector<unsigned> *triangleIndicesVector,
		std::vector<int> *forbiddenVector,
		std::vector<int> *neighborsVector,
		std::vector<int> *positionsOfNormalVectorsVector,
		std::vector<double> *pointsVector,
		std::vector<double> *betaCellsVector,
		std::vector<float> *surfacesVector,
		std::vector<double> *xOfTriangleCenterVector,
		std::vector<double> *yOfTriangleCenterVector,
		float hostNTot,
		float hostSigmaA,
		float hostSigmaE,
		unsigned hostNumberOfPoints,
		unsigned hostNumberOfTriangles,
		unsigned hostNumberOfLevels,
		float hostThicknessOfPrism,
		float hostCrystalFluorescence);

#endif
