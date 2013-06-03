#ifndef RAY_PROPAGATION_GPU_KERNEL_H
#define RAY_PROPAGATION_GPU_KERNEL_H
#include <vector>

float runRayPropagationGpu(
		std::vector<double> *ase, 
		unsigned &threads, 
		unsigned &blocks, 
		unsigned &totalNumberOfRays,
		std::vector<double> *betaValuesVector,
		std::vector<double> *xOfNormalsVector,
		std::vector<double> *yOfNormalsVector,
		std::vector<unsigned> *cellTypesVector,
		std::vector<unsigned> *triangleIndicesVector,
		std::vector<int> *forbiddenVector,
		std::vector<int> *neighborsVector,
		std::vector<int> *positionsOfNormalVectorsVector,
		std::vector<double> *pointsVector,
		std::vector<double> *betaCellsVector,
		float hostCladAbsorption,
		float hostCladNumber,
		float hostNTot,
		float hostSigmaA,
		float hostSigmaE,
		unsigned hostNumberOfPoints,
		unsigned hostNumberOfTriangles,
		unsigned hostNumberOfLevels,
		float hostThicknessOfPrism,
		float hostCrystalFluorescence);

#endif
