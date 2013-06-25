#ifndef test_functions_H
#define test_functions_H

__global__ void testKernel (
			    double *points,
			    double *xOfNormals,
			    double *yOfNormals,
			    int *neighbors,
			    int *forbidden,
			    int *positionsOfNormalVectors,
			    int* triangleIndices,
			    double* betaValues,
					float *phiAse,
					double *importance,
					unsigned *indicesOfPrisms,
					float hostNTot,
					float hostSigmaA,
					float hostSigmaE,
					unsigned hostNumberOfPoints,
					unsigned hostNumberOfTriangles,
					unsigned hostNumberOfLevels,
					float hostThicknessOfPrism,
					float hostCrystalFluorescence,
					unsigned limit
					);


#endif /* test_functions_H */
