#ifndef test_functions_H
#define test_functions_H

__global__ void testKernel (
			    double *points,
			    double *xOfNormals,
			    double *yOfNormals,
			    int *positionsOfNormalVectors,
			    int *neighbors,
			    int *forbidden,
			    int* triangleIndices,
			    int* cellTypes,
			    double* betaValues,
			    float *surfaces);

#endif /* test_functions_H */
