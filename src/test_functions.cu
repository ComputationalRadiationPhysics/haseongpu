#include <stdio.h>
/**
 * Prints some of the global device variables.
 * Is only used for testing
 */
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
			    float *surfaces){
  // printf("\nSigmaE=%.6e",sigmaE);
  // printf("\nSigmaA=%.6e",sigmaA);
  // printf("\nNumberOfLevels=%d",numberOfLevels);
  // printf("\nNumberOfPoints=%d",numberOfPoints);
  // printf("\nthicknessOfPrism_=%.6e",thicknessOfPrism);
  // printf("\nnumberOfTriangles=%d",numberOfTriangles);
  // printf("\nnTot=%.6e",nTot);
  // printf("\ncladAbsorption=%.6e",cladAbsorption);
  // printf("\ncladNumber=%d\n\n",cladNumber);

  unsigned limit = 5;
  for(int i=0;i<limit;++i){
    printf("points[%d]: %e\n",i,points[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("xOfNormals[%d]: %e\n",i,xOfNormals[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("yOfNormals[%d]: %e\n",i,yOfNormals[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("positionsOfNormalVectors[%d]: %d\n",i,positionsOfNormalVectors[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("neighbors[%d]: %d\n",i,neighbors[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("forbidden[%d]: %d\n",i,forbidden[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("triangleIndices[%d]: %d\n",i,triangleIndices[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("cellTypes[%d]: %d\n",i,cellTypes[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("betaValues[%d]: %e\n",i,betaValues[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("surfaces[%d]: %f\n",i,surfaces[i]);
  }
  printf("\n\n");
} 

