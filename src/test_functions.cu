#include "test_functions.h"
#include <stdio.h>
/**
 * Prints some of the global device variables.
 * Is only used for testing
 */
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
		float	hostCrystalFluorescence,
		unsigned limit
){
  printf("\nSigmaE=%.6e",hostSigmaE);
  printf("\nSigmaA=%.6e",hostSigmaA);
  printf("\nNumberOfLevels=%d",hostNumberOfLevels);
  printf("\nNumberOfPoints=%d",hostNumberOfPoints);
  printf("\nthicknessOfPrism_=%.6e",hostThicknessOfPrism);
  printf("\nnumberOfTriangles=%d",hostNumberOfTriangles);
	printf("\ncrystalFluorescence=%.6e",hostCrystalFluorescence);
  printf("\nnTot=%.6e\n\n",hostNTot);

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
    printf("neighbors[%d]: %d\n",i,neighbors[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("forbidden[%d]: %d\n",i,forbidden[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("positionsOfNormalVectors[%d]: %d\n",i,positionsOfNormalVectors[i]);
  }
  printf("\n\n");


  for(int i=0;i<limit;++i){
    printf("triangleIndices[%d]: %d\n",i,triangleIndices[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("betaValues[%d]: %e\n",i,betaValues[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("phiAse[%d]: %.f\n",i,phiAse[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("importance[%d]:	%f\n",i,importance[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("indicesOfPrisms[%d]: %d\n",i,indicesOfPrisms[i]);
  }
  printf("\n\n");
} 

