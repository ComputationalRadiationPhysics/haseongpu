#ifndef importance_sampling_H
#define importance_sampling_H

#include <mesh.h>

void importanceSamplingNew(
			   Point samplePoint, 
			   Mesh mesh,
			   unsigned raysPerSample,
			   double sigmaA,
			   double sigmaE,
			   double nTot,
			   double *importance,
			   unsigned *raysPerPrism);

unsigned importanceSampling(int point,
	     int startLevel,
	     double *importance,
	     unsigned *numberOfImportantRays,
	     double *points,
	     double *xOfNormals,
	     double *yOfNormals,
	     int *positionsOfNormalVectors,
	     int *neighbors,
	     int *forbidden,
	     double *betaValues,
	     double *xOfTriangleCenter,
	     double *yOfTriangleCenter,
	     float *surface,
	     unsigned raysPerSample,
	     unsigned numberOfPoints,
	     unsigned numberOfLevels,
	     unsigned numberOfTriangles,
	     float thicknessOfPrism,
	     float sigmaA,
	     float sigmaE,
	     float nTot);

#endif /* importance_sampling_H */
