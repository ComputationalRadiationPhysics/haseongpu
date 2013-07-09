#ifndef importance_sampling_H
#define importance_sampling_H

<<<<<<< HEAD
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

void importanceSampling(
			int point,
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
=======
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
>>>>>>> c872b097b14330c8dd939cf52fada8582d7015d6
#endif /* importance_sampling_H */
