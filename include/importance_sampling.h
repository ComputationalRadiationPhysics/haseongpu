#ifndef importance_sampling_H
#define importance_sampling_H

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
