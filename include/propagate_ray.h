#ifndef propagate_ray_H
#define propagate_ray_H
__host__ __device__ double propagateRay(
		double xPos,
		double yPos,
		double zPos,
		double xDestination,
		double yDestination,
		double zDestination,
		int firstTriangle,
		int firstLevel,
		double *points,
		double *xOfNormals,
		double *yOfNormals,
		int *positionsOfNormalVectors,
    int *neighbors,
    int *forbidden,
    double* betaValues,
    double nTot,
    double sigmaE,
    double sigmaA,
    double thicknessOfPrism,
    int numberOfLevels,
    int numberOfPoints,
    int numberOfTriangles);

#endif
