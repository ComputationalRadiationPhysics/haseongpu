#include <stdio.h> /* printf */

#define TEST_VALUES true
#define SMALL 1E-06


__host__ __device__ double checkSide(int offset, double xN, double yN, double *points, double xPos, double yPos, double length, double xVec, double yVec,int *positionsOfNormalVectors,int numberOfPoints){
	double denominator = xN * xVec + yN * yVec;
	if (denominator != 0.0)
	{
		double nominator =	  xN * points[ positionsOfNormalVectors[offset] ]
							+ yN * points[ positionsOfNormalVectors[offset] + numberOfPoints ]
							- xN * xPos 
							- yN * yPos;

		double lengthHelp = nominator/denominator;

		if(lengthHelp < length && lengthHelp > 0.0){
			return lengthHelp;
		}
	}
	return 0;
}

__host__ __device__ double checkTop(int levelCurrent, double zPos, double zVec, double length,double thicknessOfPrism){
	double denominator = zVec;
	if (denominator != 0.0){
		double nominator = levelCurrent * thicknessOfPrism - zPos;
		double lengthHelp = nominator/denominator;
		if (lengthHelp < length && lengthHelp > 0.0){
			return lengthHelp;
		}
	}
	return 0;
}
/**
 * @brief Propagate a ray between 2 points and calculate the resulting ASE-Flux at the Destination
 *
 * @params x_pos		the x-coordinate where the ray starts
 *         y_pos		the y-coordinate where the ray starts
 *         z_pos		the z-coordinate where the ray starts
 *         x_dest		the destination of the ray (x-coordinate)
 *         y_dest		the destination of the ray (y-coordinate)
 *         z_dest		the destination of the ray (z-coordinate)
 *         t_start		the index of the triangle, in which the ray starts
 *         mesh_start	the level of the mesh (the slice) in which the ray starts
 *		   p_in			coordinates of the sample-points of one layer (first all x-coordinates, then all y-coordinates)
 *		   n_x			x-coordinates for the normal-vectors for the 3 rectangular sides of each prism
 *		   n_y			y-coordinates for the normal-vectors for the 3 rectangular sides of each prism
 *		   n_p			indices of the points where the normal-vectors start	
 *		   neighbors	indices of the adjacent triangles	
 *		   forbidden	sides of the new triangles which are forbidden, after coming from an adjacent triangle
 *		   cell_type	contains the material-constant for each cell/prism
 *		   beta_v		contains the beta-values for each cell/prism
 *
 */
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
    int numberOfTriangles){

	double xVec, yVec, zVec;
	double distanceTotal, distanceRemaining, length;
	double gain=1.;
	int triangleCurrent, levelCurrent, forbiddenCurrent; // for the current iteration
	int triangleNext, levelNext, forbiddenNext;		// for the next iteration
	int offset;
#if TEST_VALUES==true
	double testDistance = 0;
	int loopbreaker = 0;
#endif

	// initial positions
	triangleCurrent = firstTriangle;
	levelCurrent = firstLevel;

	// direction-vector (without reflections)
	xVec = (xDestination - xPos);
	yVec = (yDestination - yPos);
	zVec = (zDestination - zPos);

	// total distance to travel
	distanceTotal = sqrt(xVec*xVec + yVec*yVec + zVec*zVec);

	// normalized direction-vector
	xVec = xVec / distanceTotal;
	yVec = yVec / distanceTotal;
	zVec = zVec / distanceTotal;

	// remaining distance to travel
	distanceRemaining = distanceTotal;

	// at the beginning, all surfaces are possible
	forbiddenCurrent = -1;

	while(fabs(distanceRemaining) > SMALL)
	{
		// the length of the ray-part inside the current prism. We try to minimize this value
		length = distanceRemaining;
		levelNext = levelCurrent;

		// check the 3 edges
		for (unsigned i = 0; i<3 ; ++i){
			if (forbiddenCurrent != i){
				offset = triangleCurrent + i * numberOfTriangles;
				double lengthHelp = checkSide(offset, xOfNormals[offset], yOfNormals[offset], points, xPos, yPos, length, xVec, yVec, positionsOfNormalVectors,numberOfPoints);
				if(lengthHelp){
					length = lengthHelp;
					forbiddenNext = (forbidden[offset]);
					triangleNext = neighbors[offset];
				}
			}
		}

		// check the upper plane
		if (forbiddenCurrent != 3){
			double lengthHelp = checkTop(levelCurrent+1, zPos, zVec, length,thicknessOfPrism);
			if(lengthHelp){
				length = lengthHelp;
				forbiddenNext = 4; // you are not allowed to go down in the next step
				triangleNext = triangleCurrent;
				levelNext = levelCurrent + 1;
			}
		}

		// check the lower plane
		if (forbiddenCurrent != 4){
			double lengthHelp = checkTop(levelCurrent, zPos, zVec, length,thicknessOfPrism);
			if (lengthHelp){
				length = lengthHelp;
				forbiddenNext = 3; // you are not allowed to go up in the next step
				triangleNext = triangleCurrent;
				levelNext = levelCurrent - 1;
			}
		}

		gain *= (double) exp(nTot * (betaValues[triangleCurrent + (levelCurrent * numberOfTriangles)] * (sigmaE + sigmaA) - sigmaA) * length);

		// the remaining distance is decreased by the length we travelled through the prism
		distanceRemaining -= length;

		// now set the next cell and position
		xPos += length * xVec;
		yPos += length * yVec;
		zPos += length * zVec;

		triangleCurrent = triangleNext;
		levelCurrent = levelNext;
		forbiddenCurrent = forbiddenNext;

#if TEST_VALUES==true
		testDistance += length;
		if(loopbreaker>500){
			printf("Loopbreaker reached.\n");
			return 0.;
		}else{
			loopbreaker++;
		}
#endif
	}

#if TEST_VALUES==true
	if(fabs(distanceTotal - testDistance) > SMALL)
		printf("Distance wrong!\n");
	if(fabs(xPos - xDestination) > SMALL) printf("X Coordinate wrong!\n");
	if(fabs(yPos - yDestination) > SMALL) printf("Y Coordinate wrong!\n");
	if(fabs(zPos - zDestination) > SMALL) printf("Z Coordinate wrong!\n");
#endif

	return gain /= (distanceTotal * distanceTotal);
}
