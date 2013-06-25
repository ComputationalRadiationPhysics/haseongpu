#include <curand_kernel.h> /* curand_uniform */
#include <stdio.h> /* printf */
#include <mesh.h>

#define TEST_VALUES true
#define SMALL 1E-06

__device__ double nTot;
__device__ double sigmaE;
__device__ double sigmaA;
__device__ double thicknessOfPrism;
__device__ int numberOfLevels;
__device__ int numberOfPoints;
__device__ int numberOfTriangles;


__device__ int getNextEdge(Triangle triangle,  Point point, Point destinationPoint, int forbiddenEdge){
  return 0;

}

__device__ double calcTriangleIntersection(Triangle triangle, Point point, Point destinationPoint, Edge edge){
  return 0;

}

__device__ Point calcNextPoint(Point point, Point destinationPoint, double length){
  return {0,0,0,0};

}

__device__ double calcPrismGain(Triangle triangle, Point edgePoint, double length){
  // return (double) exp(nTot * (triangle.betaValues[(int)edgePoint[2]] * ( sigmaE + sigmaA ) - sigmaA ) * length);
  return 0;

}

__device__ double propagateRayDeviceNew(normalRay ray, Triangle startTriangle, Triangle *triangles){
  double distanceRemaining = 0;
  double distanceTotal = ray.length;
  double length = 0;
  double gain = 1;

  Triangle nextTriangle = startTriangle;
  Point nextPoint = ray.point;
  int nextForbiddenEdge = -1;
  int nextEdge = -1;

  distanceRemaining = distanceTotal;
  
  while(fabs(distanceRemaining) < SMALL){
    nextEdge          = getNextEdge(nextTriangle,  ray, nextForbiddenEdge);
    length            = calcTriangleIntersection(nextTriangle, nextPoint, destinationPoint, nextTriangle.edges[nextEdge]);
    nextPoint         = calcNextPoint(nextPoint, destinationPoint, length);
    nextForbiddenEdge = nextTriangle.edges[nextEdge].forbidden;
    nextTriangle      = *(nextTriangle.edges[nextEdge].neighbor);

    gain *= calcPrismGain(nextTriangle, nextPoint, length);
    distanceRemaining -= length;

  }

  return gain /= (distanceTotal * distanceTotal);

};



__device__ double checkSide(int offset, double xN, double yN, double *points, double xPos, double yPos, double length, double xVec, double yVec,int *positionsOfNormalVectors){
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

__device__ double checkTop(int levelCurrent, double zPos, double zVec, double length){
	double denominator = zPos*zVec;
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
__device__ double propagateRayDevice(
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
		double* betaValues){

	double xVec, yVec,zVec;
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
		//        definition for decider
		//        0,1,2: int for the neighbors
		//        3: hor plane up
		//        4: hor plane down
		//        try the triangle faces
		//        remember the correlation between the normals and the points
		//        n1: p1-2, n2: p1-3, n3:p2-3
		//        the third coordinate (z) of the particpating points for the surfaces can be set to be z=0, 
		//        as everything uses triangular "prisms", as well as n_z=0 in this case!

		// forb describes the surface, from which the ray enters the prism.
		// this surface is no suitable candidate, since the length would be 0!

		levelNext = levelCurrent;

		// check the 3 edges
		for (int i = 0; i<3 ; ++i){
			if (forbiddenCurrent != i){
				offset = triangleCurrent + i * numberOfTriangles;
				double lengthHelp = checkSide(offset, xOfNormals[offset], yOfNormals[offset], points, xPos, yPos, length, xVec, yVec, positionsOfNormalVectors);
				if(lengthHelp){
					length = lengthHelp;
					forbiddenNext = (forbidden[offset]);
					triangleNext = neighbors[offset];
				}
			}
		}

		// check the upper plane
		if (forbiddenCurrent != 3){
			double lengthHelp = checkTop(levelCurrent+1, zPos, zVec, length);
			if(lengthHelp){
				length = lengthHelp;
				forbiddenNext = 4; // you are not allowed to go down in the next step
				triangleNext = triangleCurrent;
				levelNext = levelCurrent + 1;
			}
		}

		// check the lower plane
		if (forbiddenCurrent != 4){
			double lengthHelp = checkTop(levelCurrent, zPos, zVec, length);
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


#if TEST_VALUES==true
		testDistance += length;
		if(loopbreaker>500){
			printf("Loopbreaker reached. firstTriangle: %d, level: %d, length: %f, distanceTotal:%f, testDistance%f, distanceRemaining:%f\n",firstTriangle,firstLevel,length,distanceTotal,testDistance,distanceRemaining);
			return 0.;
		}else{
			loopbreaker++;
		}
#endif

		// now set the next cell and position
		xPos = xPos + length*xVec;
		yPos = yPos + length*yVec;
		zPos = zPos + length*zVec;

		triangleCurrent = triangleNext;
		levelCurrent = levelNext;
		// set the new forbidden surface
		forbiddenCurrent = forbiddenNext;

	}

#if TEST_VALUES==true
	if(fabs(distanceTotal - testDistance) > SMALL)
		printf("Distance too big! firstTriangle: %d, level: %d, length: %f, distanceTotal:%f, testDistance%f, distanceRemaining:%f\n",firstTriangle,firstLevel,length,distanceTotal,testDistance,distanceRemaining);
#endif

	return gain /= (distanceTotal * distanceTotal);
}



/**
 * Initializes the global variables of the GPU with the correct values.
 * All those values are from the original propagation-function which we ported.
 */
__global__ void setupGlobalVariablesKernel ( 
					    double hostSigmaE,
					    double hostSigmaA, 
					    double hostNTot, 
					    int hostNumberOfTriangles, 
					    double hostThicknessOfPrism, 
					    int hostNumberOfLevels, 
					    int hostNumberOfPoints )
{
  sigmaE = hostSigmaE;	
  sigmaA = hostSigmaA;
  nTot = hostNTot;
  numberOfTriangles = hostNumberOfTriangles;
  thicknessOfPrism = hostThicknessOfPrism;
  numberOfLevels = hostNumberOfLevels;
  numberOfPoints = hostNumberOfPoints;
} 


/**
 * Does the raytracing for a single Sample point (in a defined level).
 * This Kernel has to be started for each sample point with the same value for iterations
 * and the same number of blocks/threads.
 *
 * \var globalState the state of the mersenneTwister PRNG
 * 		(has a maximum of 200 positions!)
 * \var phi points to a memory region which is initialized with 0
 * 		(can hold one value for each sample point)
 * \var point2D the index of the current sample point (points to p_in)
 * \var level the level of the current sample point (how deep we are through the material)
 * \var raysPerThread the number rays which are computed by this thread
 * 		(always for the same combination of startprism+samplepoint
 */
__global__ void calcSamplePhiAse(curandStateMtgp32* globalState,
			     float* phiASE,
			     const int point2D,
			     const int level,
			     const int raysPerThread,
			     double *points,
			     double *xOfNormals,
			     double *yOfNormals,
			     int *positionsOfNormalVectors,
			     int *neighbors,
			     int *forbidden,
			     int* triangleIndices,
			     double* betaValues,
			     double* importance,
			     unsigned* indicesOfPrisms,
			     unsigned raysPerSample) {

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  const int endPointX = points[point2D];
  const int endPointY = points[numberOfPoints + point2D];
  const int endPointZ = level * thicknessOfPrism;


  // on thread can compute multiple rays
  for (int i=0; ; ++i){

	  // the current ray which we compute is based on the id and an offset (number of threads*blocks)
	  int rayNumber = id + (blockDim.x*gridDim.x * i);
	  if(rayNumber >= raysPerSample){
		  return;
	  }

	  // get a new prism to start from
	  int startPrism = indicesOfPrisms[rayNumber];
	  int startLevel = startPrism/numberOfTriangles;
	  int startTriangle = startPrism - (numberOfTriangles * startLevel);

#if TEST_VALUES==true
	  if(startPrism != (startTriangle + (startLevel * numberOfTriangles))){
		  printf("StartTriangle/StartLevel incorrect!");
	  }
	  if(startTriangle >= 600){
		  printf("StartTriangle/StartLevel incorrect!");
	  }
#endif

	  // Get triangle vertex indicies
	  int t1 = triangleIndices[startTriangle];
	  int t2 = triangleIndices[startTriangle + numberOfTriangles];
	  int t3 = triangleIndices[startTriangle + 2 * numberOfTriangles];

	  // random startpoint generation
	  double u = curand_uniform(&globalState[blockIdx.x]);
	  double v = curand_uniform(&globalState[blockIdx.x]);

	  if((u+v)>1)
	  {
		  u = 1-u;
		  v = 1-v;
	  }
	  double w = 1-u-v;

	  // convert the random startpoint into coordinates
	  double xRand = (points[t1] * u) + (points[t2] * v) + (points[t3] * w);
	  double yRand = (points[numberOfPoints + t1] * u) + (points[numberOfPoints + t2] * v) + (points[numberOfPoints + t3] * w);
	  double zRand = (startLevel + curand_uniform(&globalState[blockIdx.x])) * thicknessOfPrism;

	  __syncthreads();
	  // propagate the ray
	  double gain = propagateRayDevice(xRand, yRand, zRand, endPointX, endPointY, endPointZ, 
				   startTriangle, startLevel, points, xOfNormals, yOfNormals, 
				   positionsOfNormalVectors, neighbors, forbidden,  betaValues);

	  gain *= betaValues[startPrism];
	  gain *= importance[startPrism];

	  atomicAdd(&(phiASE[point2D + (level * numberOfPoints)]), float(gain));
  }
}
