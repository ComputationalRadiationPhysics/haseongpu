#include <mesh.h>
#include <assert.h>
#include <stdio.h>
#include <geometry.h>
#include <cuda_runtime_api.h>
#include <stdio.h> /* printf */

#define TEST_VALUES true
#define MAX_LEVEL 10

/*##########################################
  # RECONSTRUCTION                         #
  ##########################################*/

__host__ __device__ double checkSurface(const int currentLevel, const double zPos, const double zVec, const double length, const double thickness){
  double denominator = zVec;
  if (denominator != 0.0){
    double nominator = currentLevel * thickness - zPos;
    double lengthTmp = nominator/denominator;
    // DEBUG
    //printf("C lengthTmp : %f nominator: %f denominaor %f level %d\n", lengthTmp, nominator, denominator, currentLevel);
    if (lengthTmp <= length && lengthTmp > 0.0){
      return lengthTmp;
    }

  }
  return 0;
}

__host__ __device__ double checkEdge(const Triangle triangle, const int edge, const Ray ray, const double length){
  double denominator = triangle.edges[edge].normal.dir.x * ray.dir.x + triangle.edges[edge].normal.dir.y * ray.dir.y;
  if (denominator != 0.0)
    {
      double nominator =	  
	triangle.edges[edge].normal.dir.x * triangle.edges[edge].normal.p.x
	+ triangle.edges[edge].normal.dir.y * triangle.edges[edge].normal.p.y
	- triangle.edges[edge].normal.dir.x * ray.p.x 
	- triangle.edges[edge].normal.dir.y * ray.p.y; 

      double lengthTmp = nominator/denominator;
      // DEBUG
      //printf("C lengthTmp : %f length: %f nominator: %f denominaor %f\n", lengthTmp, length, nominator, denominator);
      if(lengthTmp <= length && lengthTmp > 0.0){
	return lengthTmp;
      }

    }
  else{
    printf("denominaor %f : %f %f %f %f\n", denominator, triangle.edges[edge].normal.dir.x, ray.dir.x, triangle.edges[edge].normal.dir.y, ray.dir.y);
  }

  return 0;
}

__host__ __device__ int getNextEdge(const Triangle triangle,  const Ray ray, const unsigned level, double length, const int forbiddenEdge, const double thickness){
  assert(forbiddenEdge >= -1 && forbiddenEdge <= 4);
  int edge = -1;
  // Check 3 edges of triangle
  for(int edge_i = 0; edge_i < 3; ++edge_i){
    if(edge_i != forbiddenEdge){
      double lengthTmp = checkEdge(triangle, edge_i, ray, length);
      if(lengthTmp){
	length = lengthTmp;
	edge = edge_i;
      }
    }
  }
  
  // check the upper surface
  if (forbiddenEdge != 3){
    double lengthTmp = checkSurface(level + 1, ray.p.z, ray.dir.z, length, thickness);
    if(lengthTmp){
      length = lengthTmp;
      edge = 3;
    }
  }

  // check the lower surface
  if (forbiddenEdge != 4){
    double lengthTmp = checkSurface(level, ray.p.z, ray.dir.z, length, thickness);
    if (lengthTmp){
      length = lengthTmp;
      edge = 4;
    }
  }
  return edge;
}

__host__ __device__ unsigned getNextLevel(const unsigned level, const int edge){
  switch(edge){
  case 3:
    return level + 1;
  case 4:
    if(level == 0)
      return level;
    return level - 1;
  default:
    return level;
  }

}

__host__ __device__ double calcTriangleIntersection(const Triangle triangle, const Ray ray, const int edge, const double length, const unsigned level, const double thickness, const double distanceRemaining){
  switch(edge){
  case -1:
    return distanceRemaining;
  case 0:
  case 1:
  case 2:
    return checkEdge(triangle, edge, ray, length);
  case 3:
    return checkSurface(level + 1, ray.p.z, ray.dir.z, length, thickness);
  case 4:
    return checkSurface(level, ray.p.z, ray.dir.z, length, thickness);

  }
  return 0;

}

__host__ __device__ Ray calcNextRay(Ray ray, const double length){
  ray.p.x = ray.p.x + length * ray.dir.x;
  ray.p.y = ray.p.y + length * ray.dir.y;
  ray.p.z = ray.p.z + length * ray.dir.z;
  ray.length = ray.length - length;

  return ray;

}

__host__ __device__ double calcPrismGain(const Triangle triangle, const unsigned level, const double length, const double sigmaA, const double sigmaE, const double nTot){
  assert(level < MAX_LEVEL);
  return (double) exp(nTot * (triangle.betaValues[level] * ( sigmaE + sigmaA ) - sigmaA ) * length);
 
}

__host__ __device__ Triangle getNextTriangle(Triangle triangle, int edge){

  if(edge == 3 || edge == 4)
    return triangle;
  return *(triangle.edges[edge].neighbor);

}

__host__ __device__ int getNextForbiddenEdge(Triangle triangle, int edge){
 switch(edge){
  case 0:
  case 1:
  case 2:
    return triangle.edges[edge].forbidden;
  case 3:
    return 4;
  case 4:
    return 3;

  }
  return 0;

}

__host__ __device__ double propagateRayNew(Ray ray, unsigned startLevel, Triangle startTriangle, const double sigmaA, const double sigmaE, const double nTot, const double thickness){
  double distanceTotal = ray.length;
  double distanceRemaining = distanceTotal;
  double length = 0;
  double gain = 1;
  
  Triangle nextTriangle = startTriangle;
  Ray nextRay = normalizeRay(ray, distanceTotal);
  int nextForbiddenEdge = -1;
  int nextEdge = -1;
  unsigned nextLevel = startLevel;

  unsigned debugLoopCount = 0;
  
  while(fabs(distanceRemaining) > SMALL){
    // DEBUG
    //    printf("\nC Calc next triangle, distanceRemaining: %f nextForbiddenEdge: %d nextLevel: %d\n", distanceRemaining, nextForbiddenEdge, nextLevel);
    //    printf("C nextRay : POS(%f,%f,%f) VEC(%f, %f, %f)\n", nextRay.p.x, nextRay.p.y, nextRay.p.z, nextRay.dir.x, nextRay.dir.y, nextRay.dir.z);
    //    printf("C nextTriangle: A(%f,%f), B(%f,%f), C(%f,%f)\n", 
    //    	     nextTriangle.A.x, nextTriangle.A.y,
    //    	     nextTriangle.B.x, nextTriangle.B.y,
    //    	     nextTriangle.C.x, nextTriangle.C.y);
    debugLoopCount++;
    assert(debugLoopCount <= 1000);
    nextEdge            = getNextEdge(nextTriangle, nextRay, nextLevel, distanceRemaining, nextForbiddenEdge, thickness);
    length              = calcTriangleIntersection(nextTriangle, nextRay, nextEdge, distanceRemaining, nextLevel, thickness, distanceRemaining);
    gain               *= calcPrismGain(nextTriangle, nextLevel, length, sigmaA, sigmaE, nTot);
    distanceRemaining  -= length;

    // Calc structs for next step
    if(nextEdge != -1){
      nextLevel         = getNextLevel(nextLevel, nextEdge);
      nextRay           = calcNextRay(nextRay, length);
      nextForbiddenEdge = getNextForbiddenEdge(nextTriangle, nextEdge);
      nextTriangle      = getNextTriangle(nextTriangle, nextEdge);
    }

  }

  return gain /= (distanceTotal * distanceTotal);
}
 
/* ############################################
   # OLD CODE                                 #
   ############################################*/

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
		  if (forbiddenCurrent != (int)i){
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
