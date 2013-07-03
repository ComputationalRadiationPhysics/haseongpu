#include <stdio.h>
#define SMALL 1E-06

double checkSide(int offset, double xN, double yN, double *points, double xPos, double yPos, double length, double xVec, double yVec,int *positionsOfNormalVectors,unsigned numberOfPoints){
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

double checkTop(int levelCurrent, double zPos, double zVec, double length,float thicknessOfPrism){
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

double propagateRayHost(
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
    unsigned numberOfPoints,
    unsigned numberOfTriangles,
    float thicknessOfPrism,
    float sigmaA,
    float sigmaE,
    float nTot){

      double xVec, yVec, zVec;
      double distanceTotal, distanceRemaining, length;
      double gain=1.;
      int triangleCurrent, levelCurrent, forbiddenCurrent; // for the current iteration
      int triangleNext, levelNext, forbiddenNext;		// for the next iteration
      int offset;

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
      xVec /= distanceTotal;
      yVec /= distanceTotal;
      zVec /= distanceTotal;

      // remaining distance to travel
      distanceRemaining = distanceTotal;

      // at the beginning, all surfaces are possible
      forbiddenCurrent = -1;

      while(fabs(distanceRemaining) > SMALL)
      {
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
      }
      return gain /= (distanceTotal * distanceTotal);
    }

/**
 * calculate the gain from the centers of each of the boxes to the observed point
 * calculate the gain and make a "mapping"
 * receipt: pick the point in the center of one cell, 
 * calculate the gain from this point to the observed point,
 * estimate the inner part of the Phi_ASE - Integral,
 * scale the amount of rays proportionally with it
 * sum the amount of rays and scale it to Int=1, which gives the inverse weights
 * the number of rays is determined via floor(), with ceil(), zero-redions could be added
 * use the routine "propagation"!, test: no reflections, just exponential
 *
 **/
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
    float nTot
    )
{
  int raysLeft;
  unsigned raysDump;
  double sumPhi;
  double surfaceTotal;
  double xPos, yPos, zPos;
  double prop;

  raysDump = 0;
  sumPhi = 0.0;
  surfaceTotal = 0.0;
  xPos = points[point];
  yPos = points[point + numberOfPoints];
  zPos = startLevel * thicknessOfPrism;

  // Calculate importance by propagation from trianglecenter to every other center
  for (int i_t=0; i_t < numberOfTriangles; ++i_t){
    for (int i_z=0; i_z < (numberOfLevels-1); ++i_z){
      prop = propagateRayHost(xOfTriangleCenter[i_t], yOfTriangleCenter[i_t], 
          thicknessOfPrism * (i_z+0.5),  xPos, yPos, zPos, i_t, i_z, 
          points, xOfNormals, yOfNormals, positionsOfNormalVectors, 
          neighbors, forbidden, betaValues,
          numberOfPoints, numberOfTriangles, thicknessOfPrism,
          sigmaA, sigmaE, nTot
          );

      importance[i_t + i_z * numberOfTriangles] = betaValues[i_t + i_z * numberOfTriangles]*(prop);
      sumPhi += importance[i_t + i_z * numberOfTriangles];

    }
    surfaceTotal += surface[i_t];

  }

  // Calculate number of rays/prism
  for (int i_t=0; i_t < numberOfTriangles; ++i_t){
    for (int i_z=0; i_z < (numberOfLevels-1); ++i_z){
      numberOfImportantRays[i_t + i_z * numberOfTriangles] = (unsigned)(floor(importance[i_t + i_z * numberOfTriangles] / sumPhi * raysPerSample));
      raysDump +=  numberOfImportantRays[i_t + i_z * numberOfTriangles];
    }

  }
  raysLeft = raysPerSample - raysDump;

  // TODO What happens with random failure ?
  // Distribute the remaining rays randomly
  //  for (int i_r=0; i_r < raysLeft; i_r++){
  //    int rand_t = (int )(rand() % numberOfTriangles);
  //    int rand_z = (int )(rand() % (numberOfLevels-1));
  //    numberOfImportantRays[rand_t + rand_z * numberOfTriangles]++;
  //
  //  }

  //  Now think about the mount of rays which would come out of this volume(surface)
  //  dividing this number with the new amount of rays gives the final importance weight for this area!
  for (int i_t=0; i_t<numberOfTriangles; ++i_t){
    for (int i_z=0; i_z<(numberOfLevels-1); ++i_z){
      if (numberOfImportantRays[i_t + i_z*numberOfTriangles] > 0){
        importance[i_t + i_z*numberOfTriangles] = raysPerSample * surface[i_t] / surfaceTotal / numberOfImportantRays[i_t + i_z*numberOfTriangles];
      }
      else{
        importance[i_t + i_z*numberOfTriangles] = 0; 
      }
    }
  }
  return raysDump;
}
