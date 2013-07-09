#include <mesh.h>
#include <stdio.h>
#include <propagate_ray.h>
#include <geometry.h>
#include <assert.h>

// ##############################################################
// # Reconstruction                                             #
// ##############################################################
void importanceSamplingNew(Point samplePoint, Mesh mesh, unsigned raysPerSample, double sigmaA, double sigmaE, double nTot,  double *importance, unsigned *raysPerPrism){
  Triangle *triangles = mesh.triangles;

  int raysLeft = 0;
  int raysDump = 0;
  double sumPhi = 0;
  double surfaceTotal = 0;
  double gain = 0;
  Ray ray;
  Point startPoint;
  Triangle startTriangle;

  // Calculate importance by propagation from trianglecenter to every other center
  for(unsigned triangle_i = 0; triangle_i < mesh.numberOfTriangles; ++triangle_i){
    for(unsigned level_i = 0; level_i < mesh.numberOfLevels - 1; ++level_i){
      startTriangle = triangles[triangle_i];
      startPoint.x = startTriangle.center.x;
      startPoint.y = startTriangle.center.y;
      startPoint.z = (level_i + 0.5) * mesh.thickness;
      ray = generateRay(startPoint, samplePoint);

      gain = propagateRayNew(ray, level_i, startTriangle, triangles, sigmaA, sigmaE, nTot, mesh.thickness);

      importance[triangle_i + level_i * mesh.numberOfTriangles] = startTriangle.betaValues[level_i] * gain;
      sumPhi += importance[triangle_i + level_i * mesh.numberOfTriangles];

    }
    surfaceTotal += triangles[triangle_i].surface;
  }

  // Calculate number of rays/prism
  for(unsigned triangle_i = 0; triangle_i < mesh.numberOfTriangles; ++triangle_i){
    for(unsigned level_i = 0; level_i < mesh.numberOfLevels - 1; ++level_i){
      raysPerPrism[triangle_i + level_i * mesh.numberOfTriangles] =  (unsigned)(floor(importance[triangle_i + level_i * mesh.numberOfTriangles] / sumPhi * raysPerSample));
      raysDump +=  raysPerPrism[triangle_i + level_i * mesh.numberOfTriangles];
    }

  }

  raysLeft = raysPerSample - raysDump;

  // TODO What happens with random failure ?
  // TODO Distribute the remaining rays randomly
  for (int i_r=0; i_r < raysLeft; i_r++){
    int rand_t = (int )(rand() % mesh.numberOfTriangles);
    int rand_z = (int )(rand() % (mesh.numberOfLevels-1));
    raysPerPrism[rand_t + rand_z * mesh.numberOfTriangles]++;

  }

  // Now think about the mount of rays which would come out of this volume(surface)
  // dividing this number with the new amount of rays gives the final importance weight for this area!
  for (int triangle_i=0; triangle_i < mesh.numberOfTriangles; ++triangle_i){
    for (int level_i=0; level_i < mesh.numberOfLevels - 1; ++level_i){
      if (raysPerPrism[triangle_i + level_i * mesh.numberOfTriangles] > 0){
  	importance[triangle_i + level_i * mesh.numberOfTriangles] = raysPerSample * triangles[triangle_i].surface / surfaceTotal / raysPerPrism[triangle_i + level_i * mesh.numberOfTriangles];

      }
      else{
  	importance[triangle_i + level_i * mesh.numberOfTriangles] = 0; 

      }

    }

  }

}


// #################################################
// # Old Code                                      #
// #################################################


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
      prop = propagateRay(xOfTriangleCenter[i_t], yOfTriangleCenter[i_t], 
          thicknessOfPrism * (i_z+0.5),  xPos, yPos, zPos, i_t, i_z, 
          points, xOfNormals, yOfNormals, positionsOfNormalVectors, 
          neighbors, forbidden, betaValues,
          nTot, sigmaE, sigmaA, thicknessOfPrism, numberOfLevels, numberOfPoints, numberOfTriangles);

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
    for (int i_r=0; i_r < raysLeft; i_r++){
      int rand_t = (int )(rand() % numberOfTriangles);
      int rand_z = (int )(rand() % (numberOfLevels-1));
      numberOfImportantRays[rand_t + rand_z * numberOfTriangles]++;
  
    }

  // Now think about the mount of rays which would come out of this volume(surface)
  // dividing this number with the new amount of rays gives the final importance weight for this area!
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
  //return raysDump;
  return raysPerSample;
}
