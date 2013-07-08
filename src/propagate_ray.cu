#include <mesh.h>
#include <assert.h>
#include <stdio.h>
#include <geometry.h>
#include <cuda_runtime_api.h>

__host__ __device__ double checkSurface(const int currentLevel, const double zPos, const double zVec, const float length, const double thickness){
  double denominator = zVec;
  if (denominator != 0.0){
    double nominator = currentLevel * thickness - zPos;
    double lengthTmp = nominator/denominator;
    // DEBUG
    //printf("C lengthTmp : %f nominator: %f denominaor %f level %d\n", lengthTmp, nominator, denominator, currentLevel);
    if (lengthTmp < length && lengthTmp > SMALL){
      return lengthTmp;
    }
    if(fabs(length - lengthTmp) < SMALL){
      return lengthTmp;
    }

  }
  else{
    printf("denominaor %f zPos %f zVec %f\n", denominator, zPos, zVec);
  }
  return 0;
}

__host__ __device__ double checkEdge(const Triangle triangle, const int edge, const Ray ray, const float length){
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
      if(lengthTmp < length && lengthTmp > SMALL){
	return lengthTmp;
      }
      // Should be in 
      if(fabs(length - lengthTmp) < SMALL){
      	return lengthTmp;
      }

    }
  else{
    printf("denominaor %f : %f %f %f %f\n", denominator, triangle.edges[edge].normal.dir.x, ray.dir.x, triangle.edges[edge].normal.dir.y, ray.dir.y);
  }

  return 0;
}

__host__ __device__ int getNextEdge(const Triangle triangle,  const Ray ray, const unsigned level, float length, const int forbiddenEdge, const double thickness){
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

__host__ __device__ unsigned getNextLevel(unsigned level, const int edge){
  switch(edge){
  case 3:
    return ++level;
  case 4:
    return --level;
  default:
    return level;
  }

}

__host__ __device__ double calcTriangleIntersection(const Triangle triangle, const Ray ray, const int edge, const float length, const unsigned level, const double thickness){
  switch(edge){
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

__host__ __device__ Ray calcNextRay(Ray ray, const float length){
  ray.p.x = ray.p.x + length * ray.dir.x;
  ray.p.y = ray.p.y + length * ray.dir.y;
  ray.p.z = ray.p.z + length * ray.dir.z;
  ray.length = ray.length - length;

  return ray;

}

__host__ __device__ double calcPrismGain(const Triangle triangle, const unsigned level, const float length, const double sigmaA, const double sigmaE, const double nTot){
  //return (double) exp(nTot * (triangle.betaValues[level] * ( sigmaE + sigmaA ) - sigmaA ) * length); // SEG FAULT
  return (double) exp(nTot * (triangle.betaValues[0] * ( sigmaE + sigmaA ) - sigmaA ) * length); // SEG FAULT
  //return (double) exp(nTot * (0.1 * ( sigmaE + sigmaA ) - sigmaA ) * length);
 
}

__host__ __device__ Triangle getNextTriangle(Triangle triangle, int edge){
  assert(edge >= 0 && edge <=4);
  if(edge == 3 || edge == 4)
    return triangle;
  return *(triangle.edges[edge].neighbor);

}

__host__ __device__ double propagateRay(Ray ray, unsigned startLevel, Triangle startTriangle, Triangle *triangles, const double sigmaA, const double sigmaE, const double nTot, const double thickness){
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

    // NOTICE was failure in surface intersection
    // SEG FAULT in calcPrismGain triangle.BetaValues[i];
    // DEBUG
    assert(nextForbiddenEdge >= -1 && nextForbiddenEdge <= 2);
    debugLoopCount++;
    assert(debugLoopCount <= 1000);
    // printf("\nC Calc next triangle, distanceRemaining: %f nextForbiddenEdge: %d nextLevel: %d\n", distanceRemaining, nextForbiddenEdge, nextLevel);
    // printf("C nextRay : POS(%f,%f,%f) VEC(%f, %f, %f)\n", nextRay.p.x, nextRay.p.y, nextRay.p.z, nextRay.dir.x, nextRay.dir.y, nextRay.dir.z);
    // printf("C nextTriangle: A(%f,%f), B(%f,%f), C(%f,%f)\n", 
    // 	   nextTriangle.A.x, nextTriangle.A.y,
    // 	   nextTriangle.B.x, nextTriangle.B.y,
    // 	   nextTriangle.C.x, nextTriangle.C.y);
    nextEdge          = getNextEdge(nextTriangle, nextRay, nextLevel, distanceRemaining, nextForbiddenEdge, thickness);
    length            = calcTriangleIntersection(nextTriangle, nextRay, nextEdge, distanceRemaining, nextLevel, thickness);
    nextLevel         = getNextLevel(nextLevel, nextEdge);
    nextRay           = calcNextRay(nextRay, length);
    nextForbiddenEdge = nextTriangle.edges[nextEdge].forbidden;
    nextTriangle      = getNextTriangle(nextTriangle, nextEdge);
    gain *= calcPrismGain(nextTriangle, nextLevel, length, sigmaA, sigmaE, nTot);

    distanceRemaining -= length;


    


  }

  return gain /= (distanceTotal * distanceTotal);

}
