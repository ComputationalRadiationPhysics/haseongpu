#include <mesh.h>
#include <stdio.h>
#include <geometry.h>
#include <cuda_runtime_api.h>
#include <stdio.h> /* printf */

__host__ __device__ double checkSurface(const int currentLevel, const double zPos, const double zVec, const double length, const double thickness){
  double denominator = zVec;
  if (denominator != 0.0){
    double nominator = currentLevel * thickness - zPos;
    double lengthTmp = nominator/denominator;
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
      if(lengthTmp <= length && lengthTmp > 0.0){
	return lengthTmp;
      }

    }
  
  return 0;
}

__host__ __device__ int calcTriangleRayIntersection(double *length, const Triangle triangle,  const Ray ray, const unsigned level, const int forbiddenEdge, const double thickness){
  int edge = -1;
  // Check 3 edges of triangle
  for(int edge_i = 0; edge_i < 3; ++edge_i){
    if(edge_i != forbiddenEdge){
      double lengthTmp = checkEdge(triangle, edge_i, ray, *length);
      if(lengthTmp){
	*length = lengthTmp;
	edge = edge_i;
      }
    }
  }
  
  // check the upper surface
  if (forbiddenEdge != 3){
    double lengthTmp = checkSurface(level + 1, ray.p.z, ray.dir.z, *length, thickness);
    if(lengthTmp){
      *length = lengthTmp;
      edge = 3;
    }
  }

  // check the lower surface
  if (forbiddenEdge != 4){
    double lengthTmp = checkSurface(level, ray.p.z, ray.dir.z, *length, thickness);
    if (lengthTmp){
      *length = lengthTmp;
      edge = 4;
    }
  }
  return edge;
}


__host__ __device__ Ray calcNextRay(Ray ray, const double length){
  ray.p.x = ray.p.x + length * ray.dir.x;
  ray.p.y = ray.p.y + length * ray.dir.y;
  ray.p.z = ray.p.z + length * ray.dir.z;

  return ray;

}

__host__ __device__ double calcPrismGain(const Triangle triangle, const unsigned level, const double length, const double sigmaA, const double sigmaE, const double nTot){
  return (double) exp(nTot * (triangle.betaValues[level] * ( sigmaE + sigmaA ) - sigmaA ) * length);
 
}

__host__ __device__ void updateFromEdge(Triangle *triangle, int *forbiddenEdge, unsigned *level, const int edge){
   switch(edge){
   case 0:
   case 1:
   case 2:
     // One of three edges
     *forbiddenEdge = triangle->edges[edge].forbidden;
     *triangle = *(triangle->edges[edge].neighbor);
     break;

   case 3:
     // Upper surface
     *forbiddenEdge = 4;
     (*level)++;
     break;

   case 4:
     // Lower surface
     *forbiddenEdge = 3;
     (*level--);
     break;

  }

}

__host__ __device__ double propagateRay(Ray nextRay, unsigned nextLevel, Triangle nextTriangle, const double sigmaA, const double sigmaE, const double nTot, const double thickness){
  double distanceTotal     = nextRay.length;
  double distanceRemaining = nextRay.length;
  double length  = 0;
  double gain    = 1;
  int nextForbiddenEdge = -1;
  int nextEdge          = -1;

  nextRay = normalizeRay(nextRay);
  while(fabs(distanceRemaining) > SMALL){
    // Calc gain for triangle intersection
    length             = distanceRemaining;
    nextEdge           = calcTriangleRayIntersection(&length, nextTriangle, nextRay, nextLevel, nextForbiddenEdge, thickness);
    nextRay            = calcNextRay(nextRay, length);
    gain              *= calcPrismGain(nextTriangle, nextLevel, length, sigmaA, sigmaE, nTot);
    distanceRemaining -= length;

    // Calc nextTriangle, nextForbiddenEdge and nextLevel
    if(nextEdge != -1){
      updateFromEdge(&nextTriangle, &nextForbiddenEdge, &nextLevel, nextEdge);
    }

  }

  return gain /= (distanceTotal * distanceTotal);
}
