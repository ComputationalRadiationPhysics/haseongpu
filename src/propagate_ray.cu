#include <mesh.h>
#include <stdio.h>
#include <geometry.h>
#include <cuda_runtime_api.h>
#include <stdio.h> /* printf */

/**
 * @brief Checks a level-plane(currentLevel * thickness) for intersection with an ray (zPos, zVec).
 *        If the intersection-length is greater then length. 
 *        Than the intersection-length will be returned. 
 *        Otherwise 0 will be returned.
 *
 * @return intersection-length if intersection-length <= length
 * @return 0 if intersection-length > length
 *
 **/
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

/**
 * @brief Checks an edges of the given triangle/prism for an intersection
 *        with ray and calculates the intersection-length. If the intersection-length
 *        is greater then length. Than the intersection-length will be
 *        returned. Otherwise 0 will be returned.
 *
 * @return intersection-length if intersection-length <= length
 * @return 0 if intersection-length > length
 **/
__host__ __device__ double checkEdge(const unsigned triangle, const int edge, const Ray ray, Mesh *mesh, const double length){
  NormalRay normal = mesh->getNormal(triangle, edge);
  double denominator = normal.dir.x * ray.dir.x + normal.dir.y * ray.dir.y;

  if (denominator != 0.0)
    {
      double nominator =	  
	normal.dir.x * normal.p.x
	+ normal.dir.y * normal.p.y
	- normal.dir.x * ray.p.x 
	- normal.dir.y * ray.p.y; 

      double lengthTmp = nominator/denominator;
      if(lengthTmp <= length && lengthTmp > 0.0){
	return lengthTmp;
      }

    }
  
  return 0;
}

/**
 * @brief Calculates the intersection-length for the propagated ray and
 *        the current triangle.
 *
 * @return edge number of the intesected edge (-1 for no intersection)
 *
 **/
__host__ __device__ int calcTriangleRayIntersection(double *length, const unsigned triangle,  const Ray ray, const unsigned level, const int forbiddenEdge, Mesh *mesh, const double thickness){
  int edge = -1;
  // Check 3 edges of triangle
  for(int edge_i = 0; edge_i < 3; ++edge_i){
    if(edge_i != forbiddenEdge){
      double lengthTmp = checkEdge(triangle, edge_i, ray, mesh, *length);
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

/**
 * @brief This is simple vector calculation. The startpoint
 *        of ray will be moved by length.
 * 
 * @return ray is the ray with moved startpoint
 *
 **/
__host__ __device__ Ray calcNextRay(Ray ray, const double length){
  ray.p.x = ray.p.x + length * ray.dir.x;
  ray.p.y = ray.p.y + length * ray.dir.y;
  ray.p.z = ray.p.z + length * ray.dir.z;

  return ray;

}

/**
 * @brief Calculates the gain for the given prism(triangle and level) and 
 *        the intersection-length of the ray.
 *
 * @return gain
 *
 **/
__host__ __device__ double calcPrismGain(const unsigned triangle, const unsigned level, const double length, Mesh *mesh, const double sigmaA, const double sigmaE, const double nTot){
  return (double) exp(nTot * (mesh->getBetaValue(triangle, level) * ( sigmaE + sigmaA ) - sigmaA ) * length);
 
}

/**
 * @brief Sets the next triangle, next forbiddenEdge 
 *        and next level depending on the cutted edge of 
 *        the current triangle and the propagated ray.
 *
 **/
__host__ __device__ void updateFromEdge(unsigned *triangle, int *forbiddenEdge, unsigned *level, Mesh *mesh, const int edge){
   switch(edge){
   case 0:
   case 1:
   case 2:
     // One of three edges
     *forbiddenEdge = mesh->getForbiddenEdge(*triangle, edge);
     *triangle = mesh->getNeighbor(*triangle, edge);
     break;

   case 3:
     // Upper surface
     *forbiddenEdge = 4;
     (*level)++;
     break;

   case 4:
     // Lower surface
     *forbiddenEdge = 3;
     (*level)--;
     break;

  }

}

__host__ __device__ double propagateRay(Ray nextRay, unsigned nextLevel, unsigned nextTriangle, Mesh *mesh, 
					const double sigmaA, const double sigmaE, const double nTot, const double thickness){
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
    nextEdge           = calcTriangleRayIntersection(&length, nextTriangle, nextRay, nextLevel, nextForbiddenEdge, mesh, thickness);
    nextRay            = calcNextRay(nextRay, length);
    gain              *= calcPrismGain(nextTriangle, nextLevel, length, mesh, sigmaA, sigmaE, nTot);
    distanceRemaining -= length;

    // Calc nextTriangle, nextForbiddenEdge and nextLevel
    if(nextEdge != -1){
      updateFromEdge(&nextTriangle, &nextForbiddenEdge, &nextLevel, mesh, nextEdge);
    }

  }

  return gain /= (distanceTotal * distanceTotal);
}
