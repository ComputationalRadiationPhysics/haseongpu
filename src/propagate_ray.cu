#include <mesh.h>
#include <geometry.h>
#include <cuda_runtime_api.h>

__host__ __device__ double checkSurface(int currentLevel, double zPos, double zVec, double length, const double thickness){
	double denominator = zPos * zVec;
	if (denominator != 0.0){
		double nominator = currentLevel * thickness - zPos;
		double lengthTmp = nominator/denominator;
		if (lengthTmp < length && lengthTmp > 0.0){
			return lengthTmp;
		}
	}
	return 0;
}

__host__ __device__ double checkEdge(Triangle triangle, int edge, Ray ray, double length){
  double denominator = triangle.edges[edge].normal.dir.x * ray.dir.x +  triangle.edges[edge].normal.dir.y * ray.dir.y;
  if (denominator != 0.0)
    {
      double nominator =	  
	triangle.edges[edge].normal.dir.x * triangle.edges[edge].normal.p.x
	+ triangle.edges[edge].normal.dir.y * triangle.edges[edge].normal.p.y
	- triangle.edges[edge].normal.dir.x * ray.p.x 
	- triangle.edges[edge].normal.dir.y * ray.p.y; 

      double lengthTmp = nominator/denominator;

      if(lengthTmp < length && lengthTmp > 0.0){
	return lengthTmp;
      }
    }
  return 0;
}

__host__ __device__ int getNextEdge(Triangle triangle,  Ray ray, unsigned level, double length, int forbiddenEdge, double thickness){
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

__host__ __device__ unsigned getNextLevel(unsigned level, int edge){
  switch(edge){
  case 3:
    return ++level;
  case 4:
    return --level;
  default:
    return level;
  }

}

__host__ __device__ double calcTriangleIntersection(Triangle triangle, Ray ray, int edge, int length, unsigned level, double thickness){
  switch(edge){
  case 0:
  case 1:
  case 2:
    return checkEdge(triangle, edge, ray, length);
  case 3:
  case 4:
    return checkSurface(level + 1, ray.p.z, ray.dir.z, length, thickness);

  }
  return 0;

}

__host__ __device__ Ray calcNextRay(Ray ray, double length){
  ray.p.x = ray.p.x + ray.length * ray.dir.x;
  ray.p.y = ray.p.y + ray.length * ray.dir.y;
  ray.p.z = ray.p.z + ray.length * ray.dir.z;
  ray.length = ray.length - length;

  return ray;

}

__host__ __device__ double calcPrismGain(const Triangle triangle, const unsigned level, double length, const double sigmaA, const double sigmaE, const double nTot){
  return (double) exp(nTot * (triangle.betaValues[level] * ( sigmaE + sigmaA ) - sigmaA ) * length);
 
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

  while(fabs(distanceRemaining) < SMALL){
    nextEdge          = getNextEdge(nextTriangle, nextRay, nextLevel, distanceRemaining, nextForbiddenEdge, thickness);
    nextLevel         = getNextLevel(nextLevel, nextEdge);
    length            = calcTriangleIntersection(nextTriangle, nextRay, nextEdge, distanceRemaining, nextLevel, thickness);
    nextRay           = calcNextRay(nextRay, length);
    nextForbiddenEdge = nextTriangle.edges[nextEdge].forbidden;
    nextTriangle      = *(nextTriangle.edges[nextEdge].neighbor);

    gain *= calcPrismGain(nextTriangle, nextLevel, length, sigmaA, sigmaE, nTot);
    distanceRemaining -= length;

  }

  return gain /= (distanceTotal * distanceTotal);

}
