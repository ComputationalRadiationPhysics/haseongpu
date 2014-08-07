/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASENonGPU
 *
 * HASENonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASENonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASENonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


#include <stdio.h>
#include <vector>
#include <string>
#include <assert.h>
#include <cfloat>
#include <cmath>
#include <algorithm>

#include <cudachecks.hpp>
#include <mesh.hpp>
#include <parser.hpp>
#include <reflection.hpp>


template <class T, class B, class E>
void assertRange(const std::vector<T> &v, const B minElement,const E maxElement, const bool equals){
  if(equals){
    assert(*std::min_element(v.begin(),v.end()) == minElement);
    assert(*std::max_element(v.begin(),v.end()) == maxElement);
  }else{
    assert(*std::min_element(v.begin(),v.end()) >= minElement);
    assert(*std::max_element(v.begin(),v.end()) <= maxElement);
  }
}

template <class T, class B>
void assertMin(const std::vector<T>& v,const  B minElement,const bool equals){
  if(equals){
    assert(*std::min_element(v.begin(),v.end()) == minElement);
  }else{
    assert(*std::min_element(v.begin(),v.end()) >= minElement);
  }
}

Mesh::Mesh(// Constants
	   double claddingAbsorption,
	   float surfaceTotal,
	   float thickness,
	   float nTot,
	   float crystalTFluo,
	   unsigned numberOfTriangles,
	   unsigned numberOfLevels,
	   unsigned numberOfPrisms,
	   unsigned numberOfPoints,
	   unsigned numberOfSamples,
	   unsigned claddingNumber,
	   // Vectors
	   std::vector<double> points,
	   std::vector<double> normalVec,
	   std::vector<double> betaVolume,
	   std::vector<double> centers,
	   std::vector<float> triangleSurfaces,
	   std::vector<int> forbiddenEdge,
	   std::vector<double> betaCells,
	   std::vector<unsigned> claddingCellTypes,
	   std::vector<float> refractiveIndices,
	   std::vector<float> reflectivities,
	   std::vector<float> totalReflectionAngles,
	   std::vector<unsigned> trianglePointIndices,
	   std::vector<int> triangleNeighbors,
	   std::vector<unsigned> triangleNormalPoint) :
  claddingAbsorption(claddingAbsorption),
  surfaceTotal(surfaceTotal),
  thickness(thickness),
  nTot(nTot),
  crystalTFluo(crystalTFluo),
  numberOfTriangles(numberOfTriangles),
  numberOfLevels(numberOfLevels),
  numberOfPrisms(numberOfPrisms),
  numberOfPoints(numberOfPoints),
  numberOfSamples(numberOfSamples),
  claddingNumber(claddingNumber),
  points(points),
  normalVec(normalVec),
  betaVolume(betaVolume),
  centers(centers),
  triangleSurfaces(triangleSurfaces),
  forbiddenEdge(forbiddenEdge),
  betaCells(betaCells),
  claddingCellTypes(claddingCellTypes),
  refractiveIndices(refractiveIndices),
  reflectivities(reflectivities),
  totalReflectionAngles(totalReflectionAngles),
  trianglePointIndices(trianglePointIndices),
  triangleNeighbors(triangleNeighbors),
  triangleNormalPoint(triangleNormalPoint){

}
Mesh::~Mesh() {
}

/**
 * @brief fetch the id of an adjacent triangle
 *
 * @param triangle the index of the triangle of which you want the neighbor
 * @param edge the side of the triangle, for whih you want the neighbor
 *
 * @return the index of the neighbor triangle
 */
__device__ int Mesh::getNeighbor(unsigned triangle, int edge) const{
  return triangleNeighbors[triangle + edge*numberOfTriangles];
}

/**
 * @brief generate a random point within a prism
 *
 * @param triangle the triangle to describe the desired prism
 * @param the level of the desired prism
 * @param *globalState a global state for the Mersenne Twister PRNG
 *
 * @return random 3D point inside the desired prism
 *
 * Uses a Mersenne Twister PRNG and Barycentric coordinates to generate a
 * random position inside a given triangle in a specific depth
 */
__device__ Point Mesh::genRndPoint(unsigned triangle, unsigned level, curandStateMtgp32 *globalState) const{
  Point startPoint = {0,0,0};
  double u = curand_uniform_double(&globalState[blockIdx.x]);
  double v = curand_uniform_double(&globalState[blockIdx.x]);

  if((u+v)>1)
  {
    u = 1-u;
    v = 1-v;
  }
  double w = 1-u-v;
  int t1 = trianglePointIndices[triangle];
  int t2 = trianglePointIndices[triangle + numberOfTriangles];
  int t3 = trianglePointIndices[triangle + 2 * numberOfTriangles];

  // convert the random startpoint into coordinates
  startPoint.z = (level + curand_uniform_double(&globalState[blockIdx.x])) * thickness;
  startPoint.x = (points[t1] * u) + (points[t2] * v) + (points[t3] * w);
  startPoint.y = (points[t1+numberOfPoints] * u) + (points[t2+numberOfPoints] * v) + (points[t3+numberOfPoints] * w);

  return startPoint;
}


/**
 * @brief get a betaVolume for a specific triangle and level
 *
 * @param triangle the id of the desired triangle
 * @param level the level of the desired prism
 *
 * @return a beta value
 */
__device__ double Mesh::getBetaVolume(unsigned triangle, unsigned level) const{
  return betaVolume[triangle + level * numberOfTriangles];
}

/**
 * @brief get a betaVolume for a specific prism
 *
 * @param the id of the desired prism
 *
 * @return a beta value
 */
__device__ double Mesh::getBetaVolume(unsigned prism) const{
  return betaVolume[prism];
}

/**
 * @brief generates a normal vector for a given side of a triangle
 *
 * @param triangle the desired triangle
 * @param the edge (0,1,2) of the triangle
 *
 * @return a normal vector with length 1
 */
__device__ NormalRay Mesh::getNormal(unsigned triangle, int edge) const{
  NormalRay ray = { {0,0},{0,0}};
  int offset =  edge*numberOfTriangles + triangle;
  ray.p.x = points[ triangleNormalPoint [offset] ];
  ray.p.y = points[ triangleNormalPoint [offset] + numberOfPoints ];

  ray.dir.x = normalVec[offset];
  ray.dir.y = normalVec[offset + 3*numberOfTriangles];

  return ray;
}	

/**
 * @brief genenerates a point with the coordinates of a given vertex
 *
 * @param sample the id of the desired samplepoint
 *
 * @return the Point with correct 3D coordinates
 */
__device__ Point Mesh::getSamplePoint(unsigned sample_i) const{
  Point p = {0,0,0};
  unsigned level = sample_i/numberOfPoints;
  p.z = level*thickness;
  unsigned pos = sample_i - (numberOfPoints * level);
  p.x = points[pos];
  p.y = points[pos + numberOfPoints];
  return p;
}

/**
 * @brief get a Point in the center of a prism
 *
 * @param triangle the id of the desired triangle
 * @param level the level of the desired prism
 *
 * @return a point with the coordinates (3D) of the prism center
 */
__device__ Point Mesh::getCenterPoint(unsigned triangle,unsigned level) const{
  Point p = {0,0,(level+0.5)*thickness};
  p.x = centers[triangle];
  p.y = centers[triangle + numberOfTriangles];
  return p;
}

/**
 * @brief gets the edge-id which will be forbidden
 *
 * @param trianle the index of the triangle you are currently in
 * @param edge the index of the edge through which you are leaving the triangle
 *
 * @return the id of the edge, which will be forbidden in the new triangle
 *
 * The forbidden edge corresponds to the edge through which you left the
 * previous triangle (has a different index in the new triangle)
 */
__device__ int Mesh::getForbiddenEdge(unsigned triangle,int edge) const{
  return forbiddenEdge[edge * numberOfTriangles + triangle];
}


__device__ unsigned Mesh::getCellType(unsigned triangle) const{
  return claddingCellTypes[triangle];
}


double distance2D(const TwoDimPoint p1, const TwoDimPoint p2) {
	return abs(sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)));
}

double getMaxDistance(std::vector<TwoDimPoint> points) {
	double maxDistance = -1;

	for(unsigned p1=0 ; p1 < points.size() ; ++p1)
		for(unsigned p2 = p1; p2 < points.size() ; ++p2)
			maxDistance = max(maxDistance,distance2D(points[p1],points[p2]));

	return maxDistance;
}

double calculateMaxDiameter(const double* points, const unsigned offset) {
	TwoDimPoint minX = {DBL_MAX,0};
	TwoDimPoint minY = {0,DBL_MAX};
	TwoDimPoint maxX = {DBL_MIN,0};
	TwoDimPoint maxY = {0,DBL_MIN};

	for(unsigned p=0; p<offset; ++p){
		TwoDimPoint np = {points[p],points[p+offset]};
		minX = (points[p] < minX.x) ? np : minX;
		maxX = (points[p] > maxX.x) ? np : maxX;
	}
	for(unsigned p=offset;p<2*offset;++p){
        TwoDimPoint np = {points[p-offset],points[p]};
		minY = points[p]<minY.y ? np : minY;
		maxY = points[p]>maxY.y ? np : maxY;
	}

	std::vector<TwoDimPoint> extrema;
	extrema.push_back(minX);
	extrema.push_back(minY);
	extrema.push_back(maxX);
	extrema.push_back(maxY);
	

	return getMaxDistance(extrema);
}

unsigned Mesh::getMaxReflections (ReflectionPlane reflectionPlane) const{
  double d = calculateMaxDiameter(points.toArray(),numberOfPoints);
  float alpha = getReflectionAngle(reflectionPlane) * M_PI / 180.;
  double h = numberOfLevels * thickness; 
  double z = d/tan(alpha);
  return ceil(z/h);
}

unsigned Mesh::getMaxReflections() const{
	unsigned top = getMaxReflections(TOP_REFLECTION);
	unsigned bottom = getMaxReflections(BOTTOM_REFLECTION);
	return max(top,bottom);
}

__device__ __host__ float Mesh::getReflectivity(ReflectionPlane reflectionPlane, unsigned triangle) const{
	switch(reflectionPlane){
		case BOTTOM_REFLECTION:
			return reflectivities[triangle];
		case TOP_REFLECTION:
			return reflectivities[triangle + numberOfTriangles];
	}
	return 0;
}

__device__ __host__ float Mesh::getReflectionAngle(ReflectionPlane reflectionPlane) const{
	switch(reflectionPlane){
		case BOTTOM_REFLECTION:
      return totalReflectionAngles[0];
			//return asin(refractiveIndices[1]/refractiveIndices[0]);
		case TOP_REFLECTION:
      return totalReflectionAngles[1];
			//return asin(refractiveIndices[3]/refractiveIndices[2]);
	}
	return  0;
}

__device__ __host__ void Mesh::test() const{
  printf("Constants:\n");
  printf("claddingAbsorption: %f\n", claddingAbsorption);
  printf("surfaceTotal: %f\n", surfaceTotal);
  printf("thickness: %f\n", thickness);
  printf("nTot: %f\n", nTot);
  printf("crystalTFluo: %f\n", crystalTFluo);
  printf("numberOfTriangles: %u\n", numberOfTriangles);
  printf("numberOfLevels: %u\n", numberOfLevels);
  printf("numberOfPrisms: %u\n", numberOfPrisms);
  printf("numberOfPoints: %u\n", numberOfPoints);
  printf("numberOfSamples: %u\n", numberOfSamples);
  printf("claddingNumber: %u\n", claddingNumber);
  
}
