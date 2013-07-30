#ifndef MESH_H
#define MESH_H 

#include<vector>
#include <curand_kernel.h> /* curand_uniform */

#define SMALL 1E-06
#define VERY_SMALL 0.0

struct TwoDimPoint {
  double x;
  double y;
};

typedef TwoDimPoint TwoDimDir;

struct Point {
  double x;
  double y;
  double z;

};

typedef Point Vector;

/**
 * @brief a Ray, defined by a startpoint, direction and length
 */
struct Ray {
  Point p;
  Vector dir;
  float length;
  
};

struct NormalRay {
  TwoDimPoint p;
  TwoDimDir dir;
};



/**
 * @brief Contains the structure of the crystal
 *
 * All the fixed values of how the crystal is meshed
 **/
struct Mesh {

  // values
  double *points;
  double *betaValues;
  double *normalVec;
  double *centers;
  float  *surfaces;
  int	 *forbidden;

  //indexstructs
  unsigned *triangles;
  int *neighbors;
  unsigned *normalPoint;

  //constants
  float surfaceTotal;
  float thickness;
  unsigned numberOfTriangles;
  unsigned numberOfLevels;
  unsigned numberOfPrisms;
  unsigned numberOfPoints;
  unsigned numberOfSamples;

  ~Mesh();

  __device__ int getNeighbor(unsigned triangle, int edge);
  __device__ Point genRndPoint(unsigned triangle, unsigned level, curandStateMtgp32 *globalState);
  __device__ double getBetaValue(unsigned triangle, unsigned level);
  __device__ double getBetaValue(unsigned prism);
  __device__ NormalRay getNormal(unsigned triangle, int edge);
  __device__ Point getSamplePoint(unsigned sample);
  __device__ Point getCenterPoint(unsigned triangle, unsigned level);
  __device__ int getForbiddenEdge(unsigned triangle, int edge);


  static void parseMultiGPU(Mesh *hMesh, 
			    Mesh **dMesh, 
			    std::vector<unsigned> *triangleIndices, 
			    unsigned numberOfTriangles, 
			    unsigned numberOfLevels, 
			    unsigned numberOfPoints, 
			    float thicknessOfPrism, 
			    std::vector<double> *pointXY, 
			    std::vector<double> *betaValues, 
			    std::vector<double> *xOfTriangleCenter, 
			    std::vector<double> *yOfTriangleCenter, 
			    std::vector<int> *positionsOfNormalVectors, 
			    std::vector<double> *xOfNormals, 
			    std::vector<double> *yOfNormals, 
			    std::vector<int> *forbidden, 
			    std::vector<int> *neighbors, 
			    std::vector<float> *surfaces,
			    unsigned numberOfDevices,
			    unsigned *devices);
};

#endif /* MESH_H */
