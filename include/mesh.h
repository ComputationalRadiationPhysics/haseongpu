#ifndef MESH_H
#define MESH_H 

#include<vector>

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
}



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

  int getNeighbor(unsigned triangle, int edge);
  Point genRndPoint(unsigned triangle, unsigned level, curandStateMtgp32 *globalState);
  double getBetaValue(unsigned triangle, unsigned level);
  double getBetaValue(unsigned prism);
  NormalRay getNormal(unsigned triangle, int edge);
  Point getSamplePoint(unsigned sample);


  static void parse(Mesh *hMesh, Mesh *dMesh, std::vector<unsigned> *triangleIndices, unsigned numberOfTriangles, unsigned numberOfLevels, unsigned numberOfPoints, float thicknessOfPrism, std::vector<double> *pointXY, std::vector<double> *betaValues, std::vector<double> *xOfTriangleCenter, std::vector<double> *yOfTriangleCenter, std::vector<int> *positionsOfNormalVectors, std::vector<double> *xOfNormals, std::vector<double> *yOfNormals, std::vector<int> *forbidden, std::vector<int> *neighbors, std::vector<float> *surfaces);

  static void parseMultiGPU(Mesh *hMesh, Mesh **dMesh, std::vector<unsigned> *triangleIndices, unsigned numberOfTriangles, unsigned numberOfLevels, unsigned numberOfPoints, float thicknessOfPrism, std::vector<double> *pointXY, std::vector<double> *betaValues, std::vector<double> *xOfTriangleCenter, std::vector<double> *yOfTriangleCenter, std::vector<int> *positionsOfNormalVectors, std::vector<double> *xOfNormals, std::vector<double> *yOfNormals, std::vector<int> *forbidden, std::vector<int> *neighbors, std::vector<float> *surfaces,unsigned numberOfDevices,unsigned *devices);
};

#endif /* MESH_H */
