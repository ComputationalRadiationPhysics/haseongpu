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

struct NormalRay {
  TwoDimPoint p;
  TwoDimDir dir;

  __host__ __device__ double length();
  __host__ __device__ void normalize();
};


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


struct Triangle;

/**
 * @brief One of the 3 sides of a triangle
 **/
struct Edge {
  NormalRay normal;
  Triangle *neighbor;
  int forbidden;
};


/**
 * @brief A triangle of the meshed structure
 *
 **/
struct Triangle {
  TwoDimPoint A;
  TwoDimPoint B;
  TwoDimPoint C;

  double *betaValues; // array of betaValues, one for each prism/level above the triangle
  Edge edges[3]; // edges of the triangle size: 3
  TwoDimPoint center;
  float surface;
};

/**
 * @brief Contains the structure of the crystal
 *
 * All the fixed values of how the crystal is meshed
 **/
struct Mesh {
  Triangle *triangles;
  Point *samples;
  unsigned numberOfTriangles;
  unsigned numberOfLevels;
  unsigned numberOfPrisms;
  unsigned numberOfPoints;
  unsigned numberOfSamples;
  float thickness;
  float surface;

  ~Mesh();

  static void parse(Mesh *hMesh, Mesh *dMesh, std::vector<unsigned> *triangleIndices, unsigned numberOfTriangles, unsigned numberOfLevels, unsigned numberOfPoints, float thicknessOfPrism, std::vector<double> *pointXY, std::vector<double> *betaValues, std::vector<double> *xOfTriangleCenter, std::vector<double> *yOfTriangleCenter, std::vector<int> *positionsOfNormalVectors, std::vector<double> *xOfNormals, std::vector<double> *yOfNormals, std::vector<int> *forbidden, std::vector<int> *neighbors, std::vector<float> *surfaces);

  static void parseMultiGPU(Mesh *hMesh, Mesh **dMesh, std::vector<unsigned> *triangleIndices, unsigned numberOfTriangles, unsigned numberOfLevels, unsigned numberOfPoints, float thicknessOfPrism, std::vector<double> *pointXY, std::vector<double> *betaValues, std::vector<double> *xOfTriangleCenter, std::vector<double> *yOfTriangleCenter, std::vector<int> *positionsOfNormalVectors, std::vector<double> *xOfNormals, std::vector<double> *yOfNormals, std::vector<int> *forbidden, std::vector<int> *neighbors, std::vector<float> *surfaces,unsigned numberOfDevices,unsigned *devices);
};

#endif /* MESH_H */
