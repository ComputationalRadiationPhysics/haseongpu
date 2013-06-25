#ifndef MESH_H
#define MESH_H 

#include<vector>

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

struct Ray {
  Point p;
  Vector dir;
  unsigned length;
  
};


struct Triangle;

struct Edge {
  NormalRay normal;
  Triangle *neighbor; // index in global triangle list
  int forbidden;
};

struct Triangle {
  TwoDimPoint A;
  TwoDimPoint B;
  TwoDimPoint C;

  double *betaValues; // array of betaValues, one for each prism/level above the triangle
  Edge edges[3]; // edges of the triangle size: 3
  TwoDimPoint center;
  float surface;
};

class Mesh {
  Triangle *triangles;
  unsigned numberOfTriangles;

  public:

  ~Mesh();

  static Mesh parse(std::vector<unsigned> *triangleIndices, unsigned numberOfTriangles, unsigned numberOfLevels, unsigned numberOfPoints, std::vector<double> *pointXY, std::vector<double> *betaValues, std::vector<double> *xOfTriangleCenter, std::vector<double> *yOfTriangleCenter, std::vector<int> *positionsOfNormalVectors, std::vector<double> *xOfNormals, std::vector<double> *yOfNormals, std::vector<int> *forbidden, std::vector<int> *neighbors, std::vector<float> *surfaces);

  void copyToGpu();
  void copyFromGpu();
};

#endif /* MESH_H */
