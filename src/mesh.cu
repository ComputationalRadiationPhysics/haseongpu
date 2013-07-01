#include "mesh.h"

double NormalRay::length() {
  return sqrt(dir.x*dir.x + dir.y*dir.y);
}

void NormalRay::normalize() {
  double l = length();
  p.x = p.x/l;
  p.y = p.y/l;
}

std::vector<TwoDimPoint> parsePoints(std::vector<double> *points, int numPoints) {
  std::vector<TwoDimPoint> p;

  for(unsigned i=0; i < numPoints; ++i) {
    TwoDimPoint point;
    point.x = points->at(i);
    point.y = points->at(numPoints + i);
    p.push_back(point);
  }

  return p;
}

Mesh::~Mesh() {
  if(!triangles) delete triangles;
}

void Mesh::parse(Mesh *hMesh, Mesh *dMesh, std::vector<unsigned> *triangleIndices, unsigned numberOfTriangles, unsigned numberOfLevels, unsigned numberOfPoints, std::vector<double> *pointXY, std::vector<double> *betaValues, std::vector<double> *xOfTriangleCenter, std::vector<double> *yOfTriangleCenter, std::vector<int> *positionsOfNormalVectors, std::vector<double> *xOfNormals, std::vector<double> *yOfNormals, std::vector<int> *forbidden, std::vector<int> *neighbors, std::vector<float> *surfaces) {
  hMesh->numberOfTriangles = numberOfTriangles;
  dMesh->numberOfTriangles = numberOfTriangles;
  hMesh->numberOfLevels = numberOfLevels;
  dMesh->numberOfLevels = numberOfLevels;

  std::vector<TwoDimPoint> points = parsePoints(pointXY, numberOfPoints);

  hMesh->triangles = new Triangle[numberOfTriangles];
  Triangle *trianglesForDevice = new Triangle[numberOfTriangles];
  cudaMalloc((void**) &dMesh->triangles, numberOfTriangles*sizeof(Triangle));

  for(unsigned i=0; i<numberOfTriangles; ++i) {
    Triangle triangle;
    triangle.A = points.at( triangleIndices->at(i) );
    triangle.B = points.at( triangleIndices->at(numberOfTriangles + i) );
    triangle.C = points.at( triangleIndices->at(2*numberOfTriangles + i) );

    TwoDimPoint center = {xOfTriangleCenter->at(i), yOfTriangleCenter->at(i)};
    triangle.center = center;
    triangle.surface = surfaces->at(i);

    hMesh->triangles[i] = triangle;
    trianglesForDevice[i] = triangle;

    for(unsigned e=0; e<3; ++e) {
      NormalRay normal;
      normal.p = points.at( positionsOfNormalVectors->at(e*numberOfTriangles + i) );
      normal.dir.x = xOfNormals->at( e*numberOfTriangles + i );
      normal.dir.y = yOfNormals->at( e*numberOfTriangles + i );

      Edge edge;
      edge.normal = normal;
      edge.forbidden = forbidden->at( e*numberOfTriangles + i);

      edge.neighbor = &hMesh->triangles[neighbors->at( e*numberOfTriangles + i)];
      hMesh->triangles[i].edges[e] = edge;

      edge.neighbor = &dMesh->triangles[neighbors->at( e*numberOfTriangles + i)];
      trianglesForDevice[i].edges[e] = edge;
    }

    hMesh->triangles[i].betaValues = new double[numberOfLevels-1];
    for(unsigned l=0; l<(numberOfLevels-1); ++l) {
      hMesh->triangles[i].betaValues[l] = betaValues->at(l*numberOfTriangles + i);
    }
    cudaMalloc((void**) &trianglesForDevice[i].betaValues, (numberOfLevels-1)*sizeof(double));
    cudaMemcpy(trianglesForDevice[i].betaValues, hMesh->triangles[i].betaValues, (numberOfLevels-1)*sizeof(double), cudaMemcpyHostToDevice);
  }

  cudaMemcpy(dMesh->triangles, trianglesForDevice, numberOfTriangles*sizeof(Triangle), cudaMemcpyHostToDevice);
}
