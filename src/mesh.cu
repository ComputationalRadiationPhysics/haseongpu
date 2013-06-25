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
}

Mesh Mesh::parse(std::vector<unsigned> *triangleIndices, unsigned numberOfTriangles, unsigned numberOfLevels, unsigned numberOfPoints, std::vector<double> *pointXY, std::vector<double> *betaValues, std::vector<double> *xOfTriangleCenter, std::vector<double> *yOfTriangleCenter, std::vector<int> *positionsOfNormalVectors, std::vector<double> *xOfNormals, std::vector<double> *yOfNormals, std::vector<int> *forbidden, std::vector<int> *neighbors, std::vector<float> *surfaces) {
  Mesh mesh;

  std::vector<TwoDimPoint> points = parsePoints(pointXY, numberOfPoints);

  std::vector<Triangle> triangles;
  for(unsigned i=0; i<numberOfTriangles; ++i) {
    Triangle triangle;

    triangle.A = points.at( triangleIndices->at(i) );
    triangle.B = points.at( triangleIndices->at(numberOfTriangles + i) );
    triangle.C = points.at( triangleIndices->at(2*numberOfTriangles + i) );

    TwoDimPoint center;
    center.x = xOfTriangleCenter->at(i);
    center.y = yOfTriangleCenter->at(i);
    triangle.center = center;
    triangle.surface = surfaces->at(i);

    triangles.push_back(triangle);
  }

  // iterate again, because now all triangles are there (for the neighbor-pointer)
  for(unsigned i=0; i<triangles.size(); ++i) {
    Triangle *triangle = &triangles[i];
    for(unsigned e=0; e<3; ++e) {
      Edge edge;

      NormalRay normal;
      normal.p = points.at( positionsOfNormalVectors->at(e*numberOfTriangles + i) );
      normal.dir.x = xOfNormals->at( e*numberOfTriangles + i );
      normal.dir.y = yOfNormals->at( e*numberOfTriangles + i );
      edge.normal = normal;
      edge.forbidden = forbidden->at( e*numberOfTriangles + i);
      edge.neighbor = &(triangles[ neighbors->at( e*triangles.size() + i)]);

      triangle->edges[i] = edge;
    }
  }

  return mesh;
}

void Mesh::copyToGpu() {

}

void Mesh::copyFromGpu() {

}
