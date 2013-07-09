#include "mesh.h"

double NormalRay::length() {
  return sqrt(dir.x*dir.x + dir.y*dir.y);
}

void NormalRay::normalize() {
  double l = length();
  p.x = p.x/l;
  p.y = p.y/l;
}

TwoDimPoint* parsePoints(std::vector<double> *points, unsigned numPoints) {
  TwoDimPoint *p = new TwoDimPoint[numPoints];

  for(unsigned i=0; i < numPoints; ++i) {
    p[i].x = points->at(i);
    p[i].y = points->at(numPoints + i);
  }

  return p;
}

Mesh::~Mesh() {
  if(!triangles) delete triangles;
}

void Mesh::parse(Mesh *hMesh, Mesh *dMesh, std::vector<unsigned> *triangleIndices, unsigned numberOfTriangles, unsigned numberOfLevels, unsigned numberOfPoints, float thicknessOfPrism, std::vector<double> *pointXY, std::vector<double> *betaValues, std::vector<double> *xOfTriangleCenter, std::vector<double> *yOfTriangleCenter, std::vector<int> *positionsOfNormalVectors, std::vector<double> *xOfNormals, std::vector<double> *yOfNormals, std::vector<int> *forbidden, std::vector<int> *neighbors, std::vector<float> *surfaces) {
  hMesh->numberOfTriangles = numberOfTriangles;
  dMesh->numberOfTriangles = numberOfTriangles;
  hMesh->numberOfLevels = numberOfLevels;
  dMesh->numberOfLevels = numberOfLevels;
  hMesh->numberOfPrisms = numberOfTriangles*(numberOfLevels-1);
  dMesh->numberOfPrisms = numberOfTriangles*(numberOfLevels-1);
  hMesh->numberOfPoints = numberOfPoints;
  dMesh->numberOfPoints = numberOfPoints;
  hMesh->numberOfSamples = numberOfPoints * numberOfLevels;
  dMesh->numberOfSamples = numberOfPoints*numberOfLevels;
  hMesh->thickness = thicknessOfPrism;
  dMesh->thickness = thicknessOfPrism;

  TwoDimPoint *points = parsePoints(pointXY, numberOfPoints);

  hMesh->samples = new Point[numberOfPoints * numberOfLevels];
  for(unsigned l=0; l<numberOfLevels; ++l) {
    for(unsigned i=0; i<numberOfPoints; ++i) {
      Point p;
      p.x = points[i].x;
      p.y = points[i].y;
      p.z = l*thicknessOfPrism;
      hMesh->samples[l*numberOfPoints + i] = p;
    }
  }

  cudaMalloc((void**) &dMesh->samples, numberOfPoints*sizeof(points));
  cudaMemcpy(dMesh->samples, hMesh->samples, numberOfPoints*sizeof(points), cudaMemcpyHostToDevice);

  hMesh->triangles = new Triangle[numberOfTriangles];
  Triangle *trianglesForDevice = new Triangle[numberOfTriangles];
  cudaMalloc((void**) &dMesh->triangles, numberOfTriangles*sizeof(Triangle));

  double totalSurface = 0;
  for(unsigned i=0; i<numberOfTriangles; ++i) {
    Triangle triangle;
    triangle.A = points[triangleIndices->at(i)];
    triangle.B = points[triangleIndices->at(numberOfTriangles + i)];
    triangle.C = points[triangleIndices->at(2*numberOfTriangles + i)];
    
    TwoDimPoint center = {xOfTriangleCenter->at(i), yOfTriangleCenter->at(i)};
    triangle.center = center;
    triangle.surface = surfaces->at(i);
    totalSurface += triangle.surface;

    hMesh->triangles[i] = triangle;
    trianglesForDevice[i] = triangle;

    for(unsigned e=0; e<3; ++e) {
      NormalRay normal;
      normal.p = points[positionsOfNormalVectors->at(e*numberOfTriangles + i)];
      normal.dir.x = xOfNormals->at( e*numberOfTriangles + i );
      normal.dir.y = yOfNormals->at( e*numberOfTriangles + i );

      Edge edge;
      edge.normal = normal;
      edge.forbidden = forbidden->at( e*numberOfTriangles + i);

      edge.neighbor = &(hMesh->triangles[neighbors->at( e*numberOfTriangles + i)]);
      hMesh->triangles[i].edges[e] = edge;

      edge.neighbor = &(dMesh->triangles[neighbors->at( e*numberOfTriangles + i)]);
      trianglesForDevice[i].edges[e] = edge;
    }

    hMesh->triangles[i].betaValues = new double[numberOfLevels-1];
    for(unsigned l=0; l<(numberOfLevels-1); ++l) {
      hMesh->triangles[i].betaValues[l] = betaValues->at(l*numberOfTriangles + i);
    }
    cudaMalloc((void**) &trianglesForDevice[i].betaValues, (numberOfLevels-1)*sizeof(double));
    cudaMemcpy(trianglesForDevice[i].betaValues, hMesh->triangles[i].betaValues, (numberOfLevels-1)*sizeof(double), cudaMemcpyHostToDevice);
  }
  hMesh->surface = totalSurface;
  dMesh->surface = totalSurface;

  cudaMemcpy(dMesh->triangles, trianglesForDevice, numberOfTriangles*sizeof(Triangle), cudaMemcpyHostToDevice);
}
