#include "mesh.h"
#include <stdio.h>
#include <cudachecks.h>

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

  cudaMalloc((void**) &dMesh->samples, numberOfPoints*numberOfLevels*sizeof(Point));
  cudaMemcpy(dMesh->samples, hMesh->samples, numberOfPoints*numberOfLevels*sizeof(Point), cudaMemcpyHostToDevice);

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

void fillHMesh(
    Mesh *hMesh,
    std::vector<unsigned> *triangleIndices, 
    unsigned numberOfTriangles, 
    unsigned numberOfLevels,
    unsigned numberOfPoints, 
    float thicknessOfPrism,
    TwoDimPoint *points, 
    std::vector<double> *betaValues, 
    std::vector<double> *xOfTriangleCenter, 
    std::vector<double> *yOfTriangleCenter, 
    std::vector<int> *positionsOfNormalVectors,
    std::vector<double> *xOfNormals, 
    std::vector<double> *yOfNormals,
    std::vector<int> *forbidden, 
    std::vector<int> *neighbors, 
    std::vector<float> *surfaces
    ) {

  hMesh->numberOfTriangles = numberOfTriangles;
  hMesh->numberOfLevels = numberOfLevels;
  hMesh->numberOfPrisms = numberOfTriangles*(numberOfLevels-1);
  hMesh->numberOfPoints = numberOfPoints;
  hMesh->numberOfSamples = numberOfPoints * numberOfLevels;
  hMesh->thickness = thicknessOfPrism;

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

  hMesh->triangles = new Triangle[numberOfTriangles];
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
    }

    hMesh->triangles[i].betaValues = new double[numberOfLevels-1];
    for(unsigned l=0; l<(numberOfLevels-1); ++l) {
      hMesh->triangles[i].betaValues[l] = betaValues->at(l*numberOfTriangles + i);
    }
  }
  hMesh->surface = totalSurface;
}
void fillDMesh(
    Mesh *hMesh,
    Mesh *dMesh, 
    std::vector<unsigned> *triangleIndices, 
    unsigned numberOfTriangles, 
    unsigned numberOfLevels,
    unsigned numberOfPoints, 
    float thicknessOfPrism,
    TwoDimPoint *points, 
    std::vector<double> *xOfTriangleCenter, 
    std::vector<double> *yOfTriangleCenter, 
    std::vector<int> *positionsOfNormalVectors,
    std::vector<double> *xOfNormals, 
    std::vector<double> *yOfNormals,
    std::vector<int> *forbidden, 
    std::vector<int> *neighbors, 
    std::vector<float> *surfaces
    ) {


  dMesh->numberOfTriangles = numberOfTriangles;
  dMesh->numberOfLevels = numberOfLevels;
  dMesh->numberOfPrisms = numberOfTriangles*(numberOfLevels-1);
  dMesh->numberOfPoints = numberOfPoints;
  dMesh->numberOfSamples = numberOfPoints*numberOfLevels;
  dMesh->thickness = thicknessOfPrism;

  CUDA_CHECK_RETURN( cudaMalloc((void**) &dMesh->samples, numberOfPoints*numberOfLevels*sizeof(Point)) );
  CUDA_CHECK_RETURN( cudaMemcpy(dMesh->samples, hMesh->samples, numberOfPoints*numberOfLevels*sizeof(Point), cudaMemcpyHostToDevice) );

  Triangle *trianglesForDevice = new Triangle[numberOfTriangles];
  CUDA_CHECK_RETURN( cudaMalloc((void**) &dMesh->triangles, numberOfTriangles*sizeof(Triangle)) );

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

    trianglesForDevice[i] = triangle;

    for(unsigned e=0; e<3; ++e) {
      NormalRay normal;
      normal.p = points[positionsOfNormalVectors->at(e*numberOfTriangles + i)];
      normal.dir.x = xOfNormals->at( e*numberOfTriangles + i );
      normal.dir.y = yOfNormals->at( e*numberOfTriangles + i );

      Edge edge;
      edge.normal = normal;
      edge.forbidden = forbidden->at( e*numberOfTriangles + i);

      edge.neighbor = &(dMesh->triangles[neighbors->at( e*numberOfTriangles + i)]);
      trianglesForDevice[i].edges[e] = edge;
    }

    CUDA_CHECK_RETURN( cudaMalloc((void**) &trianglesForDevice[i].betaValues, (numberOfLevels-1)*sizeof(double)) );
    CUDA_CHECK_RETURN( cudaMemcpy(trianglesForDevice[i].betaValues, hMesh->triangles[i].betaValues, (numberOfLevels-1)*sizeof(double), cudaMemcpyHostToDevice) );
  }
  dMesh->surface = totalSurface;

  CUDA_CHECK_RETURN( cudaMemcpy(dMesh->triangles, trianglesForDevice, numberOfTriangles*sizeof(Triangle), cudaMemcpyHostToDevice) );

}

void Mesh::parseMultiGPU(
    Mesh *hMesh,
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
    unsigned *devices) {

  TwoDimPoint *points = parsePoints(pointXY, numberOfPoints);

  fillHMesh(
      hMesh,
      triangleIndices, 
      numberOfTriangles, 
      numberOfLevels,
      numberOfPoints, 
      thicknessOfPrism,
      points, 
      betaValues, 
      xOfTriangleCenter, 
      yOfTriangleCenter, 
      positionsOfNormalVectors,
      xOfNormals, 
      yOfNormals,
      forbidden, 
      neighbors, 
      surfaces
      );

 for( unsigned i=0;i<numberOfDevices;i++){
  CUDA_CHECK_RETURN( cudaSetDevice(devices[i]) );
  fillDMesh(
      hMesh,
      &((*dMesh)[i]),
      triangleIndices, 
      numberOfTriangles, 
      numberOfLevels,
      numberOfPoints, 
      thicknessOfPrism,
      points, 
      xOfTriangleCenter, 
      yOfTriangleCenter, 
      positionsOfNormalVectors,
      xOfNormals, 
      yOfNormals,
      forbidden, 
      neighbors, 
      surfaces
      );
  cudaDeviceSynchronize();
 }
}
