#include "mesh.h"
#include <stdio.h>
#include <cudachecks.h>

/**
 * @brief converts a vector of points into a vector of TwoDimPoint
 *
 * @param *points an array of points, containing numPoints x-values, followed by numPoints y-values
 *
 * @param numPoints the number of points which are stored
 *
 * @return an array of TwoDimPoint with the length numPoints 
 *
 */
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

/**
 * @brief creates the Mesh datastructures on device and host for the propagation
 *
 * @param *hMesh the host mesh
 *
 * @param *dMesh the mesh on the device
 *
 * @param *triangleIndices indices of the points which form a triangle
 *
 * @param numberOfTriangles the number of triangles
 *
 * @param numberOfLeves the number of layers of the mesh
 *
 * @param numberOfPoints the number of vertices in one layer of the mesh
 *
 * @param thicknessOfPrism  the thickness of one layer of the mesh
 *
 * @param *pointXY coordinates of the vertices in one layer of the mesh
 * 
 * @param *betaValues constant values for each meshed prism
 *
 * @param *xOfTriangleCenter the x coordinates of each triangle's center
 *
 * @param *yOfTriangleCenter the y coordinates of each triangle's center
 *
 * @param *positionsOfNormalVectors indices to the points (pointXY), where the normals xOfNormals,yOfNormals start
 *
 * @param *xOfNormals the x components of a normal vector for each of the 3 sides of a triangle
 *
 * @param *yOfNormals the y components of a normal vector for each of the 3 sides of a triangle
 *
 * @param *forbidden the sides of the triangle from which a ray "entered" the triangle
 *
 * @param *neighbors indices to the adjacent triangles in triangleIndices
 *
 * @param *surfaces the sizes of the surface of each prism
 *
 */
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

/**
 * @brief fills the host mesh with the correct datastructures
 *
 * See parseMultiGPU for details on the parameters
 */
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


/**
 * @brief fills a device mesh with the correct datastructures
 *
 * See parseMultiGPU for details on the parameters
 */
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

/**
 * @brief creates the Mesh datastructures on the host and on all possible devices for the propagation
 *
 * @param *hMesh the host mesh
 *
 * @param **dMesh an array of device meshes (one for each device) 
 *
 * @param *triangleIndices indices of the points which form a triangle
 *
 * @param numberOfTriangles the number of triangles
 *
 * @param numberOfLeves the number of layers of the mesh
 *
 * @param numberOfPoints the number of vertices in one layer of the mesh
 *
 * @param thicknessOfPrism  the thickness of one layer of the mesh
 *
 * @param *pointXY coordinates of the vertices in one layer of the mesh
 * 
 * @param *betaValues constant values for each meshed prism
 *
 * @param *xOfTriangleCenter the x coordinates of each triangle's center
 *
 * @param *yOfTriangleCenter the y coordinates of each triangle's center
 *
 * @param *positionsOfNormalVectors indices to the points (pointXY), where the normals xOfNormals,yOfNormals start
 *
 * @param *xOfNormals the x components of a normal vector for each of the 3 sides of a triangle
 *
 * @param *yOfNormals the y components of a normal vector for each of the 3 sides of a triangle
 *
 * @param *forbidden the sides of the triangle from which a ray "entered" the triangle
 *
 * @param *neighbors indices to the adjacent triangles in triangleIndices
 *
 * @param *surfaces the sizes of the surface of each prism
 *
 * @param numberOfDevices number of devices in *devices
 *
 * @param *devices array of device indices for all possible devices 
 *
 */
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
