#include <stdio.h>
#include <vector>
#include <string>
#include <assert.h>
#include <cfloat>
#include <cmath>

#include <cudachecks.h>
#include <mesh.h>
#include <parser.h>


Mesh::~Mesh() {
//  if(points) delete points;
//  if(betaValues) delete betaValues;
//  if(normalVec) delete normalVec;
//  if(centers) delete centers;
//  if(surfaces) delete surfaces;
//  if(forbidden) delete forbidden;
//  if(betaCells) delete betaCells;
//  if(triangles) delete triangles;  
//  if(neighbors) delete neighbors;
//  if(normalPoint) delete normalPoint;
//  if(cellTypes) delete cellTypes;
//  if(refractiveIndices) delete refractiveIndices;
//  if(reflectivities) delete reflectivities;
}







/**
 * @brief fills the host mesh with the correct datastructures
 *
 * See parseMultiGPU for details on the parameters
 */
void fillHMesh(
    Mesh& hMesh,
    unsigned numberOfTriangles, 
    unsigned numberOfLevels,
    unsigned numberOfPoints, 
    float thicknessOfPrism,
    std::vector<unsigned> *triangleIndices,
    std::vector<double> *points,
    std::vector<double> *xOfTriangleCenter, 
    std::vector<double> *yOfTriangleCenter, 
    std::vector<unsigned> *positionsOfNormal,
    std::vector<double> *xOfNormals, 
    std::vector<double> *yOfNormals,
    std::vector<int> *forbidden, 
    std::vector<int> *neighbors, 
    std::vector<float> *surfaces,
    std::vector<double> *betaValues,
    std::vector<double> *betaCells,
    std::vector<unsigned> * cellTypes,
    std::vector<float> * refractiveIndices,
    std::vector<float> * reflectivities,
    float nTot,
    float crystalFluorescence,
    unsigned cladNumber,
    double cladAbsorption
) {

  hMesh.numberOfTriangles = numberOfTriangles;
  hMesh.numberOfLevels = numberOfLevels;
  hMesh.numberOfPrisms = numberOfTriangles*(numberOfLevels-1);
  hMesh.numberOfPoints = numberOfPoints;
  hMesh.numberOfSamples = numberOfPoints * numberOfLevels;
  hMesh.thickness = thicknessOfPrism;
  hMesh.crystalFluorescence = crystalFluorescence;
  hMesh.nTot = nTot;
  hMesh.cladNumber = cladNumber;
  hMesh.cladAbsorption = cladAbsorption;

  std::vector<double> *hostCenters = new std::vector<double>(xOfTriangleCenter->begin(), xOfTriangleCenter->end());
  hostCenters->insert(hostCenters->end(),yOfTriangleCenter->begin(),yOfTriangleCenter->end());

  std::vector<double> *hostNormalVec = new std::vector<double>(xOfNormals->begin(), xOfNormals->end());
  hostNormalVec->insert(hostNormalVec->end(),yOfNormals->begin(),yOfNormals->end());

  hMesh.points = &(points->at(0));
  hMesh.triangles = &(triangleIndices->at(0));
  hMesh.betaValues = &(betaValues->at(0));
  hMesh.normalVec = &(hostNormalVec->at(0));
  hMesh.centers = &(hostCenters->at(0));
  hMesh.surfaces = &(surfaces->at(0));
  hMesh.forbidden = &(forbidden->at(0));
  hMesh.neighbors = &(neighbors->at(0));
  hMesh.normalPoint = &(positionsOfNormal->at(0));
  hMesh.betaCells = &(betaCells->at(0));
  hMesh.cellTypes = &(cellTypes->at(0));
  hMesh.refractiveIndices = &(refractiveIndices->at(0));
  hMesh.reflectivities = &(reflectivities->at(0));

  std::vector<float> *totalReflectionAngles = new std::vector<float>(refractiveIndices->size()/2,0);
  for(unsigned i=0;i<refractiveIndices->size();i+=2){
    totalReflectionAngles->at(i/2) = (180. / M_PI * asin(refractiveIndices->at(i+1) / refractiveIndices->at(i)));
  }
  hMesh.totalReflectionAngles = &(totalReflectionAngles->at(0));

}

/**
 * @brief fills a device mesh with the correct datastructures
 *
 * See parseMultiGPU for details on the parameters
 */
void fillDMesh(
    Mesh& hMesh,
    Mesh& dMesh, 
    std::vector<unsigned> *triangleIndices, 
    unsigned numberOfTriangles, 
    unsigned numberOfLevels,
    unsigned numberOfPoints, 
    float thicknessOfPrism,
    std::vector<double> *pointsVector, 
    std::vector<double> *xOfTriangleCenter, 
    std::vector<double> *yOfTriangleCenter, 
    std::vector<unsigned> *positionsOfNormalVectors,
    std::vector<double> *xOfNormals, 
    std::vector<double> *yOfNormals,
    std::vector<int> *forbiddenVector, 
    std::vector<int> *neighborsVector, 
    std::vector<float> *surfacesVector,
    std::vector<double> *betaValuesVector,
    std::vector<double> *betaCells,
    std::vector<unsigned> *cellTypes,
    std::vector<float> * refractiveIndices,
    std::vector<float> * reflectivities,
    float nTot,
    float crystalFluorescence,
    unsigned cladNumber,
    double cladAbsorption
) {


  // GPU variables
  double totalSurface = 0.;

  // constants
  dMesh.numberOfTriangles = numberOfTriangles;
  dMesh.numberOfLevels = numberOfLevels;
  dMesh.numberOfPrisms = numberOfTriangles*(numberOfLevels-1);
  dMesh.numberOfPoints = numberOfPoints;
  dMesh.numberOfSamples = numberOfPoints*numberOfLevels;
  dMesh.thickness = thicknessOfPrism;
  dMesh.crystalFluorescence = crystalFluorescence;
  dMesh.nTot = nTot;
  dMesh.cladAbsorption = cladAbsorption;
  dMesh.cladNumber = cladNumber;

  for(unsigned i=0;i<numberOfTriangles;++i){
    totalSurface+=double(surfacesVector->at(i));	
  }
  dMesh.surfaceTotal = float(totalSurface);


  // values
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.points), 2 * hMesh.numberOfPoints * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.normalVec), 2 * 3 * hMesh.numberOfTriangles * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.betaValues), hMesh.numberOfPrisms * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.centers), 2 * hMesh.numberOfTriangles * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.surfaces), hMesh.numberOfTriangles * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.forbidden), 3 * hMesh.numberOfTriangles * sizeof(int)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.betaCells), hMesh.numberOfPoints * hMesh.numberOfLevels * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.cellTypes), hMesh.numberOfTriangles * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.refractiveIndices), refractiveIndices->size() * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.reflectivities), (refractiveIndices->size()/2)*(hMesh.numberOfTriangles) *sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.totalReflectionAngles), (refractiveIndices->size()/2) *sizeof(float)));

  // indexStructs
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.triangles), 3 * hMesh.numberOfTriangles * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.neighbors), 3 * hMesh.numberOfTriangles * sizeof(int)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.normalPoint), 3 * hMesh.numberOfTriangles * sizeof(unsigned)));


  /// fill values
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.points, (double*) &(pointsVector->at(0)), 2 * hMesh.numberOfPoints * sizeof(double), cudaMemcpyHostToDevice));

  std::vector<double> *hostNormalVec = new std::vector<double>(xOfNormals->begin(), xOfNormals->end());
  hostNormalVec->insert(hostNormalVec->end(),yOfNormals->begin(),yOfNormals->end());
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.normalVec, (double*) &(hostNormalVec->at(0)), 2 * 3 * hMesh.numberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
  free(hostNormalVec);

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.betaValues, (double*) &(betaValuesVector->at(0)), hMesh.numberOfPrisms * sizeof(double), cudaMemcpyHostToDevice));

  std::vector<double> *hostCenters = new std::vector<double>(xOfTriangleCenter->begin(), xOfTriangleCenter->end());
  hostCenters->insert(hostCenters->end(),yOfTriangleCenter->begin(),yOfTriangleCenter->end());
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.centers, (double*) &(hostCenters->at(0)), 2 * hMesh.numberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
  free(hostCenters);

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.surfaces, (float*) &(surfacesVector->at(0)), hMesh.numberOfTriangles * sizeof(float), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.forbidden, (int*) &(forbiddenVector->at(0)), 3 * hMesh.numberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.betaCells, (double*) &(betaCells->at(0)), hMesh.numberOfPoints * hMesh.numberOfLevels * sizeof(double), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.cellTypes, (unsigned*) &(cellTypes->at(0)), hMesh.numberOfTriangles * sizeof(unsigned), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.refractiveIndices, (float*) &(refractiveIndices->at(0)), refractiveIndices->size() * sizeof(float), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.reflectivities, (float*) &(reflectivities->at(0)), (hMesh.numberOfTriangles)* (refractiveIndices->size()/2) * sizeof(float), cudaMemcpyHostToDevice));

  std::vector<float> *totalReflectionAngles = new std::vector<float>(refractiveIndices->size()/2,0);
  for(unsigned i=0;i<refractiveIndices->size();i+=2){
    totalReflectionAngles->at(i/2) = (180. / M_PI *  asin(refractiveIndices->at(i+1) / refractiveIndices->at(i)));
  }
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.totalReflectionAngles, (float*) &(totalReflectionAngles->at(0)), refractiveIndices->size()/2 * sizeof(float), cudaMemcpyHostToDevice));
  free(totalReflectionAngles);

  // fill indexStructs
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.triangles, (unsigned*) &(triangleIndices->at(0)), 3 * hMesh.numberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.neighbors,(int*) &(neighborsVector->at(0)), 3 * hMesh.numberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.normalPoint, (unsigned*) &(positionsOfNormalVectors->at(0)), 3 * hMesh.numberOfTriangles * sizeof(unsigned), cudaMemcpyHostToDevice));

}

/**
 * @brief fetch the id of an adjacent triangle
 *
 * @param triangle the index of the triangle of which you want the neighbor
 *
 * @param edge the side of the triangle, for whih you want the neighbor
 *
 * @return the index of the neighbor triangle
 */
__device__ int Mesh::getNeighbor(unsigned triangle, int edge){
  return neighbors[triangle + edge*numberOfTriangles];
}

/**
 * @brief generate a random point within a prism
 *
 * @param triangle the triangle to describe the desired prism
 *
 * @param the level of the desired prism
 *
 * @param *globalState a global state for the Mersenne Twister PRNG
 *
 * @return random 3D point inside the desired prism
 *
 * Uses a Mersenne Twister PRNG and Barycentric coordinates to generate a
 * random position inside a given triangle in a specific depth
 */
__device__ Point Mesh::genRndPoint(unsigned triangle, unsigned level, curandStateMtgp32 *globalState){
  Point startPoint = {0,0,0};
  double u = curand_uniform(&globalState[blockIdx.x]);
  double v = curand_uniform(&globalState[blockIdx.x]);

  if((u+v)>1)
  {
    u = 1-u;
    v = 1-v;
  }
  double w = 1-u-v;
  int t1 = triangles[triangle];
  int t2 = triangles[triangle + numberOfTriangles];
  int t3 = triangles[triangle + 2 * numberOfTriangles];

  // convert the random startpoint into coordinates
  startPoint.z = (level + curand_uniform(&globalState[blockIdx.x])) * thickness;
  startPoint.x = (points[t1] * u) + (points[t2] * v) + (points[t3] * w);
  startPoint.y = (points[t1+numberOfPoints] * u) + (points[t2+numberOfPoints] * v) + (points[t3+numberOfPoints] * w);

  return startPoint;
}


/**
 * @brief get a betaValue for a specific triangle and level
 *
 * @param triangle the id of the desired triangle
 *
 * @param level the level of the desired prism
 *
 * @return a beta value
 */
__device__ double Mesh::getBetaValue(unsigned triangle, unsigned level){
  return betaValues[triangle + level*numberOfTriangles];
}

/**
 * @brief get a betaValue for a specific prism
 *
 * @param the id of the desired prism
 *
 * @return a beta value
 */
__device__ double Mesh::getBetaValue(unsigned prism){
  return betaValues[prism];
}

/**
 * @brief generates a normal vector for a given side of a triangle
 *
 * @param triangle the desired triangle
 *
 * @param the edge (0,1,2) of the triangle
 *
 * @return a normal vector with length 1
 */
__device__ NormalRay Mesh::getNormal(unsigned triangle, int edge){
  NormalRay ray = { {0,0},{0,0}};
  int offset =  edge*numberOfTriangles + triangle;
  ray.p.x = points[ normalPoint [offset] ];
  ray.p.y = points[ normalPoint [offset] + numberOfPoints ];

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
__device__ Point Mesh::getSamplePoint(unsigned sample){
  Point p = {0,0,0};
  unsigned level = sample/numberOfPoints;
  p.z = level*thickness;
  unsigned pos = sample - (numberOfPoints*level);
  p.x = points[pos];
  p.y = points[pos + numberOfPoints];
  return p;
}

/**
 * @brief get a Point in the center of a prism
 *
 * @param triangle the id of the desired triangle
 *
 * @param level the level of the desired prism
 *
 * @return a point with the coordinates (3D) of the prism center
 */
__device__ Point Mesh::getCenterPoint(unsigned triangle,unsigned level){
  Point p = {0,0,(level+0.5)*thickness};
  p.x = centers[triangle];
  p.y = centers[triangle + numberOfTriangles];
  return p;
}

/**
 * @brief gets the edge-id which will be forbidden
 *
 * @param trianle the index of the triangle you are currently in
 *
 * @param edge the index of the edge through which you are leaving the triangle
 *
 * @return the id of the edge, which will be forbidden in the new triangle
 *
 * The forbidden edge corresponds to the edge through which you left the
 * previous triangle (has a different index in the new triangle)
 */
__device__ int Mesh::getForbiddenEdge(unsigned triangle,int edge){
  return forbidden[edge * numberOfTriangles + triangle];
}


__device__ unsigned Mesh::getCellType(unsigned triangle){
  return cellTypes[triangle];
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
 * @param *points coordinates of the vertices in one layer of the mesh
 * 
 * @param *betaValues constant values for each meshed prism
 *
 * @param *xOfTriangleCenter the x coordinates of each triangle's center
 *
 * @param *yOfTriangleCenter the y coordinates of each triangle's center
 *
 * @param *positionsOfNormalVectors indices to the points (points), where the normals xOfNormals,yOfNormals start
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
 * @param *devices array of device indices for all possible devices 
 *
 * @param maxGpus maximal number of devices to allocate
 */
int Mesh::parseMultiGPU(Mesh& hMesh,
			std::vector<Mesh>& dMesh,
			std::string root,
			std::vector<unsigned> devices,
			unsigned maxGpus) {

  // Experimentdata
  std::vector<double> * betaValues = new std::vector<double>;
  std::vector<double> * xOfNormals = new std::vector<double>;
  std::vector<double> * yOfNormals = new std::vector<double>;
  std::vector<unsigned> * triangleIndices = new std::vector<unsigned>;
  std::vector<int> * forbidden = new std::vector<int>;
  std::vector<int> * neighbors = new std::vector<int>;
  std::vector<unsigned> * positionsOfNormalVectors = new std::vector<unsigned>;
  std::vector<double> * points = new std::vector<double>;
  std::vector<float> * surfaces = new std::vector<float>;
  std::vector<double> *xOfTriangleCenter = new std::vector<double>;
  std::vector<double> *yOfTriangleCenter = new std::vector<double>;
  std::vector<double> * betaCells = new std::vector<double>;
  std::vector<unsigned> * cellTypes = new std::vector<unsigned>;
  std::vector<float> * refractiveIndices = new std::vector<float>;
  std::vector<float> * reflectivities = new std::vector<float>;

  unsigned numberOfPoints = 0;
  unsigned numberOfTriangles = 0;
  unsigned numberOfLevels = 0;
  unsigned cladNumber;
  float thicknessOfPrism = 1;
  float nTot = 0;
  float crystalFluorescence = 0;
  double cladAbsorption = 0;


  // Parse experimentdata from files
  if(fileToVector(root + "n_p.txt", positionsOfNormalVectors)) return 1;
  if(fileToVector(root + "beta_v.txt", betaValues)) return 1;
  if(fileToVector(root + "forbidden.txt", forbidden)) return 1;
  if(fileToVector(root + "neighbors.txt", neighbors)) return 1;
  if(fileToVector(root + "n_x.txt", xOfNormals)) return 1;
  if(fileToVector(root + "n_y.txt", yOfNormals)) return 1;
  if(fileToVector(root + "x_center.txt", xOfTriangleCenter)) return 1;
  if(fileToVector(root + "y_center.txt", yOfTriangleCenter)) return 1;
  if(fileToVector(root + "p_in.txt", points)) return 1;
  if(fileToVector(root + "t_in.txt", triangleIndices)) return 1;
  if(fileToVector(root + "surface.txt", surfaces)) return 1;
  if(fileToValue(root + "size_p.txt", numberOfPoints)) return 1;
  if(fileToValue(root + "size_t.txt", numberOfTriangles)) return 1;
  if(fileToValue(root + "mesh_z.txt", numberOfLevels)) return 1;
  if(fileToValue(root + "z_mesh.txt", thicknessOfPrism)) return 1;
  if(fileToValue(root + "n_tot.txt", nTot)) return 1;
  if(fileToValue(root + "tfluo.txt", crystalFluorescence)) return 1;
  if(fileToValue(root + "clad_num.txt", cladNumber)) return 1;
  if(fileToValue(root + "clad_abs.txt", cladAbsorption)) return 1;
  if(fileToVector(root + "beta_cell.txt", betaCells)) return 1;
  if(fileToVector(root + "clad_int.txt", cellTypes)) return 1;
  if(fileToVector(root + "refractive_indices.txt", refractiveIndices)) return 1;
  if(fileToVector(root + "reflectivities.txt", reflectivities)) return 1;

  assert(numberOfPoints == (points->size() / 2));
  assert(numberOfTriangles == triangleIndices->size() / 3);
  assert(positionsOfNormalVectors->size() == numberOfTriangles * 3);
  assert(yOfTriangleCenter->size() == numberOfTriangles);
  assert(xOfTriangleCenter->size() == numberOfTriangles);
  assert(surfaces->size() == numberOfTriangles);
  assert(betaValues->size() == numberOfTriangles * (numberOfLevels-1));
  assert(xOfNormals->size() == numberOfTriangles * 3);
  assert(yOfNormals->size() == numberOfTriangles * 3);
  assert(triangleIndices->size() == numberOfTriangles * 3);
  assert(forbidden->size() == numberOfTriangles * 3);
  assert(neighbors->size() == numberOfTriangles * 3);
  assert(betaCells->size() == numberOfPoints * numberOfLevels);
  assert(cellTypes->size()== numberOfTriangles);
  assert(refractiveIndices->size() == 4);
  assert(reflectivities->size() == (refractiveIndices->size()/2) * numberOfTriangles);


  fillHMesh(
      hMesh,
      numberOfTriangles, 
      numberOfLevels,
      numberOfPoints, 
      thicknessOfPrism,
      triangleIndices,
      points,
      xOfTriangleCenter, 
      yOfTriangleCenter, 
      positionsOfNormalVectors,
      xOfNormals, 
      yOfNormals,
      forbidden, 
      neighbors, 
      surfaces,
      betaValues,
      betaCells,
      cellTypes,
	  refractiveIndices,
	  reflectivities,
      nTot,
      crystalFluorescence,
      cladNumber,
      cladAbsorption
  );

  for( unsigned i=0; i<maxGpus; i++){
    CUDA_CHECK_RETURN(cudaSetDevice(devices.at(i)) );
    fillDMesh(
        hMesh,
        dMesh.at(i),
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
        surfaces,
        betaValues,
        betaCells,
        cellTypes,
		refractiveIndices,
		reflectivities,
        nTot,
        crystalFluorescence,
        cladNumber,
        cladAbsorption
    );
    cudaDeviceSynchronize();
  }
  return 0;

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

unsigned Mesh::getMaxReflections (int reflectionPlane) const{
	double d = calculateMaxDiameter(points,numberOfPoints);
	float alpha = getReflectionAngle(reflectionPlane) * M_PI / 180.;
	double h = numberOfLevels * thickness; 
	double z = d/tan(alpha);
	return ceil(z/h);
}

unsigned Mesh::getMaxReflections() const{
	unsigned top = getMaxReflections(1);
	unsigned bottom = getMaxReflections(-1);
	return max(top,bottom);
}

__device__ __host__ float Mesh::getReflectivity(int reflectionPlane, unsigned triangle) const{
	switch(reflectionPlane){
		case -1:
			return reflectivities[triangle];
		case 1:
			return reflectivities[triangle + numberOfTriangles];
	}
	return 0;
}

__device__ __host__ float Mesh::getReflectionAngle(int reflectionPlane) const{
	switch(reflectionPlane){
		case -1:
      return totalReflectionAngles[0];
			//return asin(refractiveIndices[1]/refractiveIndices[0]);
		case 1:
      return totalReflectionAngles[1];
			//return asin(refractiveIndices[3]/refractiveIndices[2]);
	}
	return  0;
}

