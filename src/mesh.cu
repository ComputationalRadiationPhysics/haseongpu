#include <stdio.h>
#include <vector>
#include <string>
#include <assert.h>
#include <cfloat>
#include <cmath>
#include <algorithm>

#include <cudachecks.h>
#include <mesh.h>
#include <parser.h>
#include <reflection.h>


template <class T, class B, class E>
void assertRange(const std::vector<T> &v, const B minElement,const E maxElement, const bool equals){
  if(equals){
    assert(*std::min_element(v.begin(),v.end()) == minElement);
    assert(*std::max_element(v.begin(),v.end()) == maxElement);
  }else{
    assert(*std::min_element(v.begin(),v.end()) >= minElement);
    assert(*std::max_element(v.begin(),v.end()) <= maxElement);
  }
}

template <class T, class B>
void assertMin(const std::vector<T>& v,const  B minElement,const bool equals){
  if(equals){
    assert(*std::min_element(v.begin(),v.end()) == minElement);
  }else{
    assert(*std::min_element(v.begin(),v.end()) >= minElement);
  }
}

Mesh::~Mesh() {
//  if(points) delete points;
//  if(betaVolume) delete betaVolume;
//  if(normalVec) delete normalVec;
//  if(centers) delete centers;
//  if(triangleSurfaces) delete triangleSurfaces;
//  if(forbiddenEdge) delete forbiddenEdge;
//  if(betaCells) delete betaCells;
//  if(trianglePointIndices) delete trianglePointIndices;  
//  if(triangleNeighbors) delete triangleNeighbors;
//  if(triangleNormalPoint) delete triangleNormalPoint;
//  if(claddingCellTypes) delete claddingCellTypes;
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
    float thickness,
    std::vector<unsigned> *trianglePointIndices,
    std::vector<double> *points,
    std::vector<double> *triangleCenterX, 
    std::vector<double> *triangleCenterY, 
    std::vector<unsigned> *positionsOfNormal,
    std::vector<double> *triangleNormalsX, 
    std::vector<double> *triangleNormalsY,
    std::vector<int> *forbiddenEdge, 
    std::vector<int> *triangleNeighbors, 
    std::vector<float> *triangleSurfaces,
    std::vector<double> *betaVolume,
    std::vector<double> *betaCells,
    std::vector<unsigned> * claddingCellTypes,
    std::vector<float> * refractiveIndices,
    std::vector<float> * reflectivities,
    float nTot,
    float crystalTFluo,
    unsigned claddingNumber,
    double claddingAbsorption
) {

  hMesh.numberOfTriangles = numberOfTriangles;
  hMesh.numberOfLevels = numberOfLevels;
  hMesh.numberOfPrisms = numberOfTriangles*(numberOfLevels-1);
  hMesh.numberOfPoints = numberOfPoints;
  hMesh.numberOfSamples = numberOfPoints * numberOfLevels;
  hMesh.thickness = thickness;
  hMesh.crystalTFluo = crystalTFluo;
  hMesh.nTot = nTot;
  hMesh.claddingNumber = claddingNumber;
  hMesh.claddingAbsorption = claddingAbsorption;

  std::vector<double> *hostCenters = new std::vector<double>(triangleCenterX->begin(), triangleCenterX->end());
  hostCenters->insert(hostCenters->end(),triangleCenterY->begin(),triangleCenterY->end());

  std::vector<double> *hostNormalVec = new std::vector<double>(triangleNormalsX->begin(), triangleNormalsX->end());
  hostNormalVec->insert(hostNormalVec->end(),triangleNormalsY->begin(),triangleNormalsY->end());

  hMesh.points = &(points->at(0));
  hMesh.trianglePointIndices = &(trianglePointIndices->at(0));
  hMesh.betaVolume = &(betaVolume->at(0));
  hMesh.normalVec = &(hostNormalVec->at(0));
  hMesh.centers = &(hostCenters->at(0));
  hMesh.triangleSurfaces = &(triangleSurfaces->at(0));
  hMesh.forbiddenEdge = &(forbiddenEdge->at(0));
  hMesh.triangleNeighbors = &(triangleNeighbors->at(0));
  hMesh.triangleNormalPoint = &(positionsOfNormal->at(0));
  hMesh.betaCells = &(betaCells->at(0));
  hMesh.claddingCellTypes = &(claddingCellTypes->at(0));
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
    std::vector<unsigned> *trianglePointIndices, 
    unsigned numberOfTriangles, 
    unsigned numberOfLevels,
    unsigned numberOfPoints, 
    float thickness,
    std::vector<double> *pointsVector, 
    std::vector<double> *triangleCenterX, 
    std::vector<double> *triangleCenterY, 
    std::vector<unsigned> *triangleNormalPoint,
    std::vector<double> *triangleNormalsX, 
    std::vector<double> *triangleNormalsY,
    std::vector<int> *forbiddenVector, 
    std::vector<int> *triangleNeighborsVector, 
    std::vector<float> *surfacesVector,
    std::vector<double> *betaVolumeVector,
    std::vector<double> *betaCells,
    std::vector<unsigned> *claddingCellTypes,
    std::vector<float> * refractiveIndices,
    std::vector<float> * reflectivities,
    float nTot,
    float crystalTFluo,
    unsigned claddingNumber,
    double claddingAbsorption
) {


  // GPU variables
  double totalSurface = 0.;

  // constants
  dMesh.numberOfTriangles = numberOfTriangles;
  dMesh.numberOfLevels = numberOfLevels;
  dMesh.numberOfPrisms = numberOfTriangles*(numberOfLevels-1);
  dMesh.numberOfPoints = numberOfPoints;
  dMesh.numberOfSamples = numberOfPoints*numberOfLevels;
  dMesh.thickness = thickness;
  dMesh.crystalTFluo = crystalTFluo;
  dMesh.nTot = nTot;
  dMesh.claddingAbsorption = claddingAbsorption;
  dMesh.claddingNumber = claddingNumber;

  for(unsigned i=0;i<numberOfTriangles;++i){
    totalSurface+=double(surfacesVector->at(i));	
  }
  dMesh.surfaceTotal = float(totalSurface);


  // values
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.points), 2 * hMesh.numberOfPoints * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.normalVec), 2 * 3 * hMesh.numberOfTriangles * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.betaVolume), hMesh.numberOfPrisms * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.centers), 2 * hMesh.numberOfTriangles * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.triangleSurfaces), hMesh.numberOfTriangles * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.forbiddenEdge), 3 * hMesh.numberOfTriangles * sizeof(int)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.betaCells), hMesh.numberOfPoints * hMesh.numberOfLevels * sizeof(double)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.claddingCellTypes), hMesh.numberOfTriangles * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.refractiveIndices), refractiveIndices->size() * sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.reflectivities), (refractiveIndices->size()/2)*(hMesh.numberOfTriangles) *sizeof(float)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.totalReflectionAngles), (refractiveIndices->size()/2) *sizeof(float)));

  // indexStructs
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.trianglePointIndices), 3 * hMesh.numberOfTriangles * sizeof(unsigned)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.triangleNeighbors), 3 * hMesh.numberOfTriangles * sizeof(int)));
  CUDA_CHECK_RETURN(cudaMalloc(&(dMesh.triangleNormalPoint), 3 * hMesh.numberOfTriangles * sizeof(unsigned)));


  /// fill values
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.points, (double*) &(pointsVector->at(0)), 2 * hMesh.numberOfPoints * sizeof(double), cudaMemcpyHostToDevice));

  std::vector<double> *hostNormalVec = new std::vector<double>(triangleNormalsX->begin(), triangleNormalsX->end());
  hostNormalVec->insert(hostNormalVec->end(),triangleNormalsY->begin(),triangleNormalsY->end());
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.normalVec, (double*) &(hostNormalVec->at(0)), 2 * 3 * hMesh.numberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
  free(hostNormalVec);

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.betaVolume, (double*) &(betaVolumeVector->at(0)), hMesh.numberOfPrisms * sizeof(double), cudaMemcpyHostToDevice));

  std::vector<double> *hostCenters = new std::vector<double>(triangleCenterX->begin(), triangleCenterX->end());
  hostCenters->insert(hostCenters->end(),triangleCenterY->begin(),triangleCenterY->end());
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.centers, (double*) &(hostCenters->at(0)), 2 * hMesh.numberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
  free(hostCenters);

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.triangleSurfaces, (float*) &(surfacesVector->at(0)), hMesh.numberOfTriangles * sizeof(float), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.forbiddenEdge, (int*) &(forbiddenVector->at(0)), 3 * hMesh.numberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.betaCells, (double*) &(betaCells->at(0)), hMesh.numberOfPoints * hMesh.numberOfLevels * sizeof(double), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.claddingCellTypes, (unsigned*) &(claddingCellTypes->at(0)), hMesh.numberOfTriangles * sizeof(unsigned), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.refractiveIndices, (float*) &(refractiveIndices->at(0)), refractiveIndices->size() * sizeof(float), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.reflectivities, (float*) &(reflectivities->at(0)), (hMesh.numberOfTriangles)* (refractiveIndices->size()/2) * sizeof(float), cudaMemcpyHostToDevice));

  std::vector<float> *totalReflectionAngles = new std::vector<float>(refractiveIndices->size()/2,0);
  for(unsigned i=0;i<refractiveIndices->size();i+=2){
    totalReflectionAngles->at(i/2) = (180. / M_PI *  asin(refractiveIndices->at(i+1) / refractiveIndices->at(i)));
  }
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.totalReflectionAngles, (float*) &(totalReflectionAngles->at(0)), refractiveIndices->size()/2 * sizeof(float), cudaMemcpyHostToDevice));
  free(totalReflectionAngles);

  // fill indexStructs
  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.trianglePointIndices, (unsigned*) &(trianglePointIndices->at(0)), 3 * hMesh.numberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.triangleNeighbors,(int*) &(triangleNeighborsVector->at(0)), 3 * hMesh.numberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));

  CUDA_CHECK_RETURN(cudaMemcpy(dMesh.triangleNormalPoint, (unsigned*) &(triangleNormalPoint->at(0)), 3 * hMesh.numberOfTriangles * sizeof(unsigned), cudaMemcpyHostToDevice));

}

/**
 * @brief fetch the id of an adjacent triangle
 *
 * @param triangle the index of the triangle of which you want the neighbor
 * @param edge the side of the triangle, for whih you want the neighbor
 *
 * @return the index of the neighbor triangle
 */
__device__ int Mesh::getNeighbor(unsigned triangle, int edge) const{
  return triangleNeighbors[triangle + edge*numberOfTriangles];
}

/**
 * @brief generate a random point within a prism
 *
 * @param triangle the triangle to describe the desired prism
 * @param the level of the desired prism
 * @param *globalState a global state for the Mersenne Twister PRNG
 *
 * @return random 3D point inside the desired prism
 *
 * Uses a Mersenne Twister PRNG and Barycentric coordinates to generate a
 * random position inside a given triangle in a specific depth
 */
__device__ Point Mesh::genRndPoint(unsigned triangle, unsigned level, curandStateMtgp32 *globalState) const{
  Point startPoint = {0,0,0};
  double u = curand_uniform_double(&globalState[blockIdx.x]);
  double v = curand_uniform_double(&globalState[blockIdx.x]);

  if((u+v)>1)
  {
    u = 1-u;
    v = 1-v;
  }
  double w = 1-u-v;
  int t1 = trianglePointIndices[triangle];
  int t2 = trianglePointIndices[triangle + numberOfTriangles];
  int t3 = trianglePointIndices[triangle + 2 * numberOfTriangles];

  // convert the random startpoint into coordinates
  startPoint.z = (level + curand_uniform_double(&globalState[blockIdx.x])) * thickness;
  startPoint.x = (points[t1] * u) + (points[t2] * v) + (points[t3] * w);
  startPoint.y = (points[t1+numberOfPoints] * u) + (points[t2+numberOfPoints] * v) + (points[t3+numberOfPoints] * w);

  return startPoint;
}


/**
 * @brief get a betaVolume for a specific triangle and level
 *
 * @param triangle the id of the desired triangle
 * @param level the level of the desired prism
 *
 * @return a beta value
 */
__device__ double Mesh::getBetaVolume(unsigned triangle, unsigned level) const{
  unsigned i = triangle + level * numberOfTriangles;
    if(i > 7307){
      printf("nT: %d T: %d L: %d B: %d\n", numberOfTriangles, triangle, level, i);
    }

  return betaVolume[triangle + level * numberOfTriangles];
}

/**
 * @brief get a betaVolume for a specific prism
 *
 * @param the id of the desired prism
 *
 * @return a beta value
 */
__device__ double Mesh::getBetaVolume(unsigned prism) const{
  return betaVolume[prism];
}

/**
 * @brief generates a normal vector for a given side of a triangle
 *
 * @param triangle the desired triangle
 * @param the edge (0,1,2) of the triangle
 *
 * @return a normal vector with length 1
 */
__device__ NormalRay Mesh::getNormal(unsigned triangle, int edge) const{
  NormalRay ray = { {0,0},{0,0}};
  int offset =  edge*numberOfTriangles + triangle;
  ray.p.x = points[ triangleNormalPoint [offset] ];
  ray.p.y = points[ triangleNormalPoint [offset] + numberOfPoints ];

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
__device__ Point Mesh::getSamplePoint(unsigned sample_i) const{
  Point p = {0,0,0};
  unsigned level = sample_i/numberOfPoints;
  p.z = level*thickness;
  unsigned pos = sample_i - (numberOfPoints * level);
  p.x = points[pos];
  p.y = points[pos + numberOfPoints];
  return p;
}

/**
 * @brief get a Point in the center of a prism
 *
 * @param triangle the id of the desired triangle
 * @param level the level of the desired prism
 *
 * @return a point with the coordinates (3D) of the prism center
 */
__device__ Point Mesh::getCenterPoint(unsigned triangle,unsigned level) const{
  Point p = {0,0,(level+0.5)*thickness};
  p.x = centers[triangle];
  p.y = centers[triangle + numberOfTriangles];
  return p;
}

/**
 * @brief gets the edge-id which will be forbidden
 *
 * @param trianle the index of the triangle you are currently in
 * @param edge the index of the edge through which you are leaving the triangle
 *
 * @return the id of the edge, which will be forbidden in the new triangle
 *
 * The forbidden edge corresponds to the edge through which you left the
 * previous triangle (has a different index in the new triangle)
 */
__device__ int Mesh::getForbiddenEdge(unsigned triangle,int edge) const{
  return forbiddenEdge[edge * numberOfTriangles + triangle];
}


__device__ unsigned Mesh::getCellType(unsigned triangle) const{
  return claddingCellTypes[triangle];
}


/**
 * @brief creates the Mesh datastructures on the host and on all possible devices for the propagation
 *
 * @param *hMesh the host mesh
 * @param **dMesh an array of device meshes (one for each device) 
 * @param *trianglePointIndices indices of the points which form a triangle
 * @param numberOfTriangles the number of trianglePointIndices
 * @param numberOfLeves the number of layers of the mesh
 * @param numberOfPoints the number of vertices in one layer of the mesh
 * @param thickness  the thickness of one layer of the mesh
 * @param *points coordinates of the vertices in one layer of the mesh
 * @param *betaVolume constant values for each meshed prism
 * @param *triangleCenterX the x coordinates of each triangle's center
 * @param *triangleCenterY the y coordinates of each triangle's center
 * @param *triangleNormalPoint indices to the points (points), where the normals triangleNormalsX,triangleNormalsY start
 * @param *triangleNormalsX the x components of a normal vector for each of the 3 sides of a triangle
 * @param *triangleNormalsY the y components of a normal vector for each of the 3 sides of a triangle
 * @param *forbiddenEdge the sides of the triangle from which a ray "entered" the triangle
 * @param *triangleNeighbors indices to the adjacent trianglePointIndices in trianglePointIndices
 * @param *triangleSurfaces the sizes of the surface of each prism
 * @param *devices array of device indices for all possible devices 
 * @param maxGpus maximal number of devices to allocate
 */
int Mesh::parseMultiGPU(Mesh& hMesh,
			std::vector<Mesh>& dMesh,
			std::string root,
			std::vector<unsigned> devices,
			unsigned maxGpus) {

  // Experimentdata
  std::vector<double> * betaVolume = new std::vector<double>;
  std::vector<double> * triangleNormalsX = new std::vector<double>;
  std::vector<double> * triangleNormalsY = new std::vector<double>;
  std::vector<unsigned> * trianglePointIndices = new std::vector<unsigned>;
  std::vector<int> * forbiddenEdge = new std::vector<int>;
  std::vector<int> * triangleNeighbors = new std::vector<int>;
  std::vector<unsigned> * triangleNormalPoint = new std::vector<unsigned>;
  std::vector<double> * points = new std::vector<double>;
  std::vector<float> * triangleSurfaces = new std::vector<float>;
  std::vector<double> *triangleCenterX = new std::vector<double>;
  std::vector<double> *triangleCenterY = new std::vector<double>;
  std::vector<double> * betaCells = new std::vector<double>;
  std::vector<unsigned> * claddingCellTypes = new std::vector<unsigned>;
  std::vector<float> * refractiveIndices = new std::vector<float>;
  std::vector<float> * reflectivities = new std::vector<float>;

  unsigned numberOfPoints = 0;
  unsigned numberOfTriangles = 0;
  unsigned numberOfLevels = 0;
  unsigned claddingNumber;
  float thickness = 1;
  float nTot = 0;
  float crystalTFluo = 0;
  double claddingAbsorption = 0;


  // Parse experimentdata from files
  //if(fileToVector(root + "n_p.txt", triangleNormalPoint)) return 1;
  //if(fileToVector(root + "beta_v.txt", betaVolume)) return 1;
  //if(fileToVector(root + "forbidden.txt", forbidden)) return 1;
  //if(fileToVector(root + "triangleNeighbors.txt", triangleNeighbors)) return 1;
  //if(fileToVector(root + "n_x.txt", triangleNormalsX)) return 1;
  //if(fileToVector(root + "n_y.txt", triangleNormalsY)) return 1;
  //if(fileToVector(root + "x_center.txt", triangleCenterX)) return 1;
  //if(fileToVector(root + "y_center.txt", triangleCenterY)) return 1;
  //if(fileToVector(root + "p_in.txt", points)) return 1;
  //if(fileToVector(root + "t_in.txt", trianglePointIndices)) return 1;
  //if(fileToVector(root + "surface.txt", surfaces)) return 1;
  //if(fileToValue(root + "size_p.txt", numberOfPoints)) return 1;
  //if(fileToValue(root + "size_t.txt", numberOfTriangles)) return 1;
  //if(fileToValue(root + "mesh_z.txt", numberOfLevels)) return 1;
  //if(fileToValue(root + "z_mesh.txt", thickness)) return 1;
  //if(fileToValue(root + "n_tot.txt", nTot)) return 1;
  //if(fileToValue(root + "tfluo.txt", crystalTFluo)) return 1;
  //if(fileToValue(root + "clad_num.txt", claddingNumber)) return 1;
  //if(fileToValue(root + "clad_abs.txt", claddingAbsorption)) return 1;
  //if(fileToVector(root + "beta_cell.txt", betaCells)) return 1;
  //if(fileToVector(root + "clad_int.txt", claddingCellTypes)) return 1;
  //if(fileToVector(root + "refractive_indices.txt", refractiveIndices)) return 1;
  //if(fileToVector(root + "reflectivities.txt", reflectivities)) return 1;

  //TODO
 if(fileToVector(root + "triangleNormalPoint.txt", triangleNormalPoint)) return 1;
 if(fileToVector(root + "betaVolume.txt", betaVolume)) return 1;
 if(fileToVector(root + "forbiddenEdge.txt", forbiddenEdge)) return 1;
 if(fileToVector(root + "triangleNeighbors.txt", triangleNeighbors)) return 1;
 if(fileToVector(root + "triangleNormalsX.txt", triangleNormalsX)) return 1;
 if(fileToVector(root + "triangleNormalsY.txt", triangleNormalsY)) return 1;
 if(fileToVector(root + "triangleCenterX.txt", triangleCenterX)) return 1;
 if(fileToVector(root + "triangleCenterY.txt", triangleCenterY)) return 1;
 if(fileToVector(root + "points.txt", points)) return 1;
 if(fileToVector(root + "trianglePointIndices.txt", trianglePointIndices)) return 1;
 if(fileToVector(root + "triangleSurfaces.txt", triangleSurfaces)) return 1;
 if(fileToValue(root + "numberOfPoints.txt", numberOfPoints)) return 1;
 if(fileToValue(root + "numberOfTriangles.txt", numberOfTriangles)) return 1;
 if(fileToValue(root + "numberOfLevels.txt", numberOfLevels)) return 1;
 if(fileToValue(root + "thickness.txt", thickness)) return 1;
 if(fileToValue(root + "nTot.txt", nTot)) return 1;
 if(fileToValue(root + "crystalTFluo.txt", crystalTFluo)) return 1;
 if(fileToValue(root + "claddingNumber.txt", claddingNumber)) return 1;
 if(fileToValue(root + "claddingAbsorption.txt", claddingAbsorption)) return 1;
 if(fileToVector(root + "betaCells.txt", betaCells)) return 1;
 if(fileToVector(root + "claddingCellTypes.txt", claddingCellTypes)) return 1;
 if(fileToVector(root + "refractiveIndices.txt", refractiveIndices)) return 1;
 if(fileToVector(root + "reflectivities.txt", reflectivities)) return 1;

  // assert input sizes
  assert(numberOfPoints == (points->size() / 2));
  assert(numberOfTriangles == trianglePointIndices->size() / 3);
  assert(triangleNormalPoint->size() == numberOfTriangles * 3);
  assert(triangleCenterY->size() == numberOfTriangles);
  assert(triangleCenterX->size() == numberOfTriangles);
  assert(triangleSurfaces->size() == numberOfTriangles);
  assert(betaVolume->size() == numberOfTriangles * (numberOfLevels-1));
  assert(triangleNormalsX->size() == numberOfTriangles * 3);
  assert(triangleNormalsY->size() == numberOfTriangles * 3);
  assert(trianglePointIndices->size() == numberOfTriangles * 3);
  assert(forbiddenEdge->size() == numberOfTriangles * 3);
  assert(triangleNeighbors->size() == numberOfTriangles * 3);
  assert(betaCells->size() == numberOfPoints * numberOfLevels);
  assert(claddingCellTypes->size()== numberOfTriangles);
  assert(refractiveIndices->size() == 4);
  assert(reflectivities->size() == (refractiveIndices->size()/2) * numberOfTriangles);
  assert(claddingCellTypes->size() == numberOfTriangles);

  // assert input data validity
  assertRange(*triangleNormalPoint,0u,unsigned(numberOfPoints-1),true);
  assertMin(*betaVolume,0,false);
  assertRange(*forbiddenEdge,-1,2,true);
  assertRange(*triangleNeighbors,-1,int(numberOfTriangles-1),true);
  assertRange(*triangleNormalsX,-1,1,false);
  assertRange(*triangleNormalsY,-1,1,false);
  assertRange(*triangleCenterX,*std::min_element(points->begin(),points->end()),*std::max_element(points->begin(),points->end()),false);
  assertRange(*triangleCenterY,*std::min_element(points->begin(),points->end()),*std::max_element(points->begin(),points->end()),false);
  assertRange(*trianglePointIndices,0u,unsigned(numberOfPoints-1),true);
  assertMin(*triangleSurfaces,0,false);
  assertMin(*betaCells,0,false);
  assertRange(*refractiveIndices,0,5,false);
  assertRange(*reflectivities,0,1,false);

  fillHMesh(
      hMesh,
      numberOfTriangles, 
      numberOfLevels,
      numberOfPoints, 
      thickness,
      trianglePointIndices,
      points,
      triangleCenterX, 
      triangleCenterY, 
      triangleNormalPoint,
      triangleNormalsX, 
      triangleNormalsY,
      forbiddenEdge, 
      triangleNeighbors, 
      triangleSurfaces,
      betaVolume,
      betaCells,
      claddingCellTypes,
	  refractiveIndices,
	  reflectivities,
      nTot,
      crystalTFluo,
      claddingNumber,
      claddingAbsorption
  );

  for( unsigned i=0; i < maxGpus; i++){
    CUDA_CHECK_RETURN(cudaSetDevice(devices.at(i)) );
    fillDMesh(hMesh,
	      dMesh.at(i),
	      trianglePointIndices,
	      numberOfTriangles,
	      numberOfLevels,
	      numberOfPoints,
	      thickness,
	      points,
	      triangleCenterX,
	      triangleCenterY,
	      triangleNormalPoint,
	      triangleNormalsX,
	      triangleNormalsY,
	      forbiddenEdge,
	      triangleNeighbors,
	      triangleSurfaces,
	      betaVolume,
	      betaCells,
	      claddingCellTypes,
	      refractiveIndices,
	      reflectivities,
	      nTot,
	      crystalTFluo,
	      claddingNumber,
	      claddingAbsorption
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

unsigned Mesh::getMaxReflections (ReflectionPlane reflectionPlane) const{
	double d = calculateMaxDiameter(points,numberOfPoints);
	float alpha = getReflectionAngle(reflectionPlane) * M_PI / 180.;
	double h = numberOfLevels * thickness; 
	double z = d/tan(alpha);
	return ceil(z/h);
}

unsigned Mesh::getMaxReflections() const{
	unsigned top = getMaxReflections(TOP_REFLECTION);
	unsigned bottom = getMaxReflections(BOTTOM_REFLECTION);
	return max(top,bottom);
}

__device__ __host__ float Mesh::getReflectivity(ReflectionPlane reflectionPlane, unsigned triangle) const{
	switch(reflectionPlane){
		case BOTTOM_REFLECTION:
			return reflectivities[triangle];
		case TOP_REFLECTION:
			return reflectivities[triangle + numberOfTriangles];
	}
	return 0;
}

__device__ __host__ float Mesh::getReflectionAngle(ReflectionPlane reflectionPlane) const{
	switch(reflectionPlane){
		case BOTTOM_REFLECTION:
      return totalReflectionAngles[0];
			//return asin(refractiveIndices[1]/refractiveIndices[0]);
		case TOP_REFLECTION:
      return totalReflectionAngles[1];
			//return asin(refractiveIndices[3]/refractiveIndices[2]);
	}
	return  0;
}

