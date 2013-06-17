//#include "ray_propagation_gpu.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include <vector>
#include "curand_kernel.h"
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>

#define TEST_VALUES true
#define USE_IMPORTANCE true
#define SMALL 1E-06
#define VERY_SMALL 0.0

#define CUDA_CHECK_RETURN(value) {					\
    cudaError_t _m_cudaStat = value;					\
    if (_m_cudaStat != cudaSuccess) {					\
      fprintf(stderr, "Error %s at line %d in file %s\n",		\
	      cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);	\
      exit(1);								\
    }									\
  }
#define CUDA_CALL(x) do { if((x) != cudaSuccess) {	\
      printf("Error at %s:%d\n",__FILE__,__LINE__);	\
      return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x) != CURAND_STATUS_SUCCESS) {	\
      printf("Error at %s:%d\n",__FILE__,__LINE__);		\
      return EXIT_FAILURE;}} while(0)

__device__ double cladAbsorption;
__device__ double nTot;
__device__ double sigmaE;
__device__ double sigmaA;
__device__ double thicknessOfPrism;
__device__ int numberOfLevels;
__device__ int cladNumber;
__device__ int numberOfPoints;
__device__ int numberOfTriangles;

/**
 * @brief Propagate a ray between 2 points and calculate the resulting ASE-Flux at the Destination
 *
 * @params x_pos		the x-coordinate where the ray starts
 *         y_pos		the y-coordinate where the ray starts
 *         z_pos		the z-coordinate where the ray starts
 *         x_dest		the destination of the ray (x-coordinate)
 *         y_dest		the destination of the ray (y-coordinate)
 *         z_dest		the destination of the ray (z-coordinate)
 *         t_start		the index of the triangle, in which the ray starts
 *         mesh_start	the level of the mesh (the slice) in which the ray starts
 *		   p_in			coordinates of the sample-points of one layer (first all x-coordinates, then all y-coordinates)
 *		   n_x			x-coordinates for the normal-vectors for the 3 rectangular sides of each prism
 *		   n_y			y-coordinates for the normal-vectors for the 3 rectangular sides of each prism
 *		   n_p			indices of the points where the normal-vectors start	
 *		   neighbors	indices of the adjacent triangles	
 *		   forbidden	sides of the new triangles which are forbidden, after coming from an adjacent triangle
 *		   cell_type	contains the material-constant for each cell/prism
 *		   beta_v		contains the beta-values for each cell/prism
 *
 */
__device__ double rayPropagationGpu(
				    double xPos,
				    double yPos,
				    double zPos,
				    double xDestination,
				    double yDestination,
				    double zDestination,
				    int firstTriangle,
				    int firstLevel,
				    double *points,
				    double *xOfNormals,
				    double *yOfNormals,
				    int *positionsOfNormalVectors,
				    int *neighbors,
				    int *forbidden,
				    int* cellTypes,
				    double* betaValues){
  //    no reflections
  //
  //    create the vector between both points and calculate which surface would be the shortest to reach.
  //    then get the length, make the integration, get the information about the next cell out of the array
  //    set the point to the surface (this surface is "forbidden" in the calculations)
  //    proceed until you hit a the point or the surface
	
  double xVec, yVec,zVec;
  double distanceRemaining, length, lengthHelp, distanceTotal;
  double nominator, denominator;
  double gain=1.;
  int triangleCurrent, levelCurrent; // the current triangle number and position concerning the z's
  int triangleNext, levelNext, forbiddenCurrent, forbiddenNext;
  int offset;
#if TEST_VALUES==true
  double testDistance = 0;
  int loopbreaker = 0;
#endif


  //    initial positions
  triangleCurrent = firstTriangle;
  levelCurrent = firstLevel;

  // direction-vector (without reflections)
  xVec = (xDestination - xPos);
  yVec = (yDestination - yPos);
  zVec = (zDestination - zPos);

  // total distance to travel
  distanceTotal = sqrt(xVec*xVec+yVec*yVec+zVec*zVec);
  // normalized direction-vector
  xVec = xVec/distanceTotal;
  yVec = yVec/distanceTotal;
  zVec = zVec/distanceTotal;

  // remaining distance to travel
  distanceRemaining = distanceTotal;

  // at the beginning, all surfaces are possible
  forbiddenCurrent = -1;

  for(;;)
    {
      // the length of the ray-part inside the current prism. We try to minimize this value
      length = distanceRemaining;
      lengthHelp=0;
      //        definition for decider
      //        0,1,2: int for the neighbors
      //        3: hor plane up
      //        4: hor plane down
      //        try the triangle faces
      //        remember the correlation between the normals and the points
      //        n1: p1-2, n2: p1-3, n3:p2-3
      //        the third coordinate (z) of the particpating points for the surfaces can be set to be z=0, 
      //        as everything uses triangular "prisms", as well as n_z=0 in this case!
		
      // forb describes the surface, from which the ray enters the prism.
      // this surface is no suitable candidate, since the length would be 0!
      if (forbiddenCurrent != 0){
	denominator = xOfNormals[triangleCurrent]*xVec + yOfNormals[triangleCurrent]*yVec;
	// see if we intersect at all
	if (denominator != 0.0)
	  {
	    nominator = (xOfNormals[triangleCurrent]*points[positionsOfNormalVectors[triangleCurrent]] + yOfNormals[triangleCurrent]*points[positionsOfNormalVectors[triangleCurrent]+ numberOfPoints]) - (xOfNormals[triangleCurrent]*xPos + yOfNormals[triangleCurrent]*yPos);
	    lengthHelp = nominator/denominator;
	    // if we found a new smallest length, use it
	    if (lengthHelp < length && lengthHelp > 0.0)
	      {
		length = lengthHelp;
		forbiddenNext = (forbidden[triangleCurrent]);
		triangleNext = neighbors[triangleCurrent];
		levelNext = levelCurrent;

	      }
	  }
      }

      // see forbiddenCurrent !=0 case
      if (forbiddenCurrent != 1){
	//offset, since the 3 rectangular surfaces are stored at different positions in the array
	offset = triangleCurrent+numberOfTriangles;
	denominator = xOfNormals[offset]*xVec + yOfNormals[offset]*yVec;
	if (denominator != 0.0)
	  {
	    nominator = (xOfNormals[offset]*points[positionsOfNormalVectors[offset]] + yOfNormals[offset]*points[positionsOfNormalVectors[offset]+ numberOfPoints]) - (xOfNormals[offset]*xPos + yOfNormals[offset]*yPos);
	    lengthHelp = nominator/denominator;
	    if (lengthHelp < length && lengthHelp > 0.0)
	      {
		length = lengthHelp;
		forbiddenNext = (forbidden[offset]);
		triangleNext = neighbors[offset];
		levelNext = levelCurrent;
	      }
	  }
      }

      // see forbiddenCurrent !=0 case
      if (forbiddenCurrent !=2){
	offset = triangleCurrent+2*numberOfTriangles;
	denominator = xOfNormals[offset]*xVec + yOfNormals[offset]*yVec;
	if (denominator != 0.0)
	  {
	    nominator = (xOfNormals[offset]*points[positionsOfNormalVectors[offset]] + yOfNormals[offset]*points[positionsOfNormalVectors[offset]+ numberOfPoints]) - (xOfNormals[offset]*xPos + yOfNormals[offset]*yPos);
	    lengthHelp = nominator/denominator;
	    if (lengthHelp < length && lengthHelp > 0.0)
	      {
		length = lengthHelp;
		forbiddenNext = (forbidden[offset]);
		triangleNext = neighbors[offset];
		levelNext = levelCurrent;
	      }
	  }
      }

      // if-structure "optimized"
      denominator = zPos*zVec;
      if (denominator != 0.0){
	if (forbiddenCurrent != 3){
	  {
	    nominator = (levelCurrent+1)* thicknessOfPrism - zPos;
	    lengthHelp = nominator/denominator;
	    if (lengthHelp < length && lengthHelp > 0.0)
	      {
		length = lengthHelp;
		//decider = 3;
		forbiddenNext = 4; // you are not allowed to go down in the next step
		triangleNext = triangleCurrent;
		levelNext = levelCurrent + 1;
	      }
	  }
	}

	// next is the lower plane
	if (forbiddenCurrent != 4){
	  nominator = (levelCurrent)* thicknessOfPrism - zPos;
	  lengthHelp = nominator/denominator;
	  if (lengthHelp < length && lengthHelp > 0.0)
	    {
	      length = lengthHelp;
	      //decider = 4;
	      forbiddenNext = 3; // you are not allowed to go up in the next step
	      triangleNext = triangleCurrent;
	      levelNext = levelCurrent - 1;
	    }
	}
      }

      if (cellTypes[triangleCurrent] == cladNumber){
	gain *= exp((-1)*(cladAbsorption * length));
      }
      else {
	gain *= (double) exp(nTot * (betaValues[triangleCurrent+levelCurrent*numberOfTriangles]*(sigmaE + sigmaA)-sigmaA)*length);
      }


      // the remaining distance is decreased by the length we travelled through the prism
      distanceRemaining -= length;
		

#if TEST_VALUES==true
      testDistance += length;
      if(loopbreaker>500){
	printf("Loopbreaker reached. firstTriangle: %d, level: %d, length: %f, distanceTotal:%f, testDistance%f, distanceRemaining:%f\n",firstTriangle,firstLevel,length,distanceTotal,testDistance,distanceRemaining);
	return 0.;
      }else{
	loopbreaker++;
      }
#endif
      // if the distance between the destination and our current position is small enough, we are done
      if (fabs(distanceRemaining) < SMALL){
	break;
      }

      // now set the next cell and position
      xPos = xPos + length*xVec;
      yPos = yPos + length*yVec;
      zPos = zPos + length*zVec;

      triangleCurrent = triangleNext;
      levelCurrent = levelNext;
      // set the new forbidden surface
      forbiddenCurrent = forbiddenNext;

    }

#if TEST_VALUES==true
  if(fabs(distanceTotal-testDistance) > SMALL)
    printf("Distance too big! firstTriangle: %d, level: %d, length: %f, distanceTotal:%f, testDistance%f, distanceRemaining:%f\n",firstTriangle,firstLevel,length,distanceTotal,testDistance,distanceRemaining);
#endif
	
  return gain /= (distanceTotal*distanceTotal);
}



/**
 * Initializes the global variables of the GPU with the correct values.
 * All those values are from the original propagation-function which we ported.
 */
__global__ void setupGlobalVariablesKernel ( 
					    double hostSigmaE,
					    double hostSigmaA, 
					    int hostCladNum, 
					    double hostCladAbs, 
					    double hostNTot, 
					    int hostNumberOfTriangles, 
					    double hostThicknessOfPrism, 
					    int hostNumberOfLevels, 
					    int hostNumberOfPoints )
{
  sigmaE = hostSigmaE;	
  sigmaA = hostSigmaA;
  cladNumber = hostCladNum;
  cladAbsorption = hostCladAbs;
  nTot = hostNTot;
  numberOfTriangles = hostNumberOfTriangles;
  thicknessOfPrism = hostThicknessOfPrism;
  numberOfLevels = hostNumberOfLevels;
  numberOfPoints = hostNumberOfPoints;
  //printf("Sigma_e in setup=%f\tSigma_eHost=%f\n",sigma_e,host_sigma_e);
} 

// __global__ void importanceKernel(
// 		curandState *globalState,
// 		double *points,
// 		double *xOfNormals,
// 		double *yOfNormals,
// 		int *positionsOfNormalVectors,
// 		int *neighbors,
// 		int *forbidden,
// 		int* cell_type,
// 		int hostNumberOfTriangles,
// 		double* betaValues,
// 		double *importance,
// 		int *numberOfImportantRays,
// 		double *xOfTriangleCenter,
// 		double *yOfTriangleCenter,
// 		int *surface,
// 		int totalNumberOfRays) {

// 	int id = threadIdx.x + blockIdx.x * blockDim.x;
// 	for(int i=0; i< hostNumberOfTriangles; ++i){
// 		for(int j=0; j< numberOfLevels; ++j){
// 			importf(globalState[id], i,j, importance, numberOfImportantRays, points, xOfNormals, yOfNormals, positionsOfNormalVectors, neighbors, forbidden, cell_type, betaValues, xOfTriangleCenter, yOfTriangleCenter,surface, totalNumberOfRays);

// 		}

// 	}

// }


/**
 * Prints some of the global device variables.
 * Is only used for testing
 */
__global__ void testKernel (
			    double *points,
			    double *xOfNormals,
			    double *yOfNormals,
			    int *positionsOfNormalVectors,
			    int *neighbors,
			    int *forbidden,
			    int* triangleIndices,
			    int* cellTypes,
			    double* betaValues,
			    int *surfacesNormalized){
  printf("\nSigmaE=%.6e",sigmaE);
  printf("\nSigmaA=%.6e",sigmaA);
  printf("\nNumberOfLevels=%d",numberOfLevels);
  printf("\nNumberOfPoints=%d",numberOfPoints);
  printf("\nthicknessOfPrism_=%.6e",thicknessOfPrism);
  printf("\nnumberOfTriangles=%d",numberOfTriangles);
  printf("\nnTot=%.6e",nTot);
  printf("\ncladAbsorption=%.6e",cladAbsorption);
  printf("\ncladNumber=%d\n\n",cladNumber);

  unsigned limit = 5;
  for(int i=0;i<limit;++i){
    printf("points[%d]: %e\n",i,points[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("xOfNormals[%d]: %e\n",i,xOfNormals[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("yOfNormals[%d]: %e\n",i,yOfNormals[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("positionsOfNormalVectors[%d]: %d\n",i,positionsOfNormalVectors[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("neighbors[%d]: %d\n",i,neighbors[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("forbidden[%d]: %d\n",i,forbidden[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("triangleIndices[%d]: %d\n",i,triangleIndices[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("cellTypes[%d]: %d\n",i,cellTypes[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("betaValues[%d]: %e\n",i,betaValues[i]);
  }
  printf("\n\n");

  for(int i=0;i<limit;++i){
    printf("surfacesNormalized[%d]: %d\n",i,surfacesNormalized[i]);
  }
  printf("\n\n");
} 


/**
 * Does the raytracing for a single Sample point (in a defined level).
 * This Kernel has to be started for each sample point with the same value for iterations
 * and the same number of blocks/threads.
 *
 * \var globalState the state of the mersenneTwister PRNG
 * 		(has a maximum of 200 positions!)
 * \var phi points to a memory region which is initialized with 0
 * 		(can hold one value for each sample point)
 * \var point2D the index of the current sample point (points to p_in)
 * \var level the level of the current sample point (how deep we are through the material)
 * \var raysPerThread the number rays which are computed by this thread
 * 		(always for the same combination of startprism+samplepoint
 */
__global__ void raytraceStep(
			     curandStateMtgp32* globalState,
			     float* phiASE,
			     const int point2D,
			     const int level,
			     const int raysPerThread,
			     double *points,
			     double *xOfNormals,
			     double *yOfNormals,
			     int *positionsOfNormalVectors,
			     int *neighbors,
			     int *forbidden,
			     int* triangleIndices,
			     int* cellTypes,
			     double* betaValues,
			     double* importance,
			     int* numberOfImportantRays,
			     int* surfacesNormalized) {

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  const unsigned numberOfPrisms = (numberOfTriangles * (numberOfLevels-1));
  const unsigned threadsPerPrism = blockDim.x * gridDim.x / numberOfPrisms;
  // break, if we have more threads than we need
  if(id >= threadsPerPrism * numberOfPrisms)
    return;

  double gain = 0.;
  const int endPointX = points[point2D];
  const int endPointY = points[ numberOfPoints + point2D];
  const int endPointZ = level* thicknessOfPrism;


  // this should give the same start values multiple times (so that every thread uses the same prism, which yields
  // big benefits for the memory access (and caching!)
  unsigned startPrism = id % numberOfPrisms;
  int startLevel = (startPrism)/numberOfTriangles;
  int startTriangle = (startPrism-(numberOfTriangles*startLevel));

#if TEST_VALUES==true
  if(startPrism != (startTriangle+(startLevel*numberOfTriangles))){
    printf("StartTriangle/StartLevel incorrect!");
  }
  if(startTriangle >= 600){
	
    printf("StartTriangle/StartLevel incorrect!");
  }
  //if(startPrism == 5399){
  //	printf("startprism: %d, id=%d, threadsPerPrism=%d\n",startPrism,id,threadsPerPrism);
  //}
#endif

  // the indices of the vertices of the starttriangle
  int t1 = triangleIndices[startTriangle];
  int t2 = triangleIndices[startTriangle+ numberOfTriangles];
  int t3 = triangleIndices[startTriangle+2*numberOfTriangles];

  // do all this multiple times (we can't have more than 200 blocks due to restrictions of the Mersenne Twister)
  for (int i=0; i < numberOfImportantRays[startPrism]; ++i){
    //for (int i=0; i<blah; ++i){
    // random startpoint generation
    double u = curand_uniform(&globalState[blockIdx.x]);
    double v = curand_uniform(&globalState[blockIdx.x]);

    if((u+v)>1)
      {
	u = 1-u;
	v = 1-v;
      }
    double w = 1-u-v;

    // convert the random startpoint into coordinates
    double zRand = (startLevel + curand_uniform(&globalState[blockIdx.x]))* thicknessOfPrism;
    double xRand = points[t1]*u + points[t2]*v + points[t3]*w;
    double yRand = points[ numberOfPoints + t1]*u + points[ numberOfPoints + t2]*v + points[ numberOfPoints + t3]*w;

    __syncthreads();
    gain += rayPropagationGpu(xRand, yRand, zRand, endPointX, endPointY, endPointZ, 
			      startTriangle, startLevel ,points, xOfNormals, yOfNormals, 
			      positionsOfNormalVectors, neighbors, forbidden , cellTypes, betaValues);
    //gain += double(propagationOld(xRand, yRand, zRand, endPointX, endPointY, endPointZ, startTriangle, startLevel ,points, xOfNormals, yOfNormals, positionsOfNormalVectors, neighbors, forbidden , cellTypes, betaValues));
  }
	

  // do the multiplication just at the end of all raysPerThread
  // (gives better numeric behaviour)
  gain *= betaValues[startPrism];///surfacesNormalized[startTriangle];
#if USE_IMPORTANCE==true
  atomicAdd(&(phiASE[point2D + level*numberOfPoints]), float(gain * importance[startPrism]));
#else
  atomicAdd(&(phiASE[point2D + level*numberOfPoints]),float(gain));
#endif
  return;
}

/*********************************************************************************************
 * HOST FUNCTIONS
 *********************************************************************************************/
double rayPropagationCpu(double x_pos, 
			 double y_pos, 
			 double z_pos, 
			 double x_dest, 
			 double y_dest, 
			 double z_dest, 
			 int t_start, 
			 int mesh_start, 
			 double *p_in,
			 double *n_x,
			 double *n_y,
			 int *n_p,
			 int *neighbors,
			 int *forbidden,
			 unsigned *cell_type,
			 double *beta_v,
			 unsigned numberOfPoints,
			 unsigned numberOfTriangles,
			 float thicknessOfPrism,
			 float sigmaA,
			 float sigmaE,
			 int cladNumber,
			 float cladAbsorption,
			 float nTot
			 ){
  //    in first try no reflections
  //    calculate the vector and make the the calculation, which surface would be the shortest to reach
  //    then get the length, make the integration, get the information about the next cell out of the array
  //    set the point to the surface (this surface is "forbidden" in the calculations)
  //    proceed until you hit a the point or the surface
  //    if you are closer then "small" stop and return the value
  double vec_x, vec_y,vec_z, norm;
  double distance, length, length_help, distance_total;
  double gain=1;
  double nominator, denominator;
  int tri, cell_z; // the current triangle number and position concerning the z's
  int decider; // which one is the shortest - info
  int tri_next, cell_z_next, forb, forb_dump;
  int ct; // used to read out cell_type 
  unsigned size_p = numberOfPoints;
  unsigned N_cells = numberOfTriangles;
  float z_mesh = thicknessOfPrism;
  float sigma_a = sigmaA;
  float sigma_e = sigmaE;
  int clad_num = cladNumber;
  float clad_abs = cladAbsorption;
  float N_tot = nTot;

	
  //    initial positions
  tri = t_start;
  cell_z = mesh_start;
    
    
  //    definition of the vectors without reflections
  vec_x = (x_dest - x_pos);
  vec_y = (y_dest - y_pos);
  vec_z = (z_dest - z_pos);
    
  norm = sqrt(vec_x*vec_x+vec_y*vec_y+vec_z*vec_z);
    
  vec_x = vec_x/norm;
  vec_y = vec_y/norm;
  vec_z = vec_z/norm;
    
  //    now calculate the length to travel
  distance = sqrt((x_dest - x_pos)*(x_dest - x_pos)+(y_dest - y_pos)*(y_dest - y_pos)+(z_dest - z_pos)*(z_dest - z_pos));
  distance_total = distance;
  // does this make sense?
  length = distance;
    
  forb = -1;
	
  //	mexPrintf("Propagation called");
  //    mexEvalString("drawnow;");
    
  //    the ray has to be set to be ALIVE before!
  //    now do the unlimited for loop - break!!!
  for(;;)
    {
	
      //	  mexPrintf("Propagation for part called\n\n");
      //    mexEvalString("drawnow;");
      //        definition for decider
      //        0,1,2: int for the neighbors
      //        3: hor plane up
      //        4: hor plane down
        
      //        at first set the decider = -1;
      decider = -1;
      length = distance;
		
		
      //		  read, which type of cell it is you are propagation in
      ct = cell_type[tri];
        
      //        mexPrintf("forb: %i\n",forb);
      //        mexEvalString("drawnow;");
        
      //        try the triangle faces
      //        remember the correlation between the normals and the points
      //        n1: p1-2, n2: p1-3, n3:p2-3
      //        the third coordinate (z) of the particpating points for the surfaces can be set to be z=0, 
      //        as everything uses triangular "tubes/prisms", as well as n_z=0 in this case!
      if (forb != 0){
	nominator = (n_x[tri]*p_in[n_p[tri]] + n_y[tri]*p_in[n_p[tri]+size_p]) - (n_x[tri]*x_pos + n_y[tri]*y_pos);
	denominator = n_x[tri]*vec_x + n_y[tri]*vec_y;
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 0;
		forb_dump = (forbidden[tri]);
	      }
	  }
      }
        
      if (forb != 1){
	nominator = (n_x[tri+N_cells]*p_in[n_p[tri+N_cells]] + n_y[tri+N_cells]*p_in[n_p[tri+N_cells]+size_p]) - (n_x[tri+N_cells]*x_pos + n_y[tri+N_cells]*y_pos);
	denominator = n_x[tri+N_cells]*vec_x + n_y[tri+N_cells]*vec_y;
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 1;
		forb_dump = (forbidden[tri+N_cells]);
	      }
	  }
      }
        
      if (forb !=2){
	nominator = (n_x[tri+2*N_cells]*p_in[n_p[tri+2*N_cells]] + n_y[tri+2*N_cells]*p_in[n_p[tri+2*N_cells]+size_p]) - (n_x[tri+2*N_cells]*x_pos + n_y[tri+2*N_cells]*y_pos);
	denominator = n_x[tri+2*N_cells]*vec_x + n_y[tri+2*N_cells]*vec_y;
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 2;
		forb_dump = (forbidden[tri+2*N_cells]);
	      }
	  }
      }
        
      //        try the horizontal planes, which one is the shortest, n_x and n_y are zero!, n_z =1!
      //        at first the upper plane
      if (forb != 3){
	nominator = (cell_z+1)*z_mesh - z_pos;
	denominator = z_pos*vec_z;
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 3;
		forb_dump = 4; // you are not allowed to go down in the next step
	      }
	  }
      }
        
      //        next is the lower plane
      if (forb != 4){
	nominator = (cell_z)*z_mesh - z_pos;
	denominator = z_pos*vec_z;
            
	if (denominator != 0.0)
	  {
	    length_help = nominator/denominator;
	    if (length_help < length && length_help > VERY_SMALL)
	      {
		length = length_help;
		decider = 4;
		forb_dump = 3; // you are not allowed to go up in the next step
	      }
	  }
      }
        
      forb = forb_dump;
		
        
      //        now make a switch to differ the different cases
      switch(decider){
                
      case 0:
	//                this is the case for the intersection with the first choice triangle-surface
	tri_next = neighbors[tri];
	cell_z_next = cell_z;
	break;
                
      case 1:
	//                second triangle surface
	tri_next = neighbors[tri+N_cells];
	cell_z_next = cell_z;
	break;
                
      case 2:
	//                third triangle surface
	tri_next = neighbors[tri+2*N_cells];
	cell_z_next = cell_z;
	break;
                
      case 3:
	//                go one plane up
	tri_next = tri;
	cell_z_next = cell_z + 1;
	break;
                
      case 4:
	//                go one plane down
	tri_next = tri;
	cell_z_next = cell_z - 1;
	break;
                
      default:
	//                make an error statement
	break;
      }
        
      //        now we know where to go, let's make the integration
      //        take the beta_v[tri+cell_z*N_cells] 

      //		  at this position do the decision whether it is a gain part or cladding
      //		  it might be absorbing or amplifying, for the cladding only absorbing
      //		  a simple "if then"

      if (ct == clad_num){
	gain = gain * exp(-clad_abs * length);
      }
      else {
	gain = gain * exp(N_tot*(beta_v[tri+cell_z*N_cells]*(sigma_e + sigma_a)-sigma_a)*length);
      }
      //        gain = LineIntegralMCRK4_S(3, tri, cell_z, gain, length);
        
      //        after integration make the propagation
        
      //        mexPrintf("Distance: %f, Length: %f\n",distance, length);
      //        mexPrintf("decider: %i, forbidden: %i\n",decider, forb);
      //        mexPrintf("vec_x: %f, vec_y: %f, vec_z: %f\n", vec_x, vec_y, vec_z);
      //        mexPrintf("current_x: %f current_y: %f current_z: %f\n", x_pos, y_pos, z_pos);
      //        mexPrintf("tri: %i, tri_next: %i, cell_z: %i, cell_next: %i\n", tri, tri_next, cell_z, cell_z_next);
      //        mexEvalString("drawnow;");
      //        str=mxCreateString("Press a key");
      //        mexCallMATLAB(1,&dump,1,&str,"input"); 
      //        str and dump should be defined to be a *mxArray and don't forget to kill them at the end
        
      distance -= length;
        
      //        return 1;
      //        
        
      x_pos = x_pos + length*vec_x;
      y_pos = y_pos + length*vec_y;
      z_pos = z_pos + length*vec_z;
        
      if (abs(distance)< SMALL)
        {
	  break;
        }
        
        
      //        now set the next cell
      tri = tri_next;
      cell_z = cell_z_next;      
        
      //        break;
      //        now we should make the integration routine
    }
    
  gain /= (distance_total*distance_total);

  return gain;
}



/**
 * calculate the gain from the centers of each of the boxes to the observed point
 * calculate the gain and make a "mapping"
 * receipt: pick the point in the center of one cell, 
 * calculate the gain from this point to the observed point,
 * estimate the inner part of the Phi_ASE - Integral,
 * scale the amount of rays proportionally with it
 * sum the amount of rays and scale it to Int=1, which gives the inverse weights
 * the number of rays is determined via floor(), with ceil(), zero-redions could be added
 * use the routine "propagation"!, test: no reflections, just exponential
 *
 **/
void importf(int point,
	     int startLevel,
	     double *importance,
	     int *numberOfImportantRays,
	     double *points,
	     double *xOfNormals,
	     double *yOfNormals,
	     int *positionsOfNormalVectors,
	     int *neighbors,
	     int *forbidden,
	     unsigned *cellTypes,
	     double *betaValues,
	     double *xOfTriangleCenter,
	     double *yOfTriangleCenter,
	     float *surface,
	     unsigned raysPerSample,
	     unsigned numberOfPoints,
	     unsigned numberOfLevels,
	     unsigned numberOfTriangles,
	     float thicknessOfPrism,
	     float sigmaA,
	     float sigmaE,
	     int cladNumber,
	     float cladAbsorption,
	     float nTot
	     )
{
  int raysLeft;
  int raysDump;
  double sumPhi;
  double surfaceTotal;
  double xPos, yPos, zPos;
  double prop;

  raysDump = 0;
  sumPhi = 0.0;
  surfaceTotal = 0.0;
  xPos = points[point];
  yPos = points[point + numberOfPoints];
  zPos = startLevel * thicknessOfPrism;

  // Calculate importance by propagation from trianglecenter to every other center
  for (int i_t=0; i_t < numberOfTriangles; ++i_t){
    for (int i_z=0; i_z < (numberOfLevels-1); ++i_z){
      prop = rayPropagationCpu(xOfTriangleCenter[i_t], yOfTriangleCenter[i_t], 
      			       thicknessOfPrism * (i_z+0.5),  xPos, yPos, zPos, i_t, i_z, 
      			       points, xOfNormals, yOfNormals, positionsOfNormalVectors, 
      			       neighbors, forbidden , cellTypes, betaValues,
      			       numberOfPoints, numberOfTriangles, thicknessOfPrism,
      			       sigmaA, sigmaE, cladNumber, cladAbsorption, nTot
      			       );

      importance[i_t + i_z * numberOfTriangles] = betaValues[i_t + i_z * numberOfTriangles]*(prop);
      sumPhi += importance[i_t + i_z * numberOfTriangles];

    }
    surfaceTotal += surface[i_t];

  }

  // Calculate number of rays/prism
  for (int i_t=0; i_t < numberOfTriangles; ++i_t){
    for (int i_z=0; i_z < (numberOfLevels-1); ++i_z){
      numberOfImportantRays[i_t + i_z*numberOfTriangles] = (int)(floor(importance[i_t + i_z * numberOfTriangles] / sumPhi * raysPerSample));
      raysDump +=  numberOfImportantRays[i_t + i_z*numberOfTriangles];
      //fprintf(stderr, "[%d][%d] i: %.20f n: %d\n", i_z, i_t, importance[i_t + i_z*numberOfTriangles], numberOfImportantRays[i_t + i_z*numberOfTriangles]);
    }

  }
  raysLeft = raysPerSample - raysDump;

  // Distribute the remaining rays randomly
  for (int i_r=0; i_r < raysLeft; i_r++){
    int rand_t = (int )(rand() % numberOfTriangles);
    int rand_z = (int )(rand() % (numberOfLevels-1));
    numberOfImportantRays[rand_t + rand_z * numberOfTriangles]++;

  }

  //  Now think about the mount of rays which would come out of this volume(surface)
  //  dividing this number with the new amount of rays gives the final importance weight for this area!
  for (int i_t=0; i_t<numberOfTriangles; ++i_t){
    for (int i_z=0; i_z<(numberOfLevels-1); ++i_z){
      if (numberOfImportantRays[i_t + i_z*numberOfTriangles] > 0){
  	importance[i_t + i_z*numberOfTriangles] = raysPerSample * surface[i_t] / surfaceTotal / numberOfImportantRays[i_t + i_z*numberOfTriangles];

      }
      else{
  	importance[i_t + i_z*numberOfTriangles] = 0; 

      }


    }

  }


}


double* doubleVectorToArray(std::vector<double> *input){
  double* output;
  output = (double*) malloc(sizeof(double) * input->size());
  for(int i=0; i< input->size(); ++i){
    output[i] = input->at(i);	
  }
  return output;
}
int* intVectorToArray(std::vector<int> *input){
  int* output;
  output = (int*) malloc(sizeof(int) * input->size());
  for(int i=0; i< input->size(); ++i){
    output[i] = input->at(i);	
  }
  return output;
}
unsigned* unsignedVectorToArray(std::vector<unsigned> *input){
  unsigned* output;
  output = (unsigned*) malloc(sizeof(unsigned) * input->size());
  for(int i=0; i< input->size(); ++i){
    output[i] = input->at(i);	
  }
  return output;
}

//----------------------------------------------------
// Host Code
//----------------------------------------------------
/** GPU Kernel Variables
 * The idea is, that the number of threads is fixed (to maximize GPU occupancy)
 * and the number of blocks as well (200 is the maximum for the standard
 * Mersenne Twister implementaion). Therefore, the number of rays per sample
 * are fixed to be k*200*256.
 * That means, sometimes we have to increase the number of rays a little.
 *
 * \var raysPerThread is used to give every thread k iterations (to simulate k rays)
 *
 * note that every samplepoint receives the exact same number of rays.
 */
	
//

/** Variables for the device
 * These are on-GPU representations of the input parameters
 * of variable size.
 *
 * \var p_in: coordinates of the sample-points of one layer (first all x-coordinates, then all y-coordinates)
 * \var n_*: values of the normal-vectors for the 3 rectangular sides of each prism (described in 2D)
 * \var beta_v: the beta values of the prisms
 * \var phi: the accumulated ASE-Flux for each sample point
 * \var forbidden: the side of the prism through which the ray "entered" the prism
 * \var n_p: the points where the normals (n_x,n_y) start
 * \var neighbors: indices to the adjacent triangles in t_in
 * \var t_in: indices of the points which are considered to be a triangle (A points start from 0, B points from size_t, C points from size_t*2)
 * \var cell_type: determines which cell type we are looking at.
 * other input parameters are put to the GPU by the setupGlobalVariablesKernel
 */

float runRayPropagationGpu(
			   std::vector<double> *dndtAse, 
			   unsigned &threads, 
			   unsigned &blocks, 
			   unsigned &hostRaysPerSample,
			   std::vector<double> *betaValuesVector,
			   std::vector<double> *xOfNormalsVector,
			   std::vector<double> *yOfNormalsVector,
			   std::vector<unsigned> *cellTypesVector,
			   std::vector<unsigned> *triangleIndicesVector,
			   std::vector<int> *forbiddenVector,
			   std::vector<int> *neighborsVector,
			   std::vector<int> *positionsOfNormalVectorsVector,
			   std::vector<double> *pointsVector,
			   std::vector<double> *betaCellsVector,
			   std::vector<float> *surfacesVector,
			   std::vector<double> *xOfTriangleCenterVector,
			   std::vector<double> *yOfTriangleCenterVector,
			   float hostCladAbsorption,
			   int hostCladNumber,
			   float hostNTot,
			   float hostSigmaA,
			   float hostSigmaE,
			   unsigned hostNumberOfPoints,
			   unsigned hostNumberOfTriangles,
			   unsigned hostNumberOfLevels,
			   float hostThicknessOfPrism,
			   float hostCrystalFluorescence)
{
  fprintf(stderr, "\nConverting Vectors to Arrays\n");

  // Variable declarations
  // CPU
  // INPUT
  double* hostBetaValues;
  double* hostXOfNormals;
  double* hostYOfNormals;
  double* hostXOfTriangleCenter;
  double* hostYOfTriangleCenter;
  unsigned* hostCellTypes;
  unsigned* hostTriangleIndices;
  int* hostForbidden;
  int* hostNeighbors;
  int* hostPositionsOfNormalVectors;
  double* hostPoints;
  float* hostSurfaces;
  // TMP CALC
  unsigned hostNumberOfPrisms;
  unsigned hostRaysPerThread;
  cudaEvent_t start, stop;
  float runtimeGpu;
  float hostPhiASE[hostNumberOfPoints * (hostNumberOfLevels)];
  unsigned kernelcount;
  float surfaceTotal;
  int hostSurfacesNormalized[hostNumberOfTriangles];
  float minSurface;
  int surfacesNormalizedSum;
  // GPU
  double  *points, *xOfNormals, *yOfNormals, *betaValues;
  float *phiASE;
  int *forbidden, *positionsOfNormalVectors, *neighbors, *triangleIndices, *cellTypes, *surfacesNormalized;
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  int *numberOfImportantRays;
  double *importance;
	
  // Variables defintions
  threads = 120;
  blocks = 90;  // TODO: change number of blocks/threads dynamically to allow more flexible number of rays (increase up to 200)
  hostBetaValues = doubleVectorToArray(betaValuesVector);
  hostXOfNormals = doubleVectorToArray(xOfNormalsVector);
  hostYOfNormals = doubleVectorToArray(yOfNormalsVector);
  hostXOfTriangleCenter = doubleVectorToArray(xOfTriangleCenterVector);
  hostYOfTriangleCenter = doubleVectorToArray(yOfTriangleCenterVector);
  hostCellTypes = unsignedVectorToArray(cellTypesVector);
  hostTriangleIndices = unsignedVectorToArray(triangleIndicesVector);
  hostForbidden = intVectorToArray(forbiddenVector);
  hostNeighbors = intVectorToArray(neighborsVector);
  hostSurfaces =  (float*) &(surfacesVector[0]);
  hostPositionsOfNormalVectors = intVectorToArray(positionsOfNormalVectorsVector);
  hostPoints = doubleVectorToArray(pointsVector);
  hostNumberOfPrisms = (hostNumberOfTriangles * (hostNumberOfLevels-1));

  hostRaysPerThread = ceil(double(hostRaysPerSample) /  (blocks * threads));
  hostRaysPerSample = threads * blocks * hostRaysPerThread;
  
  runtimeGpu = 0.0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  kernelcount = 0;
  surfaceTotal=0;
  minSurface=9999999;
  surfacesNormalizedSum=0;
  double hostImportance[hostNumberOfPrisms];
  int hostNumberOfImportantRays[hostNumberOfPrisms];

  for(int i=0; i < hostNumberOfPoints * hostNumberOfLevels; ++i){
    hostPhiASE[i] = 0.0;

  }

  for(int i=0; i < hostNumberOfPrisms; ++i){
    hostNumberOfImportantRays[i] = 1;
    hostImportance[i] = 1.0;
  }

  //TODO: remove
  surfacesVector->pop_back();
  for(int i=0;i<surfacesVector->size();++i){
    surfaceTotal += surfacesVector->at(i);
    minSurface = min(minSurface,surfacesVector->at(i));
  }

  for(int i=0;i<surfacesVector->size();++i){
    hostSurfacesNormalized[i] = (surfacesVector->at(i)) / minSurface;
    surfacesNormalizedSum += hostSurfacesNormalized[i];
  }

  // Init mersenne twister PRNG
  {
    /**Allocate space for PRNG states on device */
    CUDA_CALL(cudaMalloc((void **)&devMTGPStates, blocks * sizeof(curandStateMtgp32)));

    /** Allocate space for MTGP kernel parameters */
    CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));

    /**Reformat from predefined parameter sets to kernel format,
     * and copy kernel parameters to device memory */
    CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));

    /** Initialize one state per thread block */
    /** \TODO initialize with time */
    CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, blocks, 1234));
  }



  // Allocation of memory on the GPU and setting of global GPU-variables
  {
    fprintf(stderr, "\nFilling the device Variables\n");
    //Create constant values on GPU
    setupGlobalVariablesKernel<<<1,1>>>(
					double(hostSigmaE), 
					double(hostSigmaA),
					hostCladNumber,
					double(hostCladAbsorption),
					double(hostNTot), 
					hostNumberOfTriangles, 
					double(hostThicknessOfPrism),
					hostNumberOfLevels, 
					hostNumberOfPoints); //@OPTIMIZE: initialize the constants as constants...

    cudaThreadSynchronize();

    // Memory allocation on device
    CUDA_CHECK_RETURN(cudaMalloc(&points, 2 * hostNumberOfPoints * sizeof(double)));
    CUDA_CHECK_RETURN(cudaMalloc(&xOfNormals, 3 * hostNumberOfTriangles * sizeof(double)));
    CUDA_CHECK_RETURN(cudaMalloc(&yOfNormals, 3 * hostNumberOfTriangles * sizeof(double)));
    CUDA_CHECK_RETURN(cudaMalloc(&neighbors, 3* hostNumberOfTriangles * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&forbidden, 3* hostNumberOfTriangles * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&positionsOfNormalVectors, 3* hostNumberOfTriangles * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&triangleIndices, 3* hostNumberOfTriangles * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&cellTypes,hostNumberOfTriangles * hostNumberOfLevels * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&betaValues,hostNumberOfTriangles * (hostNumberOfLevels-1) * sizeof(double)));
    CUDA_CHECK_RETURN(cudaMalloc(&phiASE,hostNumberOfPoints * hostNumberOfLevels * sizeof(float)));
    CUDA_CHECK_RETURN(cudaMalloc(&surfacesNormalized,hostNumberOfTriangles * sizeof(int)));
    // Memory importance sampling
    CUDA_CHECK_RETURN(cudaMalloc(&importance, hostNumberOfPrisms * sizeof(double)));
    CUDA_CHECK_RETURN(cudaMalloc(&numberOfImportantRays, hostNumberOfPrisms * sizeof(int)));

    /// Copy data from host to device
    CUDA_CHECK_RETURN(cudaMemcpy(points, hostPoints, 2 * hostNumberOfPoints * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(xOfNormals, hostXOfNormals, 3 * hostNumberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(yOfNormals, hostYOfNormals, 3 * hostNumberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(neighbors, hostNeighbors, 3 * hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(forbidden,hostForbidden, 3 * hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(positionsOfNormalVectors ,hostPositionsOfNormalVectors, 3 * hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(triangleIndices ,hostTriangleIndices, 3* hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(cellTypes,hostCellTypes, hostNumberOfTriangles * hostNumberOfLevels * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(betaValues, hostBetaValues, hostNumberOfTriangles * (hostNumberOfLevels-1) * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(phiASE, hostPhiASE, hostNumberOfPoints * hostNumberOfLevels * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(surfacesNormalized, hostSurfacesNormalized, hostNumberOfTriangles * sizeof(int),cudaMemcpyHostToDevice));
    // Copy importance sampling data
    CUDA_CHECK_RETURN(cudaMemcpy(importance, hostImportance, hostNumberOfPrisms * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(numberOfImportantRays, hostNumberOfImportantRays, hostNumberOfPrisms * sizeof(int), cudaMemcpyHostToDevice));
  }

  fprintf(stderr, "hostCrystalFluorescence: %e\n",hostCrystalFluorescence);
  //testKernel<<<1,1>>>(points, xOfNormals, yOfNormals, positionsOfNormalVectors, neighbors, forbidden, triangleIndices, cellTypes, betaValues,surfacesNormalized);

  // Start Kernels
  {
    fprintf(stderr, "\nStarting the propagation\n");
    cudaEventRecord(start, 0);
		
    // Every Kernel calculates one sample point
    for(int point2D = 0; point2D < hostNumberOfPoints ; ++point2D){
      //for(int level = 0; level < hostNumberOfLevels; ++level){
       for(int level = 0; level < 1; ++level){
	cudaThreadSynchronize();
	// Importance for one sample
	importf(point2D, level, hostImportance, hostNumberOfImportantRays, 
		hostPoints, hostXOfNormals, hostYOfNormals, hostPositionsOfNormalVectors, 
		hostNeighbors, hostForbidden, hostCellTypes, hostBetaValues, 
		hostXOfTriangleCenter, hostYOfTriangleCenter, hostSurfaces, hostRaysPerSample,
		hostNumberOfPoints, hostNumberOfLevels, hostNumberOfTriangles, hostThicknessOfPrism,
		hostSigmaA, hostSigmaE, hostCladNumber, hostCladAbsorption,hostNTot
		);
	CUDA_CHECK_RETURN(cudaMemcpy(importance, hostImportance, hostNumberOfPrisms * sizeof(double), cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(cudaMemcpy(numberOfImportantRays, hostNumberOfImportantRays, hostNumberOfPrisms * sizeof(int), cudaMemcpyHostToDevice));

	// Calculate for one sample
	raytraceStep<<< blocks, threads >>> ( devMTGPStates, phiASE, point2D, level, hostRaysPerThread, 
					      points, xOfNormals, yOfNormals, positionsOfNormalVectors, 
					      neighbors, forbidden, triangleIndices, cellTypes, betaValues, importance, 
					      numberOfImportantRays, surfacesNormalized );

	if(kernelcount % 200 == 0)
	  fprintf(stderr, "Sampling point %d done\n",kernelcount);
	kernelcount++;
      }
    }

    cudaThreadSynchronize();
  }


  // Final calculations
  {
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&runtimeGpu, start, stop);
    CUDA_CHECK_RETURN(cudaMemcpy(hostPhiASE, phiASE, hostNumberOfPoints * hostNumberOfLevels * sizeof(float), cudaMemcpyDeviceToHost));
    //int raysPerSampleNormalized = hostThreadsPerPrism * surfacesNormalizedSum * (hostNumberOfLevels-1) * hostRaysPerThread;
    for(int i=0; i< hostNumberOfPoints;++i){
      for(int j=0 ; j<hostNumberOfLevels ; ++j)
	{
	  int pos = i*hostNumberOfLevels+j;
	  hostPhiASE[pos] = float( (double(hostPhiASE[pos]) / (hostRaysPerSample * 4.0f * 3.14159)));
	  double gain_local = double(hostNTot)*(betaCellsVector->at(pos))*double(hostSigmaE+hostSigmaA)-double(hostNTot*hostSigmaA);
	  dndtAse->at(pos) = gain_local*hostPhiASE[pos]/hostCrystalFluorescence;
	}
    }
  }

  // Free Memory
  {
    cudaFree(points);
    cudaFree(xOfNormals);
    cudaFree(yOfNormals);
    cudaFree(neighbors);
    cudaFree(forbidden);
    cudaFree(positionsOfNormalVectors);
    cudaFree(betaValues);
    cudaFreeHost(hostBetaValues);
    cudaFreeHost(hostXOfNormals);
    cudaFreeHost(hostYOfNormals);
    cudaFreeHost(hostCellTypes);
    cudaFreeHost(hostTriangleIndices);
    cudaFreeHost(hostForbidden);
    cudaFreeHost(hostNeighbors);
    cudaFreeHost(hostPositionsOfNormalVectors);
    cudaFreeHost(hostPoints);
    cudaFree(importance);
    cudaFree(numberOfImportantRays);
  }

  cudaDeviceReset();
  return runtimeGpu;
}

