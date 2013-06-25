#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector_types.h>
#include <assert.h>
#include <vector>
#include <curand_kernel.h>
#include <cudachecks.h>
#include <importance_sampling.h>
#include <test_functions.h>
#include <calc_sample_phi_ase.h>
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>

#define TEST_VALUES true
#define SMALL 1E-06
#define VERY_SMALL 0.0


/*********************************************************************************************
 * HOST FUNCTIONS
 *********************************************************************************************/
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
float calcDndtAse(
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
		  unsigned hostCladNumber,
		  float hostNTot,
		  float hostSigmaA,
		  float hostSigmaE,
		  unsigned hostNumberOfPoints,
		  unsigned hostNumberOfTriangles,
		  unsigned hostNumberOfLevels,
		  float hostThicknessOfPrism,
		  float hostCrystalFluorescence)
{
  // Variable declarations
  // CPU
  double* hostImportance;
  unsigned* hostNumberOfImportantRays;
  int* hostIndicesOfPrisms;
  unsigned hostNumberOfPrisms;
  unsigned hostRaysPerThread;
  unsigned hostNumberOfSamples;
  cudaEvent_t start, stop;
  float runtimeGpu;
  float *hostPhiASE;
  unsigned kernelcount;
  // GPU
  double  *points, *xOfNormals, *yOfNormals, *betaValues;
  float *phiASE;
  int *forbidden, *positionsOfNormalVectors, *neighbors, *triangleIndices, *cellTypes, *surfacesNormalized;
  curandStateMtgp32 *devMTGPStates;
  mtgp32_kernel_params *devKernelParams;
  double *importance;
  unsigned *indicesOfPrisms;
	
  // Variables defintions
  threads = 50; //OPTIMIZE: find perfect number of threads
  blocks = 200;
  hostNumberOfPrisms = (hostNumberOfTriangles * (hostNumberOfLevels-1));
  hostNumberOfSamples = hostNumberOfPoints * hostNumberOfLevels;

  
  hostPhiASE = (float*) malloc(hostNumberOfSamples * sizeof(float));
  hostImportance = (double*) malloc(hostNumberOfPrisms * sizeof(double));
  hostNumberOfImportantRays = (unsigned*) malloc(hostNumberOfPrisms * sizeof(unsigned));
  hostIndicesOfPrisms = (int*) malloc(hostRaysPerSample * sizeof(int));

  runtimeGpu = 0.0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  kernelcount = 0;

  for(int i=0; i < hostRaysPerSample; ++i){
    hostIndicesOfPrisms[i] = 0;

  }

  for(int i=0; i < hostNumberOfSamples; ++i){
    hostPhiASE[i] = 0.f;

  }

  for(int i=0; i < hostNumberOfPrisms; ++i){
    hostNumberOfImportantRays[i] = 1;
    hostImportance[i] = 1.0;
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
    //fprintf(stderr, "\nC Filling the device Variables\n");
    //Create constant values on GPU
    setupGlobalVariablesKernel<<<1,1>>>(double(hostSigmaE), 
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
    CUDA_CHECK_RETURN(cudaMalloc(&neighbors, 3 * hostNumberOfTriangles * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&forbidden, 3 * hostNumberOfTriangles * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&positionsOfNormalVectors, 3 * hostNumberOfTriangles * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&triangleIndices, 3 * hostNumberOfTriangles * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&cellTypes, hostNumberOfTriangles * sizeof(int)));
    CUDA_CHECK_RETURN(cudaMalloc(&betaValues, hostNumberOfPrisms * sizeof(double)));
    CUDA_CHECK_RETURN(cudaMalloc(&phiASE, hostNumberOfSamples * sizeof(float)));
    CUDA_CHECK_RETURN(cudaMalloc(&importance, hostNumberOfPrisms * sizeof(double)));
    CUDA_CHECK_RETURN(cudaMalloc(&indicesOfPrisms, hostRaysPerSample * sizeof(unsigned)));

    /// Copy data from host to device
    CUDA_CHECK_RETURN(cudaMemcpy(points, (double*) &(pointsVector->at(0)), 2 * hostNumberOfPoints * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(xOfNormals, (double*) &(xOfNormalsVector->at(0)), 3 * hostNumberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(yOfNormals, (double*) &(yOfNormalsVector->at(0)), 3 * hostNumberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(neighbors,(int*) &(neighborsVector->at(0)), 3 * hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(forbidden, (int*) &(forbiddenVector->at(0)), 3 * hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(positionsOfNormalVectors, (int*) &(positionsOfNormalVectorsVector->at(0)), 3 * hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(triangleIndices, (unsigned*) &(triangleIndicesVector->at(0)), 3 * hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(cellTypes, (unsigned*) &(cellTypesVector->at(0)), hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(betaValues, (double*) &(betaValuesVector->at(0)), hostNumberOfPrisms * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(phiASE, hostPhiASE, hostNumberOfSamples * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK_RETURN(cudaMemcpy(importance, hostImportance, hostNumberOfPrisms * sizeof(double), cudaMemcpyHostToDevice));
  }

  //fprintf(stderr, "C hostCrystalFluorescence: %e\n",hostCrystalFluorescence);
  // testKernel<<<1,1>>>(points, xOfNormals, yOfNormals, positionsOfNormalVectors, neighbors, forbidden, triangleIndices, cellTypes, betaValues,surface);

  // Start Kernels
  {
    fprintf(stderr, "C Start Kernel\n");
    cudaEventRecord(start, 0);
		
    // Every Kernel calculates one sample point
    for(int point_i = 0; point_i < hostNumberOfPoints ; ++point_i){
      for(int level_i = 0; level_i < hostNumberOfLevels; ++level_i){
	cudaThreadSynchronize();

	// Reinit (TODO is really needed ?)
	for(int i=0; i < hostNumberOfPrisms; ++i){
	  hostNumberOfImportantRays[i] = 1;
	  hostImportance[i] = 1.0;
	}
	for(int i=0; i < hostRaysPerSample; ++i){
	  hostIndicesOfPrisms[i] = 0;

	}

	// Importance for one sample
	importanceSampling(point_i, level_i, hostImportance, hostNumberOfImportantRays, 
			   (double*) &(pointsVector->at(0)), 
			   (double*) &(xOfNormalsVector->at(0)), 
			   (double*) &(yOfNormalsVector->at(0)),
			   (int*) &(positionsOfNormalVectorsVector->at(0)), 
			   (int*) &(neighborsVector->at(0)), 
			   (int*) &(forbiddenVector->at(0)), 
			   (unsigned*) &(cellTypesVector->at(0)), 
			   (double*) &(betaValuesVector->at(0)), 
			   (double*) &(xOfTriangleCenterVector->at(0)),
			   (double*) &(yOfTriangleCenterVector->at(0)), 
			   (float*) &(surfacesVector->at(0)), 
			   hostRaysPerSample,hostNumberOfPoints, hostNumberOfLevels, hostNumberOfTriangles, 
			   hostThicknessOfPrism, hostSigmaA, hostSigmaE, hostCladNumber, hostCladAbsorption,hostNTot);

	CUDA_CHECK_RETURN(cudaMemcpy(importance, hostImportance, hostNumberOfPrisms * sizeof(double), cudaMemcpyHostToDevice));

	// Prism scheduling for gpu threads
	int index=0;
	for(int prism_i=0; prism_i < hostNumberOfPrisms; ++prism_i){
		for(int ray_i=0 ; ray_i < hostNumberOfImportantRays[prism_i]; ++ray_i){
			hostIndicesOfPrisms[index] = prism_i;
			index++;
			assert(index <= hostRaysPerSample);
		}
	}


	CUDA_CHECK_RETURN(cudaMemcpy(indicesOfPrisms, hostIndicesOfPrisms, hostRaysPerSample * sizeof(unsigned), cudaMemcpyHostToDevice));

	cudaDeviceSynchronize();
	calcSamplePhiAse<<< blocks, threads >>> ( devMTGPStates, phiASE, point_i, level_i, hostRaysPerThread, 
					      points, xOfNormals, yOfNormals, positionsOfNormalVectors, 
					      neighbors, forbidden, triangleIndices, cellTypes, betaValues, importance, 
					      surfacesNormalized,indicesOfPrisms,hostRaysPerSample );


	if(kernelcount % 200 == 0)
	   fprintf(stderr, "C Sampling point %d done\n",kernelcount);
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

    for(int sample_i=0; sample_i < hostNumberOfSamples; ++sample_i){
      hostPhiASE[sample_i] = float( (double(hostPhiASE[sample_i]) / (hostRaysPerSample * 4.0f * 3.14159)));
      double gain_local = double(hostNTot) * (betaCellsVector->at(sample_i)) * double(hostSigmaE + hostSigmaA) - double(hostNTot * hostSigmaA);
      dndtAse->at(sample_i) = gain_local * hostPhiASE[sample_i] / hostCrystalFluorescence;
      	
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
    cudaFree(importance);
  }

  cudaDeviceReset();
  return runtimeGpu;
}

