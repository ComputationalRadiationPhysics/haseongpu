// Libraries
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include "string.h"
#include <vector>
#include "curand_kernel.h"


// User header files
#include "datatypes.h"
#include "geometry.h"
#include "datatypes.h"
#include "generate_testdata.h"
#include "print.h"
#include "geometry_gpu.h"
#include "ase_bruteforce_kernel.h"
#include "ase_bruteforce_cpu.h"
//#include "testdata_transposed.h"
#include "ray_propagation_gpu.h"
#include "buildgrid.h"
#include "parser.h"
#include "write_to_vtk.h"

int main(int argc, char **argv){
  unsigned raysTotal;
  char runmode[100];
  char experimentLocation[256];
  float runtime = 0.0;
  unsigned blocks = 0;
  unsigned threads = 0;
  bool silent = false;
  
  // Experimentdata
  std::vector<double> * betaValues = new std::vector<double>;
  std::vector<double> * xOfNormals = new std::vector<double>;
  std::vector<double> * yOfNormals = new std::vector<double>;
  std::vector<unsigned> * cellTypes = new std::vector<unsigned>;
  std::vector<unsigned> * triangleIndices = new std::vector<unsigned>;
  std::vector<int> * forbidden = new std::vector<int>;
  std::vector<int> * neighbors = new std::vector<int>;
  std::vector<int> * positionsOfNormalVectors = new std::vector<int>;
  std::vector<double> * points = new std::vector<double>;
  float cladAbsorption = 0;
  float cladNumber = 0;
  float nTot = 0;
  float sigmaA = 0;
  float sigmaE = 0;
  unsigned numberOfPoints = 0;
  unsigned numberOfTriangles = 0;
  unsigned numberOfLevels = 0;
  float thicknessOfPrism = 1;

  // Parse Commandline
  if(argc <= 1){
    fprintf(stderr, "C No commandline arguments found\n");
    fprintf(stderr, "C Usage    : ./octrace --mode=[runmode] --rays=[number of rays] --experiment=[location to .zip]\n");
    fprintf(stderr, "C Runmodes : bruteforce_gpu\n");
    fprintf(stderr, "             ray_propagation_gpu\n");
    return 0;
  }
  
  // Parse number of rays
  unsigned i;
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--rays=", 6) == 0){
      const char* pos = strrchr(argv[i],'=');
      raysTotal = atoi(pos+1);
    }
  }

  // Parse location of experiements
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--experiment=", 12) == 0){
      memcpy (experimentLocation, argv[i]+13, strlen(argv[i])-13 );
    } 
  }

  // check if we want no output
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--silent", 7) == 0){
		silent=true;
    } 
  }

  if(parse(experimentLocation, betaValues, xOfNormals, yOfNormals, cellTypes, triangleIndices, forbidden, neighbors, positionsOfNormalVectors, points, &cladAbsorption, &cladNumber, &nTot, &sigmaA, &sigmaE, &thicknessOfPrism, &numberOfPoints, &numberOfTriangles, &numberOfLevels)){
    fprintf(stderr, "C Had problems while parsing experiment data\n");
    return 1;
  }

  // Debug
  // fprintf(stderr, "cladAbsorption: %f\n", cladAbsorption);
  // fprintf(stderr, "cladNumber: %f\n", cladNumber);
  // fprintf(stderr, "nTot: %e\n", nTot);
  // fprintf(stderr, "sigmaA: %e\n", sigmaA);
  // fprintf(stderr, "sigmaE: %e\n", sigmaE);
  // fprintf(stderr, "numberOfPoints: %d\n", numberOfPoints);
  // fprintf(stderr, "numberOfTriangles: %d\n", numberOfTriangles); 
  // fprintf(stderr, "numberOfLevels: %d\n\n", numberOfLevels);

  // Generate testdata
  fprintf(stderr, "C Generate Testdata\n");
  std::vector<PrismCu>  *prisms = generatePrismsFromTestdata(numberOfLevels, points, numberOfPoints, triangleIndices, numberOfTriangles, thicknessOfPrism);
  std::vector<PointCu> *samples = generateSamplesFromTestdata(numberOfLevels, points, numberOfPoints);
  std::vector<double>      *ase = new std::vector<double>(samples->size(), 0);


  // Run Experiment
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--mode=", 6) == 0){
      if(strstr(argv[i], "bruteforce_gpu") != 0){
	// threads and blocks will be set in the following function (by reference)
  	runtime = runAseBruteforceGpu(samples, 
				      prisms, 
				      betaValues, 
				      ase,
				      cladAbsorption,
				      cladNumber,
				      nTot,
				      sigmaA,
				      sigmaE,
				      threads, 
				      blocks, 
				      raysTotal);
	strcpy(runmode, "Bruteforce GPU");
	break;
      }
      else if(strstr(argv[i], "ray_propagation_gpu") != 0){
	// threads and blocks will be set in the following function (by reference)
	runtime = runRayPropagationGpu(
			ase,
			threads, 
			blocks, 
			raysTotal,
			betaValues,
			xOfNormals,
			yOfNormals,
			cellTypes,
			triangleIndices,
			forbidden,
			neighbors,
			positionsOfNormalVectors,
			points,
			cladAbsorption,
			cladNumber,
			nTot,
			sigmaA,
			sigmaE,
			numberOfPoints,
			numberOfTriangles,
			numberOfLevels,
			thicknessOfPrism);
	strcpy(runmode, "Naive Ray Propagation GPU");
	break;
      }
    
    }

  }

  // Print Solution
  unsigned sample_i;
  fprintf(stderr, "C Solutions\n");
  for(sample_i = 0; sample_i < ase->size(); ++sample_i){
    fprintf(stderr, "C ASE PHI of sample %d: %.80f\n", sample_i, ase->at(sample_i));
	if(silent){
		if(sample_i >= 10) break;
	}
  }

  // Print statistics
  fprintf(stderr, "\n");
  fprintf(stderr, "C Statistics\n");
  fprintf(stderr, "C Prism             : %d\n", (int) prisms->size());
  fprintf(stderr, "C Samples           : %d\n", (int) samples->size());
  fprintf(stderr, "C Rays/Sample       : %d\n", raysTotal / samples->size());
  fprintf(stderr, "C Rays Total        : %d\n", raysTotal);
  fprintf(stderr, "C GPU Blocks        : %d\n", blocks);
  fprintf(stderr, "C GPU Threads/Block : %d\n", threads);
  fprintf(stderr, "C GPU Threads Total : %d\n", threads * blocks);
  fprintf(stderr, "C Runmode           : %s \n", runmode);
  fprintf(stderr, "C Runtime           : %f s\n", runtime / 1000.0);
  fprintf(stderr, "\n");

  // Write experiment data to vtk
  writeToVtk(points, numberOfPoints, triangleIndices, numberOfTriangles, numberOfLevels, thicknessOfPrism, ase);

  return 0;
}


