// Libraries
#include <stdio.h> /* fprintf, memcpy, strstr, strcmp */
#include <assert.h> /* assert */
#include <string> /* string */
#include <vector> /* vector */

// User header files
#include <calc_dndt_ase.h>
#include <parser.h>
#include <write_to_vtk.h>

int main(int argc, char **argv){
  unsigned raysPerSample = 0;
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
  std::vector<unsigned> * triangleIndices = new std::vector<unsigned>;
  std::vector<int> * forbidden = new std::vector<int>;
  std::vector<int> * neighbors = new std::vector<int>;
  std::vector<int> * positionsOfNormalVectors = new std::vector<int>;
  std::vector<double> * points = new std::vector<double>;
  std::vector<double> * betaCells = new std::vector<double>;
  std::vector<float> * surfaces = new std::vector<float>;
  std::vector<double> *xOfTriangleCenter = new std::vector<double>;
  std::vector<double> *yOfTriangleCenter = new std::vector<double>;
  float nTot = 0;
  float sigmaA = 0;
  float sigmaE = 0;
  unsigned numberOfPoints = 0;
  unsigned numberOfTriangles = 0;
  unsigned numberOfLevels = 0;
  float thicknessOfPrism = 1;
  float crystalFluorescence = 0;

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
      raysPerSample = atoi(pos+1);
    }
  }

  // Parse location of experiements
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--experiment=", 12) == 0){
      memcpy (experimentLocation, argv[i]+13, strlen(argv[i])-13 );
    } 
  }

  // Check if we want no output
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--silent", 7) == 0){
		silent=true;
    } 
  }

  // Parse experimentdata
  std::string root(experimentLocation);
  
  // Add slash at the end, if missing
  if(root[root.size()-1] == 'w')
    root.erase(root.size()-1, 1);
  else if(root[root.size()-1] != '/')
    root.append("/");

  if(fileToVector(root + "n_p.txt", positionsOfNormalVectors)) return 1;
  if(fileToVector(root + "beta_v.txt", betaValues)) return 1;
  if(fileToVector(root + "forbidden.txt", forbidden)) return 1;
  if(fileToVector(root + "neighbors.txt", neighbors)) return 1;
  if(fileToVector(root + "n_x.txt", xOfNormals)) return 1;
  if(fileToVector(root + "n_y.txt", yOfNormals)) return 1;
  if(fileToVector(root + "x_center.txt", xOfTriangleCenter)) return 1;
  if(fileToVector(root + "y_center.txt", yOfTriangleCenter)) return 1;
  if(fileToVector(root + "p_in.txt", points)) return 1;
  if(fileToVector(root + "beta_cell.txt", betaCells)) return 1;
  if(fileToVector(root + "t_in.txt", triangleIndices)) return 1;
  if(fileToVector(root + "surface.txt", surfaces)) return 1;

  if(fileToValue(root + "n_tot.txt", nTot)) return 1;
  if(fileToValue(root + "sigma_a.txt", sigmaA)) return 1;
  if(fileToValue(root + "sigma_e.txt", sigmaE)) return 1;
  if(fileToValue(root + "size_p.txt", numberOfPoints)) return 1;
  if(fileToValue(root + "size_t.txt", numberOfTriangles)) return 1;
  if(fileToValue(root + "mesh_z.txt", numberOfLevels)) return 1;
  if(fileToValue(root + "z_mesh.txt", thicknessOfPrism)) return 1;
  if(fileToValue(root + "tfluo.txt", crystalFluorescence)) return 1;

  // Debug
  // fprintf(stderr, "C nTot: %e\n", nTot);
  // fprintf(stderr, "C sigmaA: %e\n", sigmaA);
  // fprintf(stderr, "C sigmaE: %e\n", sigmaE);
  // fprintf(stderr, "C numberOfPoints: %d\n", numberOfPoints);
  // fprintf(stderr, "C numberOfTriangles: %d\n", numberOfTriangles); 
  // fprintf(stderr, "C numberOfLevels: %d\n\n", numberOfLevels);

  // Test vectors
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
  assert((numberOfTriangles * (numberOfLevels-1)) <= raysPerSample);

  // Solution vector
  std::vector<double> *ase = new std::vector<double>(numberOfPoints * numberOfLevels, 0);

  // Run Experiment
  for(i=1; i < argc; ++i){
    if(strncmp(argv[i], "--mode=", 6) == 0){
      if(strstr(argv[i], "ray_propagation_gpu") != 0){
	// threads and blocks will be set in the following function (by reference)
	runtime = calcDndtAse(
			ase,
			threads, 
			blocks, 
			raysPerSample,
			betaValues,
			xOfNormals,
			yOfNormals,
			triangleIndices,
			forbidden,
			neighbors,
			positionsOfNormalVectors,
			points,
			betaCells,
			surfaces,
			xOfTriangleCenter,
			yOfTriangleCenter,
			nTot,
			sigmaA,
			sigmaE,
			numberOfPoints,
			numberOfTriangles,
			numberOfLevels,
			thicknessOfPrism,
			crystalFluorescence);
	strcpy(runmode, "Ray Propagation GPU");
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
  fprintf(stderr, "C Prism             : %d\n", (int) numberOfPoints * (numberOfLevels-1));
  fprintf(stderr, "C Samples           : %d\n", (int) ase->size());
  fprintf(stderr, "C Rays/Sample       : %d\n", raysPerSample);
  fprintf(stderr, "C Rays Total        : %d\n", raysPerSample * ase->size());
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


