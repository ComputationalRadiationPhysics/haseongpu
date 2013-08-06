// Libraries
#include <stdio.h> /* fprintf, memcpy, strstr, strcmp */
#include <assert.h> /* assert */
#include <string> /* string */
#include <vector> /* vector */

// User header files
#include <calc_dndt_ase.h>
#include <parser.h>
#include <write_to_vtk.h>
#include <write_dndt_ase.h>
#include <for_loops_clad.h>
#include <cudachecks.h>
#include <mesh.h>

#define MIN_COMPUTE_CAPABILITY_MAJOR 2
#define MIN_COMPUTE_CAPABILITY_MINOR 0

/** 
 * @brief Queries for devices on the running mashine and collects
 *        them on the devices array. Set the first device in this 
 *        array as computaion-device. On Errors the programm will
 *        be stoped by exit(). Otherwise you can set the device by command
 *        line parameter --device=
 * 
 * @param verbose > 0 prints debug output
 *        devices Array of possible devices to use
 *        device  number of device you want to set
 * 
 * @return Number of devices in devices array
 */
unsigned getCorrectDevice(int verbose,unsigned **devices, int device){
  int count = 0, candidate = 0;
  unsigned correctDevices = 0;
  cudaDeviceProp prop;
  int minMajor = MIN_COMPUTE_CAPABILITY_MAJOR;
  int minMinor = MIN_COMPUTE_CAPABILITY_MINOR;

  CUDA_CHECK_RETURN( cudaGetDeviceCount(&count));
  
  for(int i=0; i<count; ++i){
	  CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, i) );
	  if( (prop.major > minMajor) || (prop.major == minMajor && prop.minor >= minMinor) ){
		  correctDevices++;
	  }
  }

  if(correctDevices == 0){
    fprintf(stderr,"\nNone of the CUDA-capable devices is sufficient!\n");
    exit(1);
  }

  (*devices) = (unsigned*) malloc(sizeof(unsigned) * correctDevices);

  if(verbose > 0){
	  fprintf(stderr,"\nFound %d CUDA devices with Compute Capability >= %d.%d):\n", correctDevices, minMajor,minMinor); 
  }

  candidate = 0;
  for(int i=0; i<count; ++i){
    CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, i) );
    if( (prop.major > minMajor) || (prop.major == minMajor && prop.minor >= minMinor) ){
		if(verbose > 0){
			fprintf(stderr,"[%d] %s (Compute Capability %d.%d)\n", candidate, prop.name, prop.major, prop.minor); 
		}
		(*devices)[candidate]=i;
		candidate++;
    }
  }

  if(device == -1){
    CUDA_CHECK_RETURN( cudaSetDevice((*devices)[0]) );
  }
  else{
    CUDA_CHECK_RETURN( cudaSetDevice((*devices)[device]));
  }
  return correctDevices;
}

int main(int argc, char **argv){
  unsigned raysPerSample = 0;
  char runmode[100];
  char experimentLocation[256] = "";
  char compareLocation[256] = "";
  float runtime = 0.0;
  unsigned blocks;
  unsigned threads;
  bool silent = false;
  unsigned *devices; // will be assigned in getCOrrectDevice();
  unsigned numberOfDevices=0;
  int device = -1;
  
  // Constant data
  float nTot = 0;
  float crystalFluorescence = 0;
  std::vector<double> * betaCells = new std::vector<double>;
  std::vector<double> *sigmaA = new std::vector<double>;
  std::vector<double> *sigmaE = new std::vector<double>;

  // Parse Commandline
  if(argc <= 1){
    fprintf(stderr, "C No commandline arguments found\n");
    fprintf(stderr, "C Usage    : ./octrace --mode=[runmode] --rays=[number of rays] --experiment=[location to experiment-data]\n");
    fprintf(stderr, "C Runmodes : for_loops\n");
    fprintf(stderr, "             ray_propagation_gpu\n");
    return 1;
  }
  
  // Parse number of rays
  for(int i=1; i < argc; ++i){
    if(strncmp(argv[i], "--rays=", 6) == 0){
      const char* pos = strrchr(argv[i],'=');
      raysPerSample = atoi(pos+1);
      if(raysPerSample == 0){
	fprintf(stderr, "C Please specify the number of rays per sample Point with --rays=\n");
      }
    }
  }

  // Parse location of experiements
  for(int i=1; i < argc; ++i){
    if(strncmp(argv[i], "--experiment=", 12) == 0){
      memcpy (experimentLocation, argv[i]+13, strlen(argv[i])-13 );
    } 
  }

  // Parse which cuda device to choose
  for(int i=1; i < argc; ++i){
    if(strncmp(argv[i], "--device=", 8) == 0){
      const char* pos = strrchr(argv[i],'=');
      device = atoi(pos+1);
    } 
  }

  // Parse what vtk file to compare with
  for(int i=1; i < argc; ++i){
    if(strncmp(argv[i], "--compare=", 9) == 0){
      memcpy (compareLocation, argv[i]+10, strlen(argv[i])-10 );
    } 
  }


  // Check if we want no output
  for(int i=1; i < argc; ++i){
    if(strncmp(argv[i], "--silent", 7) == 0){
		silent=true;
    } 
  }
  
  std::string root(experimentLocation);

  // Add slash at the end, if missing
  if(root[root.size()-1] == 'w')
    root.erase(root.size()-1, 1);
  else if(root[root.size()-1] != '/')
    root.append("/");

  // Parse constant from files
  if(fileToValue(root + "n_tot.txt", nTot)) return 1;
  if(fileToValue(root + "tfluo.txt", crystalFluorescence)) return 1;
  if(fileToVector(root + "beta_cell.txt", betaCells)) return 1;
  if(fileToVector(root + "sigma_a.txt", sigmaA)) return 1;
  if(fileToVector(root + "sigma_e.txt", sigmaE)) return 1;
  assert(sigmaA->size() == sigmaE->size());
  
  // Set/Test device to run experiment
  numberOfDevices = getCorrectDevice(1,&devices, device);

  // Parse experiemntdata and fill mesh 
  Mesh hMesh;
  Mesh *dMesh = new Mesh[numberOfDevices];
  if(Mesh::parseMultiGPU(&hMesh, &dMesh, root, numberOfDevices, devices)) return 1;

  // Debug
  // fprintf(stderr, "C nTot: %e\n", nTot);
  // fprintf(stderr, "C sigmaA: %e\n", sigmaA);
  // fprintf(stderr, "C sigmaE: %e\n", sigmaE);
  // fprintf(stderr, "C numberOfTriangles: %d\n", hMesh.numberOfTriangles);
  // fprintf(stderr, "C numberOfLevels: %d\n", hMesh.numberOfLevels); 
  // fprintf(stderr, "C numberOfPrisms: %d\n", hMesh.numberOfPrisms);
  // fprintf(stderr, "C numberOfPoints: %d\n", hMesh.numberOfPoints); 
  // fprintf(stderr, "C numberOfSamples: %d\n\n", hMesh.numberOfSamples);

  // Solution vector
  std::vector<double> *ase = new std::vector<double>(hMesh.numberOfSamples * sigmaE->size(), 0);

  // Run Experiment
  for(int i=1; i < argc; ++i){
    if(strncmp(argv[i], "--mode=", 6) == 0){
      if(strstr(argv[i], "ray_propagation_gpu") != 0 ){
	// threads and blocks will be set in the following function (by reference)
	CUDA_CHECK_RETURN(cudaSetDevice(devices[0]));
	runtime = calcDndtAse(threads, 
			      blocks, 
			      raysPerSample,
			      dMesh[0],
			      hMesh,
			      betaCells,
			      nTot,
			      sigmaA,
			      sigmaE,
			      crystalFluorescence,
			      ase
			      );
	strcpy(runmode, "Ray Propagation New GPU");
	break;
      }
      else if(strstr(argv[i], "for_loops") != 0){
	// threads and blocks will be set in the following function (by reference)
	runtime = forLoopsClad(
			ase,
			raysPerSample,
			&hMesh,
			betaCells,
			nTot,
			sigmaA->at(0),
			sigmaE->at(0),
			hMesh.numberOfPoints,
			hMesh.numberOfTriangles,
			hMesh.numberOfLevels,
			hMesh.thickness,
			crystalFluorescence);
	strcpy(runmode, "For Loops");
	break;
      }
      else{
	fprintf(stderr, "C Please specify the runmode with --mode=\n");
	return 1;
      }
    
    }
    else{
      fprintf(stderr, "C Please specify the runmode with --mode=\n");
      return 1;
    }


  }

  // Print Solution
  for(unsigned wave_i = 0; wave_i < sigmaE->size(); ++wave_i){
    fprintf(stderr, "\n\nC Solutions %d\n", wave_i);
    for(unsigned sample_i = 0; sample_i < ase->size(); ++sample_i){
      fprintf(stderr, "C ASE PHI of sample %d: %.80f\n", sample_i, ase->at(sample_i + hMesh.numberOfSamples * wave_i));
      if(silent){
	if(sample_i >= 10) break;
      }
    }
  }

  // Print statistics
  fprintf(stderr, "\n");
  fprintf(stderr, "C Statistics\n");
  fprintf(stderr, "C Prism             : %d\n", (int) hMesh.numberOfPrisms);
  fprintf(stderr, "C Samples           : %d\n", (int) ase->size());
  fprintf(stderr, "C Rays/Sample       : %d\n", raysPerSample);
  fprintf(stderr, "C Rays Total        : %zu\n", raysPerSample * ase->size());
  fprintf(stderr, "C GPU Blocks        : %d\n", blocks);
  fprintf(stderr, "C GPU Threads/Block : %d\n", threads);
  fprintf(stderr, "C GPU Threads Total : %d\n", threads * blocks);
  fprintf(stderr, "C Runmode           : %s \n", runmode);
  fprintf(stderr, "C Runtime           : %f s\n", runtime);
  fprintf(stderr, "\n");

  // Write experiment data
  writeToVtk(&hMesh, ase, "octrace.vtk");
  compareVtk(ase, compareLocation, hMesh.numberOfSamples);
  writeToVtk(&hMesh, ase, "octrace_compare.vtk");
  writeDndtAse(ase);

  // Free memory
  delete betaCells;
  delete sigmaE;
  delete sigmaA;

  return 0;
}


