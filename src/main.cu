// Libraries
#include <stdio.h> /* fprintf, memcpy, strstr, strcmp */
#include <assert.h> /* assert */
#include <string> /* string */
#include <vector> /* vector */
#include <stdlib.h> /* atoi */

// User header files
#include <calc_dndt_ase.h>
#include <parser.h>
#include <write_to_vtk.h>
#include <write_matlab_output.h>
#include <for_loops_clad.h>
#include <cudachecks.h>
#include <mesh.h>

#define MIN_COMPUTE_CAPABILITY_MAJOR 2
#define MIN_COMPUTE_CAPABILITY_MINOR 0
#define MAX_GPUS 1

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
unsigned getCorrectDevice(int verbose,unsigned **devices){
  int count = 0, candidate = 0;
  unsigned correctDevices = 0;
  cudaDeviceProp prop;
  int minMajor = MIN_COMPUTE_CAPABILITY_MAJOR;
  int minMinor = MIN_COMPUTE_CAPABILITY_MINOR;

  // Get number of devices
  CUDA_CHECK_RETURN( cudaGetDeviceCount(&count));

  // Check devices for compute capability
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
  return correctDevices;
}

int main(int argc, char **argv){
  unsigned raysPerSample = 0;
  std::string runmode("");
  std::string compareLocation("");
  float runtime = 0.0;
  unsigned blocks;
  unsigned threads;
  bool silent = false;
  bool writeVtk = false;
  unsigned *devices; // will be assigned in getCOrrectDevice();
  unsigned numberOfDevices=0;
  unsigned maxGpus = MAX_GPUS;
  int device = -1;
  int mode = -1;
  

  std::string experimentPath;

  // Wavelength data
  std::vector<double> *sigmaA = new std::vector<double>;
  std::vector<double> *sigmaE = new std::vector<double>;

  // Set/Test device to run experiment with
  numberOfDevices = getCorrectDevice(1,&devices);

  // Parse Commandline
  parseCommandLine(argc, argv, &raysPerSample, &experimentPath, &device, &silent,
		   &writeVtk, &compareLocation, &mode);

  // sanity checks
  if(checkParameterValidity(argc, raysPerSample, experimentPath, &device, mode)) return 1;

  // Parse wavelengths from files
  if(fileToVector(experimentPath + "sigma_a.txt", sigmaA)) return 1;
  if(fileToVector(experimentPath + "sigma_e.txt", sigmaE)) return 1;
  assert(sigmaA->size() == sigmaE->size());
  assert(maxGpus <= numberOfDevices);

  // Parse experientdata and fill mesh
  Mesh hMesh;
  Mesh *dMesh = new Mesh[maxGpus];
  if(Mesh::parseMultiGPU(&hMesh, &dMesh, experimentPath, devices, maxGpus)) return 1;

  // Solution vector
  std::vector<double> *dndtAse = new std::vector<double>(hMesh.numberOfSamples * sigmaE->size(), 0);
  std::vector<float> *phiAse = new std::vector<float>(hMesh.numberOfSamples * sigmaE->size(), 0);
  std::vector<double> *expectation = new std::vector<double>(hMesh.numberOfSamples * sigmaE->size(), 0);
  CUDA_CHECK_RETURN( cudaSetDevice(devices[device])); 
  // Run Experiment
  switch(mode){
    case 0:
      // threads and blocks will be set in the following function (by reference)
      runtime = calcDndtAse(threads, 
			    blocks, 
			    raysPerSample,
			    dMesh[device],
			    hMesh,
			    *sigmaA,
			    *sigmaE,
			    dndtAse,
			    phiAse,
			    expectation
          );
      runmode="Ray Propagation New GPU";
      break;
    case 1:
      // threads and blocks will be set in the following function (by reference)
      runtime = forLoopsClad(
          dndtAse,
          raysPerSample,
          &hMesh,
          hMesh.betaCells,
          hMesh.nTot,
          sigmaA->at(0),
          sigmaE->at(0),
          hMesh.numberOfPoints,
          hMesh.numberOfTriangles,
          hMesh.numberOfLevels,
          hMesh.thickness,
          hMesh.crystalFluorescence);
      runmode = "For Loops";
      break;
  }

  // Print Solutions
  for(unsigned wave_i = 0; wave_i < sigmaE->size(); ++wave_i){
    fprintf(stderr, "\n\nC Solutions %d\n", wave_i);
    for(unsigned sample_i = 0; sample_i < dndtAse->size(); ++sample_i){
      int sampleOffset = sample_i + hMesh.numberOfSamples * wave_i;
      fprintf(stderr, "C Dndt ASE[%d]: %.80f %.10f\n", sample_i, dndtAse->at(sampleOffset), expectation->at(sampleOffset));
      if(silent){
        if(sample_i >= 10) break;
      }
    }
  }

  // Print statistics
  fprintf(stderr, "\n");
  fprintf(stderr, "C Statistics\n");
  fprintf(stderr, "C Prism             : %d\n", (int) hMesh.numberOfPrisms);
  fprintf(stderr, "C Samples           : %d\n", (int) dndtAse->size());
  fprintf(stderr, "C Rays/Sample       : %d\n", raysPerSample);
  fprintf(stderr, "C Rays Total        : %zu\n", raysPerSample * dndtAse->size());
  fprintf(stderr, "C GPU Blocks        : %d\n", blocks);
  fprintf(stderr, "C GPU Threads/Block : %d\n", threads);
  fprintf(stderr, "C GPU Threads Total : %d\n", threads * blocks);
  fprintf(stderr, "C Runmode           : %s \n", runmode.c_str());
  fprintf(stderr, "C Runtime           : %f s\n", runtime);
  fprintf(stderr, "\n");

  // Write experiment data
  std::vector<unsigned> *mockupN_rays = new std::vector<unsigned>(sigmaE->size(),1);
  writeMatlabOutput(
		  phiAse,
		  mockupN_rays,
		  expectation,
		  sigmaE->size(),
		  hMesh.numberOfSamples);

  if(writeVtk) writeToVtk(&hMesh, dndtAse, "octrace_dndt");
  if(compareLocation!="") {
	  std::vector<double> compareAse = compareVtk(*dndtAse, compareLocation, hMesh.numberOfSamples);
	  if(writeVtk) writeToVtk(&hMesh, dndtAse, "octrace_compare");
  }
  if(writeVtk) writeToVtk(&hMesh, expectation, "octrace_expectation");

  // Free memory
  delete devices;
  delete sigmaE;
  delete sigmaA;
  delete dndtAse;
  delete phiAse;
  delete expectation;
  delete mockupN_rays;
  cudaFree(dMesh);

  return 0;
}


