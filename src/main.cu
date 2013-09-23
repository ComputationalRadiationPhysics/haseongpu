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
std::vector<unsigned> getCorrectDevice(int verbose){
  cudaDeviceProp prop;
  int minMajor = MIN_COMPUTE_CAPABILITY_MAJOR;
  int minMinor = MIN_COMPUTE_CAPABILITY_MINOR;
  int count;
  std::vector<unsigned> devices;

  // Get number of devices
  CUDA_CHECK_RETURN( cudaGetDeviceCount(&count));

  // Check devices for compute capability and if device is busy
  for(int i=0; i<count; ++i){
    CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, i) );
    if( (prop.major > minMajor) || (prop.major == minMajor && prop.minor >= minMinor) ){
      cudaSetDevice(i);
      int* test;
      if(cudaMalloc((void**) &test, sizeof(int)) == cudaSuccess){
        devices.push_back(i);
        cudaFree(test);
        cudaDeviceReset();
      }
    }
  }

  if(devices.size() == 0){
    fprintf(stderr,"\nNone of the free CUDA-capable devices is sufficient!\n");
    exit(1);
  }

  cudaSetDevice(devices.at(0));

  if(verbose > 0){
    fprintf(stderr,"\nFound %d available CUDA devices with Compute Capability >= %d.%d):\n", int(devices.size()), minMajor,minMinor);
    for(unsigned i=0; i<devices.size(); ++i){
      CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, devices[i]) );
      fprintf(stderr,"[%d] %s (Compute Capability %d.%d)\n", devices[i], prop.name, prop.major, prop.minor);
    }
  }

  return devices;

}

int main(int argc, char **argv){
  unsigned raysPerSample = 0;
  unsigned maxRaysPerSample = 0;
  std::string runmode("");
  std::string compareLocation("");
  float runtime = 0.0;
  unsigned blocks;
  unsigned threads;
  bool silent = false;
  bool writeVtk = false;
  bool useReflections = false;
  float expectationThreshold = 0;
  std::vector<unsigned> devices; // will be assigned in getCOrrectDevice();
  unsigned maxGpus = MAX_GPUS;
  int device = -1;
  int mode = -1;

  std::string experimentPath;

  // Wavelength data
  std::vector<double> sigmaA;
  std::vector<double> sigmaE;

  // Set/Test device to run experiment with
  devices = getCorrectDevice(1);

  // Parse Commandline
  parseCommandLine(argc, argv, &raysPerSample, &maxRaysPerSample, &experimentPath, &device, &silent,
      &writeVtk, &compareLocation, &mode, &useReflections, &expectationThreshold);

  // sanity checks
  if(checkParameterValidity(argc, raysPerSample, &maxRaysPerSample, experimentPath, &device, devices.size(), mode, &expectationThreshold)) return 1;

  // Parse wavelengths from files
  if(fileToVector(experimentPath + "sigma_a.txt", &sigmaA)) return 1;
  if(fileToVector(experimentPath + "sigma_e.txt", &sigmaE)) return 1;
  assert(sigmaA.size() == sigmaE.size());
  assert(maxGpus <= devices.size());

  // Parse experientdata and fill mesh
  Mesh hMesh;
  std::vector<Mesh> dMesh(maxGpus);
  
  if(Mesh::parseMultiGPU(hMesh, dMesh, experimentPath, devices, maxGpus)) return 1;

  // Solution vector
  std::vector<double> dndtAse(hMesh.numberOfSamples * sigmaE.size(), 0);
  std::vector<float>  phiAse(hMesh.numberOfSamples * sigmaE.size(), 0);
  std::vector<double> expectation(hMesh.numberOfSamples * sigmaE.size(), 0);
  CUDA_CHECK_RETURN( cudaSetDevice(devices.at(device))); 

  fprintf(stderr, "reflectionAngle: %f\n",hMesh.getReflectionAngle(-1));
  fprintf(stderr, "reflectionAngle: %f\n",hMesh.getReflectionAngle(1));
  fprintf(stderr, "maxreflections: %d\n",hMesh.getMaxReflections());
  fprintf(stderr, "sigma_A: %f %f\n",sigmaA[0],sigmaA[1]);
  fprintf(stderr, "sigma_E: %f %f\n",sigmaE[0],sigmaE[1]);
  
  // Run Experiment
  switch(mode){
    case 0:
      // threads and blocks will be set in the following function (by reference)
      runtime = calcDndtAse(threads, 
			    blocks, 
			    raysPerSample,
			    maxRaysPerSample,
			    dMesh.at(device),
			    hMesh,
			    sigmaA,
			    sigmaE,
			    expectationThreshold,
			    useReflections,
			    dndtAse,
			    phiAse,
			    expectation
			    );
      cudaDeviceReset();
      runmode="Ray Propagation New GPU";
      break;
    case 1:
      // threads and blocks will be set in the following function (by reference)
      runtime = forLoopsClad(
          &dndtAse,
          raysPerSample,
          &hMesh,
          hMesh.betaCells,
          hMesh.nTot,
          sigmaA.at(0),
          sigmaE.at(0),
          hMesh.numberOfPoints,
          hMesh.numberOfTriangles,
          hMesh.numberOfLevels,
          hMesh.thickness,
          hMesh.crystalFluorescence);
      runmode = "For Loops";
      break;
  }


  // Print Solutions
  for(unsigned wave_i = 0; wave_i < sigmaE.size(); ++wave_i){
    fprintf(stderr, "\n\nC Solutions %d\n", wave_i);
    for(unsigned sample_i = 0; sample_i < hMesh.numberOfSamples; ++sample_i){
      int sampleOffset = sample_i + hMesh.numberOfSamples * wave_i;
      fprintf(stderr, "C Dndt ASE[%d]: %.80f %.10f\n", sample_i, dndtAse.at(sampleOffset), expectation.at(sampleOffset));
      if(silent){
        if(sample_i >= 10) break;
      }
    }
    for(unsigned sample_i = 0; sample_i < hMesh.numberOfSamples; ++sample_i){
      int sampleOffset = sample_i + hMesh.numberOfSamples * wave_i;
      fprintf(stderr, "C PHI ASE[%d]: %.80f %.10f\n", sample_i, phiAse.at(sampleOffset), expectation.at(sampleOffset));
      if(silent){
        if(sample_i >= 10) break;
      }
	}
  }

  // Print statistics
  fprintf(stderr, "\n");
  fprintf(stderr, "C Statistics\n");
  fprintf(stderr, "C Prism             : %d\n", (int) hMesh.numberOfPrisms);
  fprintf(stderr, "C Samples           : %d\n", (int) dndtAse.size());
  fprintf(stderr, "C Rays/Sample       : %d\n", raysPerSample);
  fprintf(stderr, "C Rays Total        : %zu\n", raysPerSample * dndtAse.size());
  fprintf(stderr, "C GPU Blocks        : %d\n", blocks);
  fprintf(stderr, "C GPU Threads/Block : %d\n", threads);
  fprintf(stderr, "C GPU Threads Total : %d\n", threads * blocks);
  fprintf(stderr, "C Runmode           : %s \n", runmode.c_str());
  fprintf(stderr, "C Runtime           : %f s\n", runtime);
  fprintf(stderr, "\n");

  // Write experiment data
  std::vector<unsigned> mockupN_rays(sigmaE.size() * hMesh.numberOfSamples, 1);
  writeMatlabOutput(
		  experimentPath,
		  phiAse,
		  mockupN_rays,
		  expectation,
		  sigmaE.size(),
		  hMesh.numberOfSamples,
		  hMesh.numberOfLevels
      );

  if(writeVtk) writeToVtk(hMesh, dndtAse, experimentPath + "octrace_dndt", raysPerSample, maxRaysPerSample, expectationThreshold, runtime);
  if(compareLocation!="") {
	  std::vector<double> compareAse = compareVtk(dndtAse, compareLocation, hMesh.numberOfSamples);
	  if(writeVtk) writeToVtk(hMesh, dndtAse, experimentPath + "octrace_compare", raysPerSample, maxRaysPerSample, expectationThreshold, runtime);
  }
  if(writeVtk) writeToVtk(hMesh, expectation, experimentPath + "octrace_expectation", raysPerSample, maxRaysPerSample, expectationThreshold, runtime);

  return 0;
}


