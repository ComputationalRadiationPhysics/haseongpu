// Libraries
#include <stdio.h> /* fprintf, memcpy, strstr, strcmp */
#include <assert.h> /* assert */
#include <string> /* string */
#include <vector> /* vector */
#include <stdlib.h> /* atoi */
#include <pthread.h> /* pthread_t, pthread_join */
#include <algorithm> /* max_element */

// User header files
#include <calc_phi_ase.h>
#include <calc_phi_ase_threaded.h>
#include <calc_phi_ase_mpi.h>
#include <parser.h>
#include <write_to_vtk.h>
#include <write_matlab_output.h>
#include <for_loops_clad.h>
#include <cudachecks.h>
#include <mesh.h>
#include <test_environment.h>

#include <write_value_to_file.h>

#include <logging.h>

#define MIN_COMPUTE_CAPABILITY_MAJOR 2
#define MIN_COMPUTE_CAPABILITY_MINOR 0
unsigned verbosity = V_ERROR | V_INFO | V_WARNING; // extern through logging.h


/** 
 * @brief Queries for devices on the running mashine and collects
 *        them on the devices array. Set the first device in this 
 *        array as computaion-device. On Errors the programm will
 *        be stoped by exit(). Otherwise you can set the device by command
 *        line parameter --device=
 * 
 * @param verbose > 0 prints debug output
 * 
 * @return vector of possible devices
 */
std::vector<unsigned> getCorrectDevice(unsigned maxGpus){
  cudaDeviceProp prop;
  int minMajor = MIN_COMPUTE_CAPABILITY_MAJOR;
  int minMinor = MIN_COMPUTE_CAPABILITY_MINOR;
  int count;
  std::vector<unsigned> devices;

  // Get number of devices
  CUDA_CHECK_RETURN( cudaGetDeviceCount(&count));

  // Check devices for compute capability and if device is busy
  int devicesAllocated = 0;
  for(int i=0; i < count; ++i){
    CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, i) );
    if( (prop.major > minMajor) || (prop.major == minMajor && prop.minor >= minMinor) ){
      cudaSetDevice(i);
      int* occupy;
      if(cudaMalloc((void**) &occupy, sizeof(int)) == cudaSuccess){
        devices.push_back(i);
	devicesAllocated++;
	if(devicesAllocated == maxGpus)
	  break;

      }

    }

  }

  if(devices.size() == 0){
    dout(V_ERROR) << "None of the free CUDA-capable devices is sufficient!" << std::endl;
    exit(1);
  }

  cudaSetDevice(devices.at(0));

  dout(V_INFO) << "Found " << int(devices.size()) << " available CUDA devices with Compute Capability >= " << minMajor << "." << minMinor << "):" << std::endl;
  for(unsigned i=0; i<devices.size(); ++i){
    CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, devices[i]) );
    dout(V_INFO) << "[" << devices[i] << "] " << prop.name << " (Compute Capability " << prop.major << "." << prop.minor << ")" << std::endl;
  }


  return devices;

}

double calcDndtAse(const Mesh& mesh, const double sigmaA, const double sigmaE, const float phiAse, const unsigned sample_i){
  double gain_local = mesh.nTot * mesh.betaCells[sample_i] * (sigmaE + sigmaA) - double(mesh.nTot * sigmaA);
  return gain_local * phiAse / mesh.crystalFluorescence;
}


int main(int argc, char **argv){
  unsigned raysPerSample = 0;
  unsigned maxRaysPerSample = 0;
  float maxExpectation = 0;
  float  avgExpectation = 0;
  unsigned highExpectation = 0;
  std::string runmode("");
  std::string compareLocation("");
  float runtime = 0.0;
  bool silent = false;
  bool writeVtk = false;
  bool useReflections = false;
  std::vector<unsigned> devices; // will be assigned in getCOrrectDevice();
  unsigned maxGpus = 0;
  RunMode mode = NONE;
  int minSampleRange = 0;
  int maxSampleRange = 0;

  std::string experimentPath;
  verbosity = 63;

  // Wavelength data
  std::vector<double> sigmaA;
  std::vector<double> sigmaE;
  std::vector<float> mseThreshold;

  // Parse Commandline
  parseCommandLine(argc, argv, &raysPerSample, &maxRaysPerSample, &experimentPath, &silent,
      &writeVtk, &compareLocation, &mode, &useReflections, &maxGpus, &minSampleRange, &maxSampleRange);

  // Set/Test device to run experiment with
  devices = getCorrectDevice(maxGpus);

  // sanity checks
  if(checkParameterValidity(argc, raysPerSample, &maxRaysPerSample, experimentPath, devices.size(), mode, &maxGpus, minSampleRange, maxSampleRange)) return 1;

  // Parse wavelengths from files
  if(fileToVector(experimentPath + "sigma_a.txt", &sigmaA)) return 1;
  if(fileToVector(experimentPath + "sigma_e.txt", &sigmaE)) return 1;
  if(fileToVector(experimentPath + "mse_threshold.txt", &mseThreshold)) return 1;
  assert(sigmaA.size() == sigmaE.size());
  assert(mseThreshold.size() == sigmaE.size());


  // Parse experientdata and fill mesh
  Mesh hMesh;
  std::vector<Mesh> dMesh(maxGpus);

  if(Mesh::parseMultiGPU(hMesh, dMesh, experimentPath, devices, maxGpus)) return 1;

  // Solution vector
  std::vector<double> dndtAse(hMesh.numberOfSamples * sigmaE.size(), 0);
  std::vector<float>  phiAse(hMesh.numberOfSamples * sigmaE.size(), 0);
  std::vector<double> mse(hMesh.numberOfSamples * sigmaE.size(), 1000);
  std::vector<unsigned> totalRays(hMesh.numberOfSamples * sigmaE.size(), 0);

  // Run Experiment
  std::vector<float> runtimes(maxGpus, 0);
  std::vector<pthread_t> threadIds(maxGpus, 0);

  unsigned samplesPerNode = maxSampleRange-minSampleRange+1;
  float samplePerGpu = samplesPerNode / (float) maxGpus;
  switch(mode){
  case RAY_PROPAGATION_GPU:
    for(unsigned gpu_i = 0; gpu_i < maxGpus; ++gpu_i){
      unsigned minSample_i = gpu_i * samplePerGpu;
      unsigned maxSample_i = min((float)samplesPerNode, (gpu_i + 1) * samplePerGpu);

      minSample_i += minSampleRange;
      maxSample_i += minSampleRange; 

      threadIds[gpu_i] = calcPhiAseThreaded( raysPerSample,
					     maxRaysPerSample,
					     dMesh.at(devices.at(gpu_i)),
					     hMesh,
					     sigmaA,
					     sigmaE,
					     mseThreshold,
					     useReflections,
					     phiAse, 
					     mse, 
					     totalRays,
					     devices.at(gpu_i),
					     minSample_i,
					     maxSample_i,
					     runtimes.at(gpu_i)
					     );

    }
    joinAll(threadIds);
    for(std::vector<float>::iterator it = runtimes.begin(); it != runtimes.end(); ++it){
      runtime = max(*it, runtime);
    }
    cudaDeviceReset();      
    runmode="Ray Propagation GPU";
    break;

  case RAY_PROPAGATION_MPI:
    runtime = calcPhiAseMPI( raysPerSample,
			     maxRaysPerSample,
			     dMesh.at(0),
			     hMesh,
			     sigmaA,
			     sigmaE,
			     mseThreshold,
			     useReflections,
			     phiAse,
			     mse,
			     totalRays,
			     devices.at(0),
			     maxSampleRange
			    );
    runmode = "RAY PROPAGATION MPI";
    break;

  case FOR_LOOPS:
    runtime = forLoopsClad( &dndtAse,
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

  case TEST:
    testEnvironment(raysPerSample,
		    maxRaysPerSample,
		    dMesh.at(0),
		    hMesh,
		    sigmaA,
		    sigmaE,
		    mseThreshold.at(0),
		    useReflections,
		    dndtAse,
		    phiAse,
		    mse
		    );
    cudaDeviceReset();
    runmode="Test Environment";
    break;
  default:
    exit(0);
  }

  // Filter maxExpectation
  for(std::vector<double>::iterator it = mse.begin(); it != mse.end(); ++it){
    maxExpectation = max(maxExpectation, *it);
    avgExpectation += *it;
    if(*it > mseThreshold.at(0))
      highExpectation++;
  }
  avgExpectation /= mse.size();


   // Print Solutions
   for(unsigned wave_i = 0; wave_i < sigmaE.size(); ++wave_i){
     dout(V_DEBUG) << "\n\nSolutions " <<  wave_i << std::endl;
     for(unsigned sample_i = 0; sample_i < hMesh.numberOfSamples; ++sample_i){
       int sampleOffset = sample_i + hMesh.numberOfSamples * wave_i;
       dndtAse.at(sampleOffset) = calcDndtAse(hMesh, sigmaA.at(wave_i), sigmaE.at(wave_i), phiAse.at(sampleOffset), sample_i);
       if(sample_i <=10)
         dout(V_DEBUG) << "Dndt ASE[" << sample_i << "]: " << dndtAse.at(sampleOffset) << " " << mse.at(sampleOffset) << std::endl;
     }
     for(unsigned sample_i = 0; sample_i < hMesh.numberOfSamples; ++sample_i){
       int sampleOffset = sample_i + hMesh.numberOfSamples * wave_i;
       dout(V_DEBUG) << "PHI ASE[" << sample_i << "]: " << phiAse.at(sampleOffset) << " " << mse.at(sampleOffset) <<std::endl;
       if(sample_i >= 10) break;
     }
   }

  // Compare with vtk
  // if(compareLocation!="") {
  //   std::vector<double> compareAse = compareVtk(dndtAse, compareLocation, hMesh.numberOfSamples);

  // }

  //Print statistics
  dout(V_STAT) << "\nStatistics\n" << std::endl;
  dout(V_STAT) << "Prism             : " << (int) hMesh.numberOfPrisms << std::endl;
  dout(V_STAT) << "Samples           : " << (int) dndtAse.size() << std::endl;
  dout(V_STAT) << "MSE threshold     : " << *(std::max_element(mseThreshold.begin(),mseThreshold.end())) << std::endl;
  dout(V_STAT) << "max. MSE          : " << maxExpectation << std::endl;
  dout(V_STAT) << "avg. MSE          : " << avgExpectation << std::endl;
  dout(V_STAT) << "too high MSE      : " << highExpectation << std::endl;
  dout(V_STAT) << "Runmode           : " << runmode.c_str() << std::endl;
  dout(V_STAT) << "Nr of GPUs        : " << maxGpus << std::endl;
  dout(V_STAT) << "Runtime           : " << runtime << "s\n\n" << std::endl;

  // Write experiment data
  // writeMatlabOutput(
  //     "output/",
  //     phiAse,
  //     totalRays,
  //     mse,
  //     sigmaE.size(),
  //     hMesh.numberOfSamples,
  //     hMesh.numberOfLevels
  //     );

  // Write output in single files
  for(unsigned wave_i=0 ; wave_i < sigmaE.size() ; ++wave_i){
    for(unsigned sample_i = minSampleRange; sample_i < minSampleRange+samplesPerNode ; sample_i++){
      unsigned sampleOffset = sample_i + hMesh.numberOfSamples * wave_i;
      writeValueToFile(phiAse.at(sampleOffset),"output/results/","wavelength",wave_i,"sample",sample_i);
    }
  }

  // FOR OUTPUT
  std::vector<double> tmpPhiAse(phiAse.begin(), phiAse.end());
  std::vector<double> tmpTotalRays(totalRays.begin(), totalRays.end());

  if(writeVtk) writeToVtk(hMesh, dndtAse, "vtk/dndt", raysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, runtime);
  if(writeVtk) writeToVtk(hMesh, tmpPhiAse, "vtk/phiase", raysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, runtime);
  if(writeVtk) writeToVtk(hMesh, mse, "vtk/mse", raysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, runtime);
  if(writeVtk) writeToVtk(hMesh, tmpTotalRays, "vtk/total_rays", raysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, runtime);

  return 0;
}


