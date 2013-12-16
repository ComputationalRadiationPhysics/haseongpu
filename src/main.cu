// Libraries
#include <stdio.h> /* fprintf, memcpy, strstr, strcmp */
#include <assert.h> /* assert */
#include <string> /* string */
#include <vector> /* vector */
#include <stdlib.h> /* atoi */
#include <pthread.h> /* pthread_t, pthread_join */
#include <algorithm> /* max_element */
#include <numeric> /* accumulate*/

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

#include <logging.h>
#include <ray_histogram.h>

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
  unsigned devicesAllocated = 0;
  for(int i=0; i < count; ++i){
    CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, i) );
    if( (prop.major > minMajor) || (prop.major == minMajor && prop.minor >= minMinor) ){
      cudaSetDevice(i);
      int* occupy; //TODO: occupy gets allocated, but never cudaFree'd -> small memory leak!
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
  unsigned maxRepetitions = 4;
  float maxMSE = 0;
  float  avgMSE = 0;
  unsigned highMSE = 0;
  std::string runmode("");
  std::string compareLocation("");
  float runtime = 0.0;
  bool writeVtk = false;
  bool useReflections = false;
  std::vector<unsigned> devices; // will be assigned in getCOrrectDevice();
  unsigned maxGpus = 0;
  RunMode mode = NONE;
  int minSampleRange = 0;
  int maxSampleRange = 0;
  time_t starttime   = time(0);
  unsigned usedGpus  = 0;

  std::string inputPath;
  std::string outputPath;
  verbosity = 31; //ALL //TODO: remove in final code

  // Wavelength data
  std::vector<double> sigmaA;
  std::vector<double> sigmaE;
  std::vector<float> mseThreshold;

  // Parse Commandline
  parseCommandLine(argc, argv, &raysPerSample, &maxRaysPerSample, &inputPath,
		   &writeVtk, &compareLocation, &mode, &useReflections, &maxGpus, &minSampleRange, &maxSampleRange, &maxRepetitions, &outputPath);

  // Set/Test device to run experiment with
  //TODO: this call takes a LOT of time (2-5s). Can this be avoided?
  //TODO: maybe move this to a place where GPUs are actually needed (for_loops_clad doesn't even need GPUs!)
  devices = getCorrectDevice(maxGpus);


  // sanity checks
  if(checkParameterValidity(argc, raysPerSample, &maxRaysPerSample, inputPath, devices.size(), mode, &maxGpus, minSampleRange, maxSampleRange, maxRepetitions, outputPath)) return 1;

  // Parse wavelengths from files
  if(fileToVector(inputPath + "sigma_a.txt", &sigmaA)) return 1;
  if(fileToVector(inputPath + "sigma_e.txt", &sigmaE)) return 1;
  if(fileToVector(inputPath + "mse_threshold.txt", &mseThreshold)) return 1;
  assert(sigmaA.size() == sigmaE.size());
  assert(mseThreshold.size() == sigmaE.size());


  // Parse experientdata and fill mesh
  Mesh hMesh;
  std::vector<Mesh> dMesh(maxGpus);

  // TODO: split into hMesh and dMesh parsing 
  // -> parse dMesh only where needed
  if(Mesh::parseMultiGPU(hMesh, dMesh, inputPath, devices, maxGpus)) return 1;

  // Solution vector
  std::vector<double> dndtAse(hMesh.numberOfSamples * sigmaE.size(), 0);
  std::vector<float>  phiAse(hMesh.numberOfSamples * sigmaE.size(), 0);
  std::vector<double> mse(hMesh.numberOfSamples * sigmaE.size(), 1000);
  std::vector<unsigned> totalRays(hMesh.numberOfSamples * sigmaE.size(), 0);

  // for(unsigned i = 0; i < hMesh.numberOfPrisms; ++i){
  //   dout(V_DEBUG) << i << " " << hMesh.betaValues[i] << std::endl;
  // }

  // Run Experiment
  std::vector<pthread_t> threadIds(maxGpus, 0);
  std::vector<float> runtimes(maxGpus, 0);
  switch(mode){
    case RAY_PROPAGATION_GPU:
      for(unsigned gpu_i = 0; gpu_i < maxGpus; ++gpu_i){
        const unsigned samplesPerNode = maxSampleRange-minSampleRange+1;
        const float samplePerGpu = samplesPerNode / (float) maxGpus;
        unsigned minSample_i = gpu_i * samplePerGpu;
        unsigned maxSample_i = min((float)samplesPerNode, (gpu_i + 1) * samplePerGpu);

        minSample_i += minSampleRange;
        maxSample_i += minSampleRange; 

        threadIds[gpu_i] = calcPhiAseThreaded( raysPerSample,
            maxRaysPerSample,
            maxRepetitions,
            dMesh.at(gpu_i),
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
      usedGpus = maxGpus;
      for(std::vector<float>::iterator it = runtimes.begin(); it != runtimes.end(); ++it){
        runtime = max(*it, runtime);
      }
      cudaDeviceReset();      
      runmode="Ray Propagation GPU";
      break;

    case RAY_PROPAGATION_MPI:
      usedGpus = calcPhiAseMPI( raysPerSample,
          maxRaysPerSample,
          maxRepetitions,
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

    case FOR_LOOPS: //Possibly deprecated!
      // TODO: make available for MPI?
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


  if(verbosity & V_DEBUG){
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
  }

  // Compare with vtk
  // if(compareLocation!="") {
  //   std::vector<double> compareAse = compareVtk(dndtAse, compareLocation, hMesh.numberOfSamples);

  // }


  // Write experiment data
  // output folder has to be the same as TMP_FOLDER in the calling MatLab script
  writeMatlabOutput(
      outputPath,
      phiAse,
      totalRays,
      mse,
      sigmaE.size(),
      hMesh.numberOfSamples,
      hMesh.numberOfLevels
      );


  // FOR OUTPUT
  if(writeVtk){
    std::vector<double> tmpPhiAse(phiAse.begin(), phiAse.end());
    std::vector<double> tmpTotalRays(totalRays.begin(), totalRays.end());

    writeToVtk(hMesh, dndtAse, outputPath + "vtk/dndt", raysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, runtime);
    writeToVtk(hMesh, tmpPhiAse, outputPath + "vtk/phiase", raysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, runtime);
    writeToVtk(hMesh, mse, outputPath + "vtk/mse", raysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, runtime);
    writeToVtk(hMesh, tmpTotalRays, outputPath + "vtk/total_rays", raysPerSample, maxRaysPerSample, mseThreshold.at(0), useReflections, runtime);
  }

  if(verbosity & V_STAT){
    // Filter maxMSE
    for(std::vector<double>::iterator it = mse.begin(); it != mse.end(); ++it){
      maxMSE = max(maxMSE, *it);
      avgMSE += *it;
      if(*it > mseThreshold.at(0))
        highMSE++;
    }
    avgMSE /= mse.size();

    //Print statistics
    std::cout.imbue(std::locale(""));
    dout(V_STAT | V_NOLABEL) << std::endl;
    dout(V_STAT) << "=== Statistics ===" << std::endl;
    dout(V_STAT) << "Runmode           : " << runmode.c_str() << std::endl;
    dout(V_STAT) << "Prisms            : " << (int) hMesh.numberOfPrisms << std::endl;
    dout(V_STAT) << "Samples           : " << (int) dndtAse.size() << std::endl;
    dout(V_STAT) << "Wavelength        : " << (int) sigmaE.size() << std::endl;
    dout(V_STAT) << "RaysPerSample     : " << raysPerSample;
    if(maxRaysPerSample > raysPerSample) { dout(V_STAT | V_NOLABEL) << " - " << maxRaysPerSample << " (adaptive)"; }
    dout(V_STAT | V_NOLABEL) << std::endl;
    dout(V_STAT) << "sum(totalRays)    : " << std::accumulate(totalRays.begin(), totalRays.end(), 0.) << std::endl;
    dout(V_STAT) << "MSE threshold     : " << *(std::max_element(mseThreshold.begin(),mseThreshold.end())) << std::endl;
    dout(V_STAT) << "max. MSE          : " << maxMSE << std::endl;
    dout(V_STAT) << "avg. MSE          : " << avgMSE << std::endl;
    dout(V_STAT) << "too high MSE      : " << highMSE << std::endl;
    dout(V_STAT) << "Nr of GPUs        : " << usedGpus << std::endl;
    dout(V_STAT) << "Runtime           : " << difftime(time(0),starttime) << "s" << std::endl;
    dout(V_STAT) << std::endl;
    if(maxRaysPerSample > raysPerSample){
      dout(V_STAT) << "=== Sampling resolution as Histogram ===" << std::endl;
      ray_histogram(totalRays,raysPerSample,maxRaysPerSample,highMSE);
    }
    dout(V_STAT) << std::endl;

  }
  return 0;

}
