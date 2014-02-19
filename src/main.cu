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
#include <parser.h> /* RunMode */
#include <write_to_vtk.h>
#include <write_matlab_output.h>
#include <for_loops_clad.h>
#include <cudachecks.h>
#include <mesh.h>

#include <logging.h>
#include <ray_histogram.h>

#define MIN_COMPUTE_CAPABILITY_MAJOR 2
#define MIN_COMPUTE_CAPABILITY_MINOR 0
#define MAX_INTERPOLATION 1000
#define LAMBDA_START 905
#define LAMBDA_STOP 1095

// default without V_DEBUG
unsigned verbosity = V_ERROR | V_INFO | V_WARNING | V_PROGRESS | V_STAT; // extern through logging.h

/** 
 * @brief Queries for devices on the running mashine and collects
 *        them on the devices array. Set the first device in this 
 *        array as computation-device. On Errors the programm will
 *        be stoped by exit(). 
 * 
 * @param maxGpus max. devices which should be allocated
 * @return vector of possible devices
 */
std::vector<unsigned> getFreeDevices(unsigned maxGpus){
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

  // Exit if no device was found
  if(devices.size() == 0){
    dout(V_ERROR) << "None of the free CUDA-capable devices is sufficient!" << std::endl;
    exit(1);
  }

  // Print device information
  cudaSetDevice(devices.at(0));
  dout(V_INFO) << "Found " << int(devices.size()) << " available CUDA devices with Compute Capability >= " << minMajor << "." << minMinor << "):" << std::endl;
  for(unsigned i=0; i<devices.size(); ++i){
    CUDA_CHECK_RETURN( cudaGetDeviceProperties(&prop, devices[i]) );
    dout(V_INFO) << "[" << devices[i] << "] " << prop.name << " (Compute Capability " << prop.major << "." << prop.minor << ")" << std::endl;
  }

  return devices;

}

/** 
 * @brief Calculates dndt ASE from phi ASE values
 * 
 * @param mesh needed for some constants
 * @param sigmaA absorption
 * @param sigmaE emission
 * @param phiAse results from calcPhiAse
 * @param sample_i index of sample point
 * @return dndtAse
 *
 */
double calcDndtAse(const Mesh& mesh, const double sigmaA, const double sigmaE, const float phiAse, const unsigned sample_i){
  double gain_local = mesh.nTot * mesh.betaCells[sample_i] * (sigmaE + sigmaA) - double(mesh.nTot * sigmaA);
  return gain_local * phiAse / mesh.crystalFluorescence;
}

/**
 * @brief Returns the index of an value in vector v,
 *        that is smaller than t.
 *        
 * @param v vector
 * @param t bigger value
 * @return index of smaller value
 *         otherwise 0
 *
 */
unsigned getNextSmallerIndex(std::vector<double> v, double t){
  unsigned index = 0;
  for(unsigned i = 0; i < v.size(); ++i){
    if(v.at(i) < t) index = i;
    else break;
  }
  return index;
}

/**
 * @brief Returns the index of an value in vector v,
 *        that is bigger than t.
 *        
 * @param v vector
 * @param t smaller value
 * @return index of smaller value
 *         otherwise 0
 *
 */
unsigned getNextBiggerIndex(std::vector<double> v, double t){
  for(unsigned i = 0; i < v.size(); ++i){
    if(v.at(i) > t)
      return i;
  }
  return 0;
}

/**
 * @brief Interpolates the values of sigma_y to n values(interpolation range) linear. 
 *        With the assumption, they are distributed between lambda_start
 *        and lambda_stop equidistant. For Example could you interpolate
 *        100 sigma values to 1000 sigma values, to reach a better resolution.
 *
 * @param sigma_y             y values
 * @param interpolation_range number of interpolated values
 * @param lambda_start        start of x range
 * @param lambda_stop         stop of x range
 *
 * @return vector of linear interpolated values
 */
std::vector<double> interpolateWavelength(const std::vector<double> sigma_y, const unsigned interpolation_range, const double lambda_start, const double lambda_stop){
  assert(interpolation_range >= sigma_y.size());
  assert(lambda_stop >= lambda_start);

  // Monochromatic case
  if(sigma_y.size() == 1){
    return std::vector<double>(1, sigma_y.at(0));
  }

  std::vector<double> y(interpolation_range, 0);
  const double lambda_range = lambda_stop - lambda_start;
  assert(sigma_y.size() >= lambda_range);

  // Generate sigma_x
  std::vector<double> sigma_x;
  for(unsigned i = lambda_start; i <= lambda_stop; ++i){
    sigma_x.push_back(i);
  }
  
  for(unsigned i = 0; i < interpolation_range; ++i){
    double x = lambda_start + (i * (lambda_range / interpolation_range));

    // Get index of points before and after x
    double y1_i = getNextSmallerIndex(sigma_x, x);
    double y2_i = getNextBiggerIndex(sigma_x, x);
    int sigma_diff = y2_i - y1_i;

    if(sigma_diff == 1){
      // First point p1=(x1/y1) before x
      double x1 = lambda_start + y1_i;
      double y1 = sigma_y.at(y1_i);

      // Second point p2=(x2/y2) after x
      double x2 = lambda_start + y2_i;
      double y2 = sigma_y.at(y2_i);
      assert(sigma_y.size() >= y1_i);

      // linear function between p1 and p2 (y=mx+b)
      double m = (y2 - y1) / (x2 / x1);
      double b = y1 - (m * x1);

      // Interpolate y from linear function
      y.at(i) = m * x + b;

    }
    else if(sigma_diff == 2){
      // No interpolation needed
      y.at(i) = sigma_y.at(y1_i + 1);
    }
    else {
      dout(V_ERROR) << "Index of smaller and bigger sigma too seperated" << std::endl;
      exit(0);
    }
    
  }
  
  return y;
}


int main(int argc, char **argv){
  unsigned minRaysPerSample = 0;
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
  std::vector<unsigned> devices; 
  unsigned maxGpus = 0;
  RunMode mode = NONE;
  int minSampleRange = 0;
  int maxSampleRange = 0;
  time_t starttime   = time(0);
  unsigned usedGpus  = 0;

  std::string inputPath;
  std::string outputPath;
  double mseThreshold = 0;
  verbosity = 63; //ALL //TODO: remove in final code

  // Wavelength data
  std::vector<double> sigmaA;
  std::vector<double> sigmaE;

  // Parse Commandline
  parseCommandLine(argc, argv, &minRaysPerSample, &maxRaysPerSample, &inputPath,
		   &writeVtk, &compareLocation, &mode, &useReflections, &maxGpus, &minSampleRange, &maxSampleRange, &maxRepetitions, &outputPath, &mseThreshold);

  // Set/Test device to run experiment with
  //TODO: this call takes a LOT of time (2-5s). Can this be avoided?
  //TODO: maybe move this to a place where GPUs are actually needed (for_loops_clad doesn't even need GPUs!)
  devices = getFreeDevices(maxGpus);

  // sanity checks
  if(checkParameterValidity(argc, minRaysPerSample, &maxRaysPerSample, inputPath, devices.size(), mode, &maxGpus, minSampleRange, maxSampleRange, maxRepetitions, outputPath, &mseThreshold)) return 1;

  // Parse wavelengths from files
  if(fileToVector(inputPath + "sigma_a.txt", &sigmaA)) return 1;
  if(fileToVector(inputPath + "sigma_e.txt", &sigmaE)) return 1;
  assert(sigmaA.size() == sigmaE.size());

  // Interpolate sigmaA / sigmaE function
  std::vector<double> sigmaAInterpolated = interpolateWavelength(sigmaA, MAX_INTERPOLATION, LAMBDA_START, LAMBDA_STOP);
  std::vector<double> sigmaEInterpolated = interpolateWavelength(sigmaE, MAX_INTERPOLATION, LAMBDA_START, LAMBDA_STOP);

  // Calc max sigmaA / sigmaE
  double maxSigmaE = 0.0;
  double maxSigmaA = 0.0;
  for(unsigned i = 0; i < sigmaE.size(); ++i){
    if(sigmaE.at(i) > maxSigmaE){
      maxSigmaE = sigmaE.at(i);
      maxSigmaA = sigmaA.at(i);
    }
  }

  // Parse experientdata and fill mesh
  Mesh hMesh;
  std::vector<Mesh> dMesh(maxGpus);

  // TODO: split into hMesh and dMesh parsing 
  // -> parse dMesh only where needed
  if(Mesh::parseMultiGPU(hMesh, dMesh, inputPath, devices, maxGpus)) return 1;

  // Solution vector
  std::vector<double> dndtAse(hMesh.numberOfSamples, 0);
  std::vector<float>  phiAse(hMesh.numberOfSamples, 0);
  std::vector<double> mse(hMesh.numberOfSamples, 100000);
  std::vector<unsigned> totalRays(hMesh.numberOfSamples, 0);

  // Run Experiment
  std::vector<pthread_t> threadIds(maxGpus, 0);
  std::vector<float> runtimes(maxGpus, 0);
  switch(mode){
    // TODO: Replace completly by MPI
    case RAY_PROPAGATION_GPU:
      for(unsigned gpu_i = 0; gpu_i < maxGpus; ++gpu_i){
        const unsigned samplesPerNode = maxSampleRange-minSampleRange+1;
        const float samplePerGpu = samplesPerNode / (float) maxGpus;
        unsigned minSample_i = gpu_i * samplePerGpu;
        unsigned maxSample_i = min((float)samplesPerNode, (gpu_i + 1) * samplePerGpu);

        minSample_i += minSampleRange;
        maxSample_i += minSampleRange; 

        threadIds[gpu_i] = calcPhiAseThreaded( minRaysPerSample,
            maxRaysPerSample,
            maxRepetitions,
            dMesh.at(gpu_i),
            hMesh,
            sigmaAInterpolated,
            sigmaEInterpolated,
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
      usedGpus = calcPhiAseMPI( minRaysPerSample,
          maxRaysPerSample,
          maxRepetitions,
          dMesh.at(0),
          hMesh,
          sigmaAInterpolated,
          sigmaEInterpolated,
          mseThreshold,
          useReflections,
          phiAse,
          mse,
          totalRays,
          devices.at(0)
          );
      runmode = "RAY PROPAGATION MPI";
      break;

    case FOR_LOOPS: //Possibly deprecated!
      // TODO: make available for MPI?
      runtime = forLoopsClad( &dndtAse,
          minRaysPerSample,
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

    default:
      exit(0);
  }

  // Print Solution
  if(verbosity & V_DEBUG){
      for(unsigned sample_i = 0; sample_i < hMesh.numberOfSamples; ++sample_i){
        dndtAse.at(sample_i) = calcDndtAse(hMesh, maxSigmaA, maxSigmaE, phiAse.at(sample_i), sample_i);
        if(sample_i <=10)
          dout(V_DEBUG) << "Dndt ASE[" << sample_i << "]: " << dndtAse.at(sample_i) << " " << mse.at(sample_i) << std::endl;
      }
      for(unsigned sample_i = 0; sample_i < hMesh.numberOfSamples; ++sample_i){
        dout(V_DEBUG) << "PHI ASE[" << sample_i << "]: " << phiAse.at(sample_i) << " " << mse.at(sample_i) <<std::endl;
        if(sample_i >= 10) break;
      }
  }

  // Compare with vtk input
  // if(compareLocation!="") {
  //   std::vector<double> compareAse = compareVtk(dndtAse, compareLocation);

  // }

  // Write experiment data
  // output folder has to be the same as TMP_FOLDER in the calling MatLab script
  writeMatlabOutput(
      outputPath,
      phiAse,
      totalRays,
      mse,
      hMesh.numberOfSamples,
      hMesh.numberOfLevels
		    );


  // Write solution to vtk files
  if(writeVtk){
    std::vector<double> tmpPhiAse(phiAse.begin(), phiAse.end());
    std::vector<double> tmpTotalRays(totalRays.begin(), totalRays.end());

    writePointsToVtk(hMesh, dndtAse, outputPath + "vtk/dndt", minRaysPerSample, maxRaysPerSample, mseThreshold, useReflections, runtime);
    writePointsToVtk(hMesh, tmpPhiAse, outputPath + "vtk/phiase", minRaysPerSample, maxRaysPerSample, mseThreshold, useReflections, runtime);
    writePointsToVtk(hMesh, mse, outputPath + "vtk/mse", minRaysPerSample, maxRaysPerSample, mseThreshold, useReflections, runtime);
    writePointsToVtk(hMesh, tmpTotalRays, outputPath + "vtk/total_rays", minRaysPerSample, maxRaysPerSample, mseThreshold, useReflections, runtime);
  }

  //Print statistics
  if(verbosity & V_STAT){
    for(std::vector<double>::iterator it = mse.begin(); it != mse.end(); ++it){
      maxMSE = max(maxMSE, *it);
      avgMSE += *it;
      if(*it >= mseThreshold)
        highMSE++;
    }
    avgMSE /= mse.size();

    std::cout.imbue(std::locale(""));
    dout(V_STAT | V_NOLABEL) << std::endl;
    dout(V_STAT) << "=== Statistics ===" << std::endl;
    dout(V_STAT) << "Runmode           : " << runmode.c_str() << std::endl;
    dout(V_STAT) << "Prisms            : " << (int) hMesh.numberOfPrisms << std::endl;
    dout(V_STAT) << "Samples           : " << (int) dndtAse.size() << std::endl;
    dout(V_STAT) << "RaysPerSample     : " << minRaysPerSample;
    if(maxRaysPerSample > minRaysPerSample) { dout(V_STAT | V_NOLABEL) << " - " << maxRaysPerSample << " (adaptive)"; }
    dout(V_STAT | V_NOLABEL) << std::endl;
    dout(V_STAT) << "sum(totalRays)    : " << std::accumulate(totalRays.begin(), totalRays.end(), 0.) << std::endl;
    dout(V_STAT) << "MSE threshold     : " << mseThreshold << std::endl;
    dout(V_STAT) << "Wavelength        : " << sigmaA.size() << std::endl;
    dout(V_STAT) << "int. Wavelength   : " << sigmaAInterpolated.size() << std::endl;
    dout(V_STAT) << "max. MSE          : " << maxMSE << std::endl;
    dout(V_STAT) << "avg. MSE          : " << avgMSE << std::endl;
    dout(V_STAT) << "too high MSE      : " << highMSE << std::endl;
    dout(V_STAT) << "Nr of GPUs        : " << usedGpus << std::endl;
    dout(V_STAT) << "Runtime           : " << difftime(time(0),starttime) << "s" << std::endl;
    dout(V_STAT) << std::endl;
    if(maxRaysPerSample > minRaysPerSample){
      dout(V_STAT) << "=== Sampling resolution as Histogram ===" << std::endl;
      ray_histogram(totalRays,maxRaysPerSample,mseThreshold,mse);
    }
    dout(V_STAT) << std::endl;

  }
  return 0;

}
