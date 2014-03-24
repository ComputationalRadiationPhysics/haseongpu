// Libraries
#include <assert.h> /* assert */
#include <string> /* string */
#include <vector> /* vector */
#include <stdlib.h> /* atoi */
#include <pthread.h> /* pthread_t, pthread_join */
#include <algorithm> /* max_element */
#include <numeric> /* accumulate*/
#include <stdexcept>

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
#include <cuda_utils.h> /* getFreeDevices */
#include <logging.h>
#include <ray_histogram.h>
#include <interpolation.h> /* interpolateWavelength*/

// default without V_DEBUG
unsigned verbosity = V_ERROR | V_INFO | V_WARNING | V_PROGRESS | V_STAT; // extern through logging.h

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
  return gain_local * phiAse / mesh.crystalTFluo;
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
  int minSampleRange = -1;
  int maxSampleRange = -1;
  time_t starttime   = time(0);
  unsigned usedGpus  = 0;

  std::string inputPath;
  std::string outputPath;
  double mseThreshold = 0;

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
  if(fileToVector(inputPath + "sigmaA.txt", &sigmaA)) return 1;
  if(fileToVector(inputPath + "sigmaE.txt", &sigmaE)) return 1;
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
  std::vector<Mesh> meshs = parseMesh(inputPath, devices, maxGpus);

  checkSampleRange(&minSampleRange,&maxSampleRange,meshs[0].numberOfSamples);

  // Solution vector
  std::vector<double> dndtAse(meshs[0].numberOfSamples, 0);
  std::vector<float>  phiAse(meshs[0].numberOfSamples, 0);
  std::vector<double> mse(meshs[0].numberOfSamples, 100000);
  std::vector<unsigned> totalRays(meshs[0].numberOfSamples, 0);

  // Run Experiment
  std::vector<pthread_t> threadIds(maxGpus, 0);
  std::vector<float> runtimes(maxGpus, 0);
  switch(mode){
    // TODO: Replace completly by MPI
    case GPU_THREADED:
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
            meshs[gpu_i],
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

    case GPU_MPI:
      usedGpus = calcPhiAseMPI( minRaysPerSample,
          maxRaysPerSample,
          maxRepetitions,
          meshs[0],
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

    case CPU: //Possibly deprecated!
      // TODO: make available for MPI?
      runtime = forLoopsClad( &dndtAse,
          minRaysPerSample,
          &meshs[0],
          meshs[0].betaCells,
          meshs[0].nTot,
          sigmaA.at(0),
          sigmaE.at(0),
          meshs[0].numberOfPoints,
          meshs[0].numberOfTriangles,
          meshs[0].numberOfLevels,
          meshs[0].thickness,
          meshs[0].crystalTFluo);
      runmode = "For Loops";
      break;

    default:
      dout(V_ERROR) << "No valid runmode!" << std::endl;
      exit(0);
  }

  // Print Solution
  if(verbosity & V_DEBUG){
      for(unsigned sample_i = 0; sample_i < meshs[0].numberOfSamples; ++sample_i){
        dndtAse.at(sample_i) = calcDndtAse(meshs[0], maxSigmaA, maxSigmaE, phiAse.at(sample_i), sample_i);
        if(sample_i <=10)
          dout(V_DEBUG) << "Dndt ASE[" << sample_i << "]: " << dndtAse.at(sample_i) << " " << mse.at(sample_i) << std::endl;
      }
      for(unsigned sample_i = 0; sample_i < meshs[0].numberOfSamples; ++sample_i){
        dout(V_DEBUG) << "PHI ASE[" << sample_i << "]: " << phiAse.at(sample_i) << " " << mse.at(sample_i) <<std::endl;
        if(sample_i >= 10) break;
      }
  }

  // Write experiment data
  // output folder has to be the same as TMP_FOLDER in the calling MatLab script
  writeMatlabOutput(outputPath,
		    phiAse,
		    totalRays,
		    mse,
		    meshs[0].numberOfSamples,
		    meshs[0].numberOfLevels);

  // Write solution to vtk files
  if(writeVtk){
    std::vector<double> tmpPhiAse(phiAse.begin(), phiAse.end());
    std::vector<double> tmpTotalRays(totalRays.begin(), totalRays.end());

    writePointsToVtk(meshs[0], dndtAse, outputPath + "vtk/dndt", minRaysPerSample, maxRaysPerSample, mseThreshold, useReflections, runtime);
    writePointsToVtk(meshs[0], tmpPhiAse, outputPath + "vtk/phiase", minRaysPerSample, maxRaysPerSample, mseThreshold, useReflections, runtime);
    writePointsToVtk(meshs[0], mse, outputPath + "vtk/mse", minRaysPerSample, maxRaysPerSample, mseThreshold, useReflections, runtime);
    writePointsToVtk(meshs[0], tmpTotalRays, outputPath + "vtk/total_rays", minRaysPerSample, maxRaysPerSample, mseThreshold, useReflections, runtime);
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

    try{ std::cout.imbue(std::locale("")); }
    catch(std::runtime_error e){}

    dout(V_STAT | V_NOLABEL) << std::endl;
    dout(V_STAT) << "=== Statistics ===" << std::endl;
    dout(V_STAT) << "Runmode           : " << runmode.c_str() << std::endl;
    dout(V_STAT) << "Prisms            : " << (int) meshs[0].numberOfPrisms << std::endl;
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
