/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @mainpage HASEonGPU - High performance Amplified Spontaneous EmissioN on GPU
 *
 * Project with HZDR for porting their ASE-code to a GPU cluster.
 *
 * @author Erik Zenker, Carlchristian Eckert, Marius Melzer
 */

// STL
#include <assert.h> /* assert */
#include <string> /* string */
#include <vector> /* vector */
#include <stdlib.h> /* atoi */
#include <pthread.h> /* pthread_t, pthread_join */
#include <algorithm> /* std::max */
#include <numeric> /* accumulate*/
#include <stdexcept>

// BOOST
#include <boost/filesystem.hpp> /* fs::path */
namespace fs = boost::filesystem;

// User header files
#include <calc_phi_ase_threaded.hpp>
#include <calc_phi_ase_mpi.hpp>
#include <calc_phi_ase_graybat.hpp>
#include <parser.hpp> /* DeviceMode, ParallelMode */
#include <write_to_vtk.hpp>
#include <write_matlab_output.hpp>
#include <for_loops_clad.hpp>
#include <cudachecks.hpp>
#include <mesh.hpp>
#include <cuda_utils.hpp> /* getFreeDevices */
#include <logging.hpp>
#include <ray_histogram.hpp>
#include <types.hpp>

// default without V_DEBUG
unsigned verbosity = V_ERROR | V_INFO | V_WARNING | V_PROGRESS | V_STAT; // extern through logging.hpp

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

    // Statistics data
    float runtime       = 0.0;
    double maxMSE       = 0;
    float  avgMSE       = 0;
    unsigned highMSE    = 0;
    time_t starttime    = time(0);

    // Simulation data
    ExperimentParameters experiment;
    ComputeParameters    compute;
    Result               result;
    std::vector<Mesh>    meshs;

    // Parse commandline and prepate all data structures
    parse(argc, argv, experiment, compute, meshs, result);

    // This declaration will hopefully deleted in future
    Mesh mesh = meshs.at(0);
    unsigned maxGpus = compute.devices.size();
    std::vector<float> runtimes(maxGpus, 0);
    unsigned usedGPUs = 0;
    std::vector<ComputeParameters> computes(maxGpus, compute);

    
    /***************************************************************************
     * COMPUTATIONS
     **************************************************************************/    
    
    if (compute.deviceMode == DeviceMode::CPU){
            
        // Possibly deprecated! (Definitely deprecated!)
        runtime = forLoopsClad( &(result.dndtAse),
                                experiment.minRaysPerSample,
                                &mesh,
                                mesh.betaCells,
                                mesh.nTot,
                                experiment.sigmaA.at(0),
                                experiment.sigmaE.at(0),
                                mesh.numberOfPoints,
                                mesh.numberOfTriangles,
                                mesh.numberOfLevels,
                                mesh.thickness,
                                mesh.crystalTFluo);

    }else if(compute.deviceMode == DeviceMode::GPU){
            
        if(compute.parallelMode == ParallelMode::THREADED){

            for(unsigned gpu_i = 0; gpu_i < maxGpus; ++gpu_i){
                const unsigned samplesPerNode = compute.maxSampleRange - compute.minSampleRange+1;
                const float samplePerGpu = samplesPerNode / (float) maxGpus;
                unsigned minSample_i = gpu_i * samplePerGpu;
                unsigned maxSample_i = std::min((float)samplesPerNode, (gpu_i + 1) * samplePerGpu);

                minSample_i += compute.minSampleRange;
                maxSample_i += compute.minSampleRange; 

                computes[gpu_i].gpu_i = gpu_i;

                calcPhiAseThreaded( experiment,
                        computes[gpu_i],
                        meshs[gpu_i],
                        result,
                        minSample_i,
                        maxSample_i,
                        runtimes.at(gpu_i));
            }

            joinAll();
            usedGPUs = maxGpus;
            runtime = *(std::max_element(runtimes.begin(),runtimes.end()));
            cudaDeviceReset();      
        }else if(compute.parallelMode == ParallelMode::MPI){
                
            usedGPUs = calcPhiAseMPI( experiment,
                                      compute,
                                      mesh,
                                      result );

        }else if(compute.parallelMode == ParallelMode::GRAYBAT){
                
            usedGPUs = calcPhiAseGrayBat( experiment,
                                          compute,
                                          mesh,
                                          result );

        }else{
            dout(V_ERROR) << "No valid parallel-mode for GPU!" << std::endl;
            exit(1);
        }

    }else{
        dout(V_ERROR) << "No valid device-mode!" << std::endl;
        exit(1);
    }



    /***************************************************************************
     * PRINT SOLUTION
     **************************************************************************/
    if(verbosity & V_DEBUG){
        for(unsigned sample_i = 0; sample_i < mesh.numberOfSamples; ++sample_i){
            result.dndtAse.at(sample_i) = calcDndtAse(mesh,
                                                      experiment.maxSigmaA,
                                                      experiment.maxSigmaE,
                                                      result.phiAse.at(sample_i), sample_i);
            if(sample_i <=10)
                dout(V_DEBUG) << "Dndt ASE[" << sample_i << "]: " << result.dndtAse.at(sample_i) << " " << result.mse.at(sample_i) << std::endl;
        }
        for(unsigned sample_i = 0; sample_i < mesh.numberOfSamples; ++sample_i){
            dout(V_DEBUG) << "PHI ASE[" << sample_i << "]: " << result.phiAse.at(sample_i) << " " << result.mse.at(sample_i) <<std::endl;
            if(sample_i >= 10) break;
        }
    }


    /***************************************************************************
     * WRITE MATLAB OUTPUT
     **************************************************************************/
    // output folder has to be the same as TMP_FOLDER in the calling MatLab script
    writeMatlabOutput(compute.outputPath,
                      result.phiAse,
                      result.totalRays,
                      result.mse,
                      mesh.numberOfSamples,
                      mesh.numberOfLevels);

    /***************************************************************************
     * WRITE VTK FILES
     **************************************************************************/
    if(compute.writeVtk){
        std::vector<double> tmpPhiAse(result.phiAse.begin(), result.phiAse.end());
        std::vector<double> tmpTotalRays(result.totalRays.begin(), result.totalRays.end());

        writePointsToVtk( mesh,
                          result.dndtAse,
                          compute.outputPath / "vtk/dndt",
                          experiment.minRaysPerSample,
                          experiment.maxRaysPerSample,
                          experiment.mseThreshold,
                          experiment.useReflections,
                          runtime );
      
        writePointsToVtk( mesh,
                          tmpPhiAse,
                          compute.outputPath / "vtk/phiase",
                          experiment.minRaysPerSample,
                          experiment.maxRaysPerSample,
                          experiment.mseThreshold,
                          experiment.useReflections,
                          runtime );
      
        writePointsToVtk( mesh,
                          result.mse,
                          compute.outputPath / "vtk/mse",
                          experiment.minRaysPerSample,
                          experiment.maxRaysPerSample,
                          experiment.mseThreshold,
                          experiment.useReflections,
                          runtime );
      
        writePointsToVtk( mesh,
                          tmpTotalRays,
                          compute.outputPath / "vtk/total_rays",
                          experiment.minRaysPerSample,
                          experiment.maxRaysPerSample,
                          experiment.mseThreshold,
                          experiment.useReflections,
                          runtime );
    }

    /***************************************************************************
     * PRINT STATISTICS
     **************************************************************************/
    if(verbosity & V_STAT){
        for(std::vector<double>::iterator it = result.mse.begin(); it != result.mse.end(); ++it){
            maxMSE = std::max(maxMSE, *it);
            avgMSE += *it;
            if(*it >= experiment.mseThreshold)
                highMSE++;
        }
        avgMSE /= result.mse.size();

        try{ std::cout.imbue(std::locale("")); }
        catch(std::runtime_error e){}

        dout(V_STAT | V_NOLABEL) << std::endl;
        dout(V_STAT) << "=== Statistics ===" << std::endl;
        dout(V_STAT) << "DeviceMode        : " << compute.deviceMode << std::endl;
        dout(V_STAT) << "ParallelMode      : " << compute.parallelMode << std::endl;
        dout(V_STAT) << "Prisms            : " << (int) mesh.numberOfPrisms << std::endl;
        dout(V_STAT) << "Samples           : " << (int) result.dndtAse.size() << std::endl;
        dout(V_STAT) << "RaysPerSample     : " << experiment.minRaysPerSample;
        if(experiment.maxRaysPerSample > experiment.minRaysPerSample) { dout(V_STAT | V_NOLABEL) << " - " << experiment.maxRaysPerSample << " (adaptive)"; }
        dout(V_STAT | V_NOLABEL) << std::endl;
        dout(V_STAT) << "sum(totalRays)    : " << std::accumulate(result.totalRays.begin(), result.totalRays.end(), 0.) << std::endl;
        dout(V_STAT) << "MSE threshold     : " << experiment.mseThreshold << std::endl;
        //dout(V_STAT) << "Wavelength        : " << experiment.sigmaA.size() << std::endl;
        dout(V_STAT) << "int. Wavelength   : " << experiment.sigmaA.size() << std::endl;
        dout(V_STAT) << "max. MSE          : " << maxMSE << std::endl;
        dout(V_STAT) << "avg. MSE          : " << avgMSE << std::endl;
        dout(V_STAT) << "too high MSE      : " << highMSE << std::endl;
        dout(V_STAT) << "Nr of GPUs        : " << usedGPUs << std::endl;
        dout(V_STAT) << "Runtime           : " << difftime(time(0), starttime) << "s" << std::endl;
        dout(V_STAT) << std::endl;
        if(experiment.maxRaysPerSample > experiment.minRaysPerSample){
            dout(V_STAT) << "=== Sampling resolution as Histogram ===" << std::endl;
            ray_histogram(result.totalRays, experiment.maxRaysPerSample, experiment.mseThreshold, result.mse);
        }
        dout(V_STAT) << std::endl;

    }
  
    return 0;

}
