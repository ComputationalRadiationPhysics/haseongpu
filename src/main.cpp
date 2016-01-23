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
#include <string>    /* string */
#include <vector>    /* vector */
#include <algorithm> /* std::max */
#include <numeric>   /* accumulate*/
#include <stdexcept>
#include <ctime>     /* time */
#include <locale>    /* std::locale */

// HASEonGPU
#include <parser.hpp> /* parse, DeviceMode, ParallelMode */
#include <write_matlab_output.hpp>
#include <logging.hpp>
#include <ray_histogram.hpp>
#include <types.hpp> /* ComputeParamets, ExperimentParameters, Result */
#include <calc_phi_ase_threaded.hpp>

#if defined(MPI_FOUND)
#include <calc_phi_ase_mpi.hpp>
#endif

#if defined(BOOST_MPI_FOUND) || defined(ZMQ_FOUND)
#include <calc_phi_ase_graybat.hpp>
#endif


// default without V_DEBUG
unsigned verbosity = V_ERROR | V_INFO | V_WARNING | V_PROGRESS | V_STAT; // extern through logging.hpp

int main(int argc, char **argv){

    // Statistics data
    double maxMSE    = 0;
    float  avgMSE    = 0;
    unsigned highMSE = 0;
    time_t starttime = time(0);

    // Simulation data
    ExperimentParameters experiment;
    ComputeParameters    compute;
    Result               result;

    // Parse commandline
    parse(argc, argv, experiment, compute, result);

    unsigned usedDevices = 0;
    
    /***************************************************************************
     * COMPUTATIONS
     **************************************************************************/    
    if(compute.parallelMode == ParallelMode::THREADED){

	usedDevices = calcPhiAseThreaded( experiment,
					  compute,
					  result );
    }
#if defined(MPI_FOUND)
    else if(compute.parallelMode == ParallelMode::MPI){
                
	usedDevices = calcPhiAseMPI( experiment,
				     compute,
				     result );

    }
#endif
#if defined(BOOST_MPI_FOUND) || defined(ZMQ_FOUND)
    else if(compute.parallelMode == ParallelMode::GRAYBAT){
	usedDevices = calcPhiAseGrayBat( experiment,
					 compute,
					 result );

    }
#endif
    else{
	dout(V_ERROR) << "No valid parallel-mode for GPU!" << std::endl;
	exit(1);
    }

    /***************************************************************************
     * WRITE MATLAB OUTPUT
     **************************************************************************/
    // output folder has to be the same as TMP_FOLDER in the calling MatLab script
    writeMatlabOutput(compute.outputPath,
                      result.phiAse,
                      result.totalRays,
                      result.mse,
                      experiment.numberOfSamples,
                      experiment.numberOfLevels);

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
        dout(V_STAT) << "Prisms            : " << experiment.numberOfPrisms << std::endl;
        dout(V_STAT) << "Samples           : " << (int) result.dndtAse.size() << std::endl;
        dout(V_STAT) << "RaysPerSample     : " << experiment.minRaysPerSample;
        if(experiment.maxRaysPerSample > experiment.minRaysPerSample) { dout(V_STAT | V_NOLABEL) << " - " << experiment.maxRaysPerSample << " (adaptive)"; }
        dout(V_STAT | V_NOLABEL) << std::endl;
        dout(V_STAT) << "sum(totalRays)    : " << std::accumulate(result.totalRays.begin(), result.totalRays.end(), 0.) << std::endl;
        dout(V_STAT) << "MSE threshold     : " << experiment.mseThreshold << std::endl;
        dout(V_STAT) << "int. Wavelength   : " << experiment.sigmaA.size() << std::endl;
        dout(V_STAT) << "max. MSE          : " << maxMSE << std::endl;
        dout(V_STAT) << "avg. MSE          : " << avgMSE << std::endl;
        dout(V_STAT) << "too high MSE      : " << highMSE << std::endl;
        dout(V_STAT) << "Nr of Devices     : " << usedDevices << std::endl;
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
