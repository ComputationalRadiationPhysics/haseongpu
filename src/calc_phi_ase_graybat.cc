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

// STL
#include <vector>
#include <algorithm> /* std::iota */


// HASEonGPU
#include <calc_phi_ase.hpp>
#include <mesh.hpp>
#include <logging.hpp>
#include <progressbar.hpp>
#include <types.hpp>

// GrayBat
#include <graybat/Cage.hpp>
#include <graybat/communicationPolicy/BMPI.hpp>
#include <graybat/graphPolicy/BGL.hpp>
#include <graybat/pattern/BidirectionalStar.hpp>
#include <graybat/mapping/Roundrobin.hpp>

// Const messages
const int abortTag   = -1;
const float requestTag = -1.0;
const std::array<float, 4> requestMsg{{requestTag, 0, 0, 0}};
const std::array<int, 1> abortMsg{{abortTag}};

template <class Vertex, class Cage>
void distributeSamples(Vertex master,
                       const std::vector<unsigned> samples,
                       Cage &cage,
                       const Mesh& mesh,
                       Result& result ){

    typedef typename Cage::Edge Edge;
    
    // Messages
    std::array<float, 4> resultMsg;
    std::array<int, 1>   sampleMsg; 
    unsigned nReceivedResults = 0;
    
    auto sample = samples.begin();
    while(nReceivedResults != samples.size()){
        // Receive request or results
        Edge inEdge = cage.recv(resultMsg);
		
        if(resultMsg[0] == requestTag){
            
            if(sample != samples.end()){
                // Send next sample
                sampleMsg = std::array<int, 1>{{ (int) *sample++ }};
                cage.send(inEdge.inverse(), sampleMsg);
                
            }

        }
        else {
            // Process result
            unsigned sample_i             = (unsigned) (resultMsg[0]);
            result.phiAse.at(sample_i)    = resultMsg[1];
            result.mse.at(sample_i)       = resultMsg[2];
            result.totalRays.at(sample_i) = (unsigned) resultMsg[3];

            // Update receive counter
            nReceivedResults++;
            
            // Update progress bar
            fancyProgressBar(mesh.numberOfSamples);

        }

    }

    // Send abort message to all vertices
    master.spread(abortMsg);
       
}

template <class Vertex, class Cage>
void processSamples(const Vertex slave,
		   const Vertex master,
		   Cage &cage,
		   const ExperimentParameters& experiment,
		   const ComputeParameters& compute,
		   const Mesh& mesh,
		   Result& result ){

    typedef typename Cage::Edge Edge;
    
    verbosity &= ~V_PROGRESS;

    // Messages
    std::array<float, 4> resultMsg;
    std::array<int, 1>   sampleMsg;

    float runtime = 0.0;

    bool abort = false;
    
    while(!abort){

	// Get edge to master
	Edge outEdge = cage.getEdge(slave, master);

	// Request new sampling point
	cage.send(outEdge, requestMsg);

	// Receive new sampling point or abort
	cage.recv(outEdge.inverse(), sampleMsg);
		    
	if(sampleMsg.at(0) == abortTag){
	  abort = true;
	}
	else {
	    calcPhiAse ( experiment,
			 compute,
			 mesh,
			 result,
			 sampleMsg.at(0),
			 sampleMsg.at(0) + 1,
			 runtime);
			
	    unsigned sample_i = sampleMsg[0];
	    resultMsg = std::array<float, 4>{{ (float) sample_i,
					       (float) result.phiAse.at(sample_i),
					       (float) result.mse.at(sample_i),
					       (float) result.totalRays.at(sample_i) }};

	    // Send simulation results
	    cage.send(outEdge, resultMsg);
						
	}

    }

 
}


float calcPhiAseGrayBat ( const ExperimentParameters &experiment,
			  const ComputeParameters &compute,
			  Mesh& mesh,
			  Result &result ){

    /***************************************************************************
     * CAGE
     **************************************************************************/
    // Configuration

    typedef typename graybat::communicationPolicy::BMPI CP;
    typedef typename graybat::graphPolicy::BGL<>        GP;
    typedef typename graybat::Cage<CP, GP>              Cage;
    typedef typename Cage::Vertex                       Vertex;
    typedef typename Cage::Edge                         Edge;
    
    // Init
    Cage cage;
    const unsigned nPeers = cage.getPeers().size();
    cage.setGraph(graybat::pattern::BidirectionalStar(nPeers));
    cage.distribute(graybat::mapping::Roundrobin());
    const Vertex master = cage.getVertex(0);
    
    /***************************************************************************
     * ASE SIMULATION
     **************************************************************************/
    // Create sample indices
    std::vector<unsigned> samples(mesh.numberOfSamples);
    std::iota(samples.begin(), samples.end(), 0);

    // Determine phi ase for each sample
    for(Vertex vertex : cage.hostedVertices) {

	/*******************************************************************
	 * MASTER
	 *******************************************************************/
	if(vertex == master){
	    distributeSamples(vertex, samples, cage, mesh, result);

	}

	/*******************************************************************
	 * SLAVES
	 *******************************************************************/
	if(vertex != master){
	  processSamples(vertex, master, cage, experiment, compute, mesh, result);
	  cage.~Cage();
          mesh.free();
	  exit(0);

	}	

    }
    
    return cage.getVertices().size() - 1;
    
}
