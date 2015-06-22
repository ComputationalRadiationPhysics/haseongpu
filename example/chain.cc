// GrayBat
#include <Cage.hpp>
#include <communicationPolicy/BMPI.hpp>
#include <graphPolicy/BGL.hpp>

// Mappings
#include <mapping/Consecutive.hpp>
#include <mapping/Random.hpp>
#include <mapping/Roundrobin.hpp>

// Pattern
#include <pattern/Chain.hpp>

// STL
#include <iostream>   /* std::cout */
#include <vector>     /* std::vector */
#include <array>      /* std::array */
#include <cmath>      /* sqrt */
#include <cstdlib>    /* atoi */
#include <numeric>    /* std::accumulate */


int exp() {
    /***************************************************************************
     * Configuration
     ****************************************************************************/

    // CommunicationPolicy
    typedef graybat::communicationPolicy::BMPI CP;
    
    // GraphPolicy
    typedef graybat::graphPolicy::BGL<>    GP;
    
    // Cage
    typedef graybat::Cage<CP, GP> Cage;
    typedef typename Cage::Event  Event;
    typedef typename Cage::Vertex Vertex;

    /***************************************************************************
     * Initialize Communication
     ****************************************************************************/
    const unsigned nChainLinks = 1000;
    auto inc = [](unsigned &a){a++;};

    
    // Create GoL Graph
    Cage cage;

    cage.setGraph(graybat::pattern::Chain(nChainLinks));
    
    // Distribute vertices
    cage.distribute(graybat::mapping::Consecutive());

    /***************************************************************************
     * Run Simulation
     ****************************************************************************/
    std::vector<Event> events;

    std::array<unsigned, 1> input {{0}};
    std::array<unsigned, 1> output {{0}};
    std::array<unsigned, 1> intermediate {{0}};

    const Vertex entry = cage.getVertex(0);
    const Vertex exit  = cage.getVertex(nChainLinks-1);


    for(Vertex v : cage.hostedVertices){

	if(v == entry){
	    v.spread(input, events);
	    std::cout << "Input: " << input[0] << std::endl;
	}

	if(v == exit){
	    v.collect(output);
	    std::cout << "Output: " << output[0] << std::endl;
	}

	if(v != entry and v != exit){
	    v.collect(intermediate);
	    inc(intermediate[0]);
	    std::cout << "Increment: " << intermediate[0] << std::endl;
	    v.spread(intermediate, events);
	    
	}
	
    }

    for(unsigned i = 0; i < events.size(); ++i){
	events.back().wait();
	events.pop_back();
    }
    
    return 0;

}

int main(){
    exp();
    return 0;
}
