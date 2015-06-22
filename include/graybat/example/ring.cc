// GrayBat
#include <Cage.hpp>
#include <communicationPolicy/BMPI.hpp>
#include <graphPolicy/BGL.hpp>


// Mappings
#include <mapping/Consecutive.hpp>
#include <mapping/Random.hpp>
#include <mapping/Roundrobin.hpp>


// Pattern
#include <pattern/Ring.hpp>


// STL
#include <iostream>   /* std::cout */
#include <vector>     /* std::vector */
#include <array>      /* std::array */


struct Function {

    void process(std::tuple<unsigned, std::string> &a){
	std::get<0>(a)++;
	std::get<1>(a) += " world";
	
    }
    
};


int exp() {
    /***************************************************************************
     * Configuration
     ****************************************************************************/

    // CommunicationPolicy
    typedef graybat::communicationPolicy::BMPI CP;
    
    // GraphPolicy
    typedef graybat::graphPolicy::BGL<Function>    GP;
    
    // Cage
    typedef graybat::Cage<CP, GP> Cage;
    typedef typename Cage::Event  Event;
    typedef typename Cage::Vertex Vertex;

    /***************************************************************************
     * Initialize Communication
     ****************************************************************************/
    const unsigned nRingLinks = 50;

    // Create GoL Graph
    Cage cage;
    assert(cage.getPeers().size() >= nRingLinks);

    // Create ring communication pattern
    cage.setGraph(graybat::pattern::Ring(nRingLinks));

    
    // Distribute vertices
    cage.distribute(graybat::mapping::Roundrobin());

    
    /***************************************************************************
     * Run Simulation
     ****************************************************************************/
    std::vector<Event> events;

    std::array<std::tuple<unsigned, std::string>, 1> input{{std::make_tuple(0, "hello")}};
    std::array<std::tuple<unsigned, std::string>, 1> output;
    std::array<std::tuple<unsigned, std::string>, 1> intermediate;

    const Vertex stimula = cage.getVertex(0);

    for(Vertex v : cage.hostedVertices){

    	// Entry Vertex
    	if(v == stimula){
    	    v.spread(input, events);
    	    std::cout << "Input: " << std::get<0>(input[0]) << std::endl;	    	    
    	}

	while(true){
	    v.collect(intermediate);
	    v().process(intermediate[0]);
	    std::cout << "Increment: " << std::get<0>(intermediate[0]) << std::endl;
	    v.spread(intermediate);
		
	}
	    
    }
    
    return 0;

}

int main(){
    exp();
    return 0;
}
