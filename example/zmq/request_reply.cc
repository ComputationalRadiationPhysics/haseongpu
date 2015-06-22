// GrayBat
<<<<<<< HEAD
#include <graybat.hpp>
#include <mapping/Roundrobin.hpp>
#include <pattern/StarBidirectional.hpp>
=======
#include <Cage.hpp>
#include <communicationPolicy/BMPI.hpp>
#include <graphPolicy/BGL.hpp>

#include <mapping/Roundrobin.hpp>
#include <pattern/BidirectionalStar.hpp>
>>>>>>> 49f4bd3... Adapted graybat to haseongpu so it can be used
#include <boost/optional.hpp> /* boost::optional */

/**
 * @brief Clients start communication with "Hello", Server answers with "World".
 *
 * This example implements the request and reply example of ZMQ[1]
 * in graybat style.
 *
 * [1] http://zguide.zeromq.org/page:all#Ask-and-Ye-Shall-Receive
 *
 */


int main() {
    /***************************************************************************
     * Configuration
     ****************************************************************************/

    // CommunicationPolicy
    typedef graybat::communicationPolicy::BMPI CP;
    
    // GraphPolicy
    typedef graybat::graphPolicy::BGL<>        GP;
    
    // Cage
    typedef graybat::Cage<CP, GP> Cage;
    typedef typename Cage::Event  Event;
    typedef typename Cage::Vertex Vertex;
    typedef typename Cage::Edge   Edge;

    /***************************************************************************
     * Initialize Communication
     ****************************************************************************/
    // Create GoL Graph
    Cage cage;

    // Set communication pattern
    cage.setGraph(graybat::pattern::BidirectionalStar(cage.getPeers().size()));
    
    // Map vertices to peers
    cage.distribute(graybat::mapping::Roundrobin());

    /***************************************************************************
     * Run 
     ****************************************************************************/
    Vertex server = cage.getVertex(0);

    std::array<std::string, 1> hello{{"Hello"}};
    std::array<std::string, 1> world{{"World"}};
    
    while(true){

	for(Vertex v: cage.hostedVertices) {

	    // Server
	    if(v == server){
		for(Edge recvEdge : cage.getInEdges(v)){
		    // Wait for next request from client
		    cage.recv(recvEdge, hello);
		    std::cout << "Received " << hello[0] << std::endl;

		    // Send reply back to client
		    cage.send(recvEdge.inverse(), world);
		    std::cout << "Send " << world[0] << std::endl;			
		}
		
	    }

	    // Clients
	    if(v != server){
		for(Edge sendEdge : cage.getOutEdges(v)){
		    // Send a hello
		    cage.send(sendEdge, hello);
		    std::cout << "Send " << hello[0] << std::endl;
		
		    // Get the reply
		    cage.recv(sendEdge.inverse(), world);
		    std::cout << "Received " << world[0] << std::endl;
		}
		
	    }

	}
	
    }
    
    return 0;
}
