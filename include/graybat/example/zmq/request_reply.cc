// GrayBat
#include <graybat.hpp>
#include <mapping/Roundrobin.hpp>
#include <pattern/StarBidirectional.hpp>
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
    typedef graybat::graphPolicy::BGL<>    GP;
    
    // Cage
    typedef graybat::Cage<CP, GP>   MyCage;
    typedef typename MyCage::Event  Event;
    typedef typename MyCage::Vertex Vertex;
    typedef typename MyCage::Edge   Edge;

    /***************************************************************************
     * Initialize Communication
     ****************************************************************************/
    // Create GoL Graph
    MyCage cage(graybat::pattern::StarBidirectional(4));
    
    // Distribute vertices
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
		    boost::optional<Edge> replyEdge = recvEdge.inverse();
		    if(replyEdge){
			cage.send(*replyEdge, world);
			std::cout << "Send " << world[0] << std::endl;			
		    }
		
		}
	    }

	    // Clients
	    if(v != server){
		// Send the hello request
		for(Edge sendEdge : cage.getOutEdges(v)){
		    cage.send(sendEdge, hello);
		    std::cout << "Send " << hello[0] << std::endl;
		
		    // Get the reply
		    boost::optional<Edge> replyEdge = sendEdge.inverse();
		    if(replyEdge){
			cage.recv(*replyEdge, world);
			std::cout << "Received " << world[0] << std::endl;
		    }

		}
	    }

	}
	
    }
    
    return 0;
}
