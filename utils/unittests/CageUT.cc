// boost 
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Cage
#include <boost/test/unit_test.hpp>

// graybat
#include <graybat.hpp>
#include <mapping/Random.hpp>
#include <mapping/Consecutive.hpp>
#include <pattern/FullyConnected.hpp>
#include <pattern/Star.hpp>
#include <pattern/Grid.hpp>


// stl
#include <array>
#include <vector>
#include <iostream>
#include <functional>

/***************************************************************************
 * Configuration
 ****************************************************************************/

// CommunicationPolicy
typedef graybat::communicationPolicy::BMPI CP;
    
// GraphPolicy
typedef graybat::graphPolicy::BGL<>   GP;
    
// Cage
typedef graybat::Cage<CP, GP>   MyCage;
typedef typename MyCage::Event  Event;
typedef typename MyCage::Vertex Vertex;
typedef typename MyCage::Edge   Edge;


/***************************************************************************
 * Test Cases
 ****************************************************************************/

BOOST_AUTO_TEST_SUITE(point_to_point)


MyCage allToAll(graybat::pattern::FullyConnected(2));
MyCage star(graybat::pattern::Star(2));

BOOST_AUTO_TEST_CASE( send_recv ){
    star.distribute(graybat::mapping::Consecutive());
    const unsigned nElements = 1000;
    
    std::vector<Event> events; 
    std::vector<unsigned> send(nElements,0);
    std::vector<unsigned> recv(nElements,0);

    for(unsigned i = 0; i < send.size();++i){
	send.at(i) = i;
    }

    // Send state to neighbor cells
    for(Vertex &v : star.hostedVertices){
	for(Edge edge : star.getOutEdges(v)){
	    star.send(edge, send);
	    
	}
    }

    // Recv state from neighbor cells
    for(Vertex &v : star.hostedVertices){
	for(Edge edge : star.getInEdges(v)){
	    star.recv(edge, recv);
	    for(unsigned i = 0; i < recv.size();++i){
		BOOST_CHECK_EQUAL(recv.at(i), i);
	    }

	}
	
    }
    
}


BOOST_AUTO_TEST_CASE( asyncSend_recv ){

    allToAll.distribute(graybat::mapping::Consecutive());

    const unsigned nElements = 1000;
    
    std::vector<Event> events; 
    std::vector<unsigned> send(nElements,0);
    std::vector<unsigned> recv(nElements,0);

    for(unsigned i = 0; i < send.size();++i){
	send.at(i) = i;
    }

    // Send state to neighbor cells
    for(Vertex &v : allToAll.hostedVertices){
	for(Edge edge : allToAll.getOutEdges(v)){
	    allToAll.send(edge, send, events);
	}
    }

    // Recv state from neighbor cells
    for(Vertex &v : allToAll.hostedVertices){
	for(Edge edge : allToAll.getInEdges(v)){
	    allToAll.recv(edge, recv);
	    for(unsigned i = 0; i < recv.size();++i){
		BOOST_CHECK_EQUAL(recv.at(i), i);
	    }

	}
	
    }
    
    // Wait to finish events
    for(unsigned i = 0; i < events.size(); ++i){
	events.back().wait();
	events.pop_back();
    }

}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( collectives )

MyCage grid(graybat::pattern::Grid(3,3));
MyCage star(graybat::pattern::Star(20));

BOOST_AUTO_TEST_CASE( reduce ){

    grid.distribute(graybat::mapping::Consecutive());

    const unsigned nElements = 10;
    
    std::vector<unsigned> send(nElements, 1);
    std::vector<unsigned> recv(nElements, 0);

    Vertex rootVertex = grid.getVertex(0);
    
    for(Vertex v: grid.hostedVertices){
	grid.reduce(rootVertex, v, std::plus<unsigned>(), send, recv);
    }

    if(grid.peerHostsVertex(rootVertex)){
	for(unsigned receivedElement: recv){
	    BOOST_CHECK_EQUAL(receivedElement, grid.getVertices().size());
	}
    }
    
}


BOOST_AUTO_TEST_CASE( allReduce ){

    grid.distribute(graybat::mapping::Consecutive());

    const unsigned nElements = 10;
    
    std::vector<unsigned> send(nElements, 1);
    std::vector<unsigned> recv(nElements, 0);

    for(Vertex v: grid.hostedVertices){
	grid.allReduce(v, std::plus<unsigned>(), send, recv);
    }

    for(unsigned receivedElement: recv){
	BOOST_CHECK_EQUAL(receivedElement, grid.getVertices().size());
    }
    
}

BOOST_AUTO_TEST_CASE( gather ){

    grid.distribute(graybat::mapping::Consecutive());

    const unsigned nElements = 10;
    const unsigned testValue = 1;
    const bool reorder = true;

    
    std::vector<unsigned> send(nElements, testValue);
    std::vector<unsigned> recv(nElements * grid.getVertices().size(), 0);

    Vertex rootVertex = grid.getVertex(0);
    
    for(Vertex v: grid.hostedVertices){
	grid.gather(rootVertex, v, send, recv, reorder);
    }

    if(grid.peerHostsVertex(rootVertex)){
	for(unsigned receivedElement: recv){
	    BOOST_CHECK_EQUAL(receivedElement, testValue);
	}
    }
    
    
}

BOOST_AUTO_TEST_CASE( allGather ){

    grid.distribute(graybat::mapping::Consecutive());

    const unsigned nElements = 10;
    const unsigned testValue = 1;
    const bool reorder = true;

    
    std::vector<unsigned> send(nElements, testValue);
    std::vector<unsigned> recv(nElements * grid.getVertices().size(), 0);

    for(Vertex v: grid.hostedVertices){
	grid.allGather(v, send, recv, reorder);
    }

    for(unsigned receivedElement: recv){
	BOOST_CHECK_EQUAL(receivedElement, testValue);
    }
        
}

BOOST_AUTO_TEST_CASE( spreadAndCollect ){

    star.distribute(graybat::mapping::Consecutive());

    const unsigned nElements = 10;
    const unsigned testValue = 1;
    std::vector<Event> events; 
    
    std::vector<unsigned> send(nElements, testValue);


    for(Vertex v: star.hostedVertices){
	v.spread(send, events);
    }

    for(Vertex v: star.hostedVertices){
	std::vector<unsigned> recv(v.nInEdges() * nElements, 0);
	v.collect(recv);
	for(unsigned receivedElement: recv){
	    BOOST_CHECK_EQUAL(receivedElement, testValue);
	}
	
    }

    for(unsigned i = 0; i < events.size(); ++i){
	events.back().wait();
	events.pop_back();
    }
  
}

BOOST_AUTO_TEST_SUITE_END()
