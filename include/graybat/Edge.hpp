#pragma once

#include <boost/optional.hpp> /* boost::optional */

template <class T_Cage>
struct CommunicationEdge {

    typedef T_Cage                               Cage;
    typedef unsigned                             EdgeID;
    typedef typename Cage::GraphPolicy           GraphPolicy;
    typedef typename Cage::Vertex                Vertex;
    typedef typename Cage::Event                 Event;
    typedef typename GraphPolicy::EdgeProperty   EdgeProperty;
    typedef typename GraphPolicy::VertexProperty VertexProperty;
	    
    EdgeID id;
    Vertex target;
    Vertex source;
    EdgeProperty &edgeProperty;
    Cage &cage;

    CommunicationEdge(const EdgeID id,
		      Vertex source,
		      Vertex target,
		      EdgeProperty &edgeProperty,
		      Cage &cage) :
	id(id),
	target(target),
	source(source),	
	edgeProperty(edgeProperty),
	cage(cage){
	    
    }

    /***************************************************************************
     * Graph Operations
     ****************************************************************************/

    EdgeProperty& operator()(){
	return edgeProperty;
    }

     CommunicationEdge inverse(){
	 boost::optional<CommunicationEdge> edge = cage.getEdge(target, source);
	 if(edge){
	     return *edge;
	 }
	 return *this;
    }

    /***************************************************************************
     * Communication Operations
     ****************************************************************************/
    
    template <class T_Send>
    Event operator<<(const T_Send &data){
	std::vector<Event> events;
	cage.send(*this, data, events);
	return events.back();
    }

    template <class T_Recv>
    void operator>>(T_Recv &data){
	cage.recv(*this, data);
    }
	
};
