#pragma once

template <class T_Cage>
struct CommunicationVertex {

    typedef unsigned                             VertexID; 
    typedef T_Cage                               Cage;
    typedef typename Cage::GraphPolicy           GraphPolicy;
    typedef typename Cage::Edge                  Edge;
    typedef typename Cage::Event                 Event;
    typedef typename GraphPolicy::VertexProperty VertexProperty;

    VertexID id;
    VertexProperty &vertexProperty;
    Cage &cage;

    CommunicationVertex(const VertexID id, VertexProperty &vertexProperty, Cage &cage) :
	id(id),
	vertexProperty(vertexProperty),
	cage(cage){
	    
    }

    /***************************************************************************
     * Graph Operations
     ****************************************************************************/
    VertexProperty& operator()(){
	return vertexProperty;
    }

    CommunicationVertex& operator=(const CommunicationVertex &other){
	id = other.id;
	vertexProperty = other.vertexProperty;

	return *this;
    }

    
    size_t nInEdges() const {
	return cage.getInEdges(*this).size();
    }

    size_t nOutEdges() const {
	return cage.getOutEdges(*this).size();
    }

    bool operator==(CommunicationVertex v){
	return (id == v.id);
    }

    bool operator!=(CommunicationVertex v){
	return (id != v.id);
    }

    /***************************************************************************
     * Communication Operations
     ****************************************************************************/

    template <typename T>
    void spread(const T& data, std::vector<Event> &events){
	cage.spread(*this, data, events);
    }

    template <typename T>
    void spread(const T& data){
	cage.spread(*this, data);
    }
    

    template <typename T>
    void collect(T& data){
	cage.collect(*this, data);
	
    }

    template <typename T_Op>
    typename T_Op::result_type accumulate(const T_Op op, const typename T_Op::result_type init){
	std::vector<typename T_Op::result_type> data (nInEdges());
	cage.collect(*this, data);
	return std::accumulate(data.begin(), data.end(), init, op);
	
    }
	
};
