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

    template <class T_Data>
    void spread(const T_Data& data, std::vector<Event> &events){
	cage.spread(*this, data, events);
    }

    template <class T_Data>
    void spread(const T_Data& data){
	cage.spread(*this, data);
    }
    

    template <class T_Data>
    void collect(T_Data& data){
	cage.collect(*this, data);
	
    }

    template <class T_Data, class T_Functor>
    void forward(T_Data& data, T_Functor f){
	cage.collect(*this, data);
	f(data);
	cage.spread(*this, data);
	
    }

    template <class T_Data, class T_Functor>
    void forward(T_Data& data){
	cage.collect(*this, data);
	cage.spread(*this, data);
	
    }    
    

    template <typename T_Op>
    typename T_Op::result_type accumulate(const T_Op op, const typename T_Op::result_type init){
	std::vector<typename T_Op::result_type> data (nInEdges());
	cage.collect(*this, data);
	return std::accumulate(data.begin(), data.end(), init, op);
	
    }
	
};
