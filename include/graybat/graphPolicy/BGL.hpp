#pragma once
#include <vector>
#include <iostream>
#include <tuple>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graphviz.hpp>
#include <tuple>
#include <string>
#include <fstream> /* std::fstream */
#include <utility> /* std::pair, std::make_pair */



namespace graybat {
    
    namespace graphPolicy {
	
	struct SimpleProperty{
	    typedef unsigned ID;
	    SimpleProperty() : id(0){}
	    SimpleProperty(ID id) : id(id){}

	    ID id;
	};

	/************************************************************************//**
         * @class BGL
	 *									   
	 * @brief A class to describe directed graphs.
	 *
	 * GraphPolicy on basis of the boost graph library.
	 *
	 ***************************************************************************/
	template <class T_VertexProperty = SimpleProperty, class T_EdgeProperty = SimpleProperty>
	class BGL {

	public:
	    // Public typedefs
	    typedef T_VertexProperty                                                VertexProperty;
	    typedef T_EdgeProperty                                                  EdgeProperty;	    
	    typedef std::pair<unsigned, unsigned>                                   EdgeDescription;
	    typedef std::pair<std::vector<unsigned>, std::vector<EdgeDescription> > GraphDescription;

	    typedef unsigned                                                        GraphID;


	    // BGL typdefs
	    typedef boost::adjacency_list<boost::vecS, 
					  boost::vecS, 
					  boost::bidirectionalS, 
					  boost::property<boost::vertex_index_t, size_t, std::pair<unsigned, VertexProperty> >,
					  boost::property<boost::edge_index_t, size_t, std::pair<unsigned, EdgeProperty> > > GraphType;

	    typedef boost::subgraph<GraphType> BGLGraph;
	    typedef typename BGLGraph::vertex_descriptor VertexID;
	    typedef typename BGLGraph::edge_descriptor   EdgeID;

	    typedef typename boost::graph_traits<BGLGraph>::in_edge_iterator   InEdgeIter;
	    typedef typename boost::graph_traits<BGLGraph>::out_edge_iterator  OutEdgeIter;
	    typedef typename boost::graph_traits<BGLGraph>::adjacency_iterator AdjacentVertexIter;
	    typedef typename boost::graph_traits<BGLGraph>::vertex_iterator    AllVertexIter;
	    


	    // Member
	    BGLGraph* graph;
	    std::vector<BGL<VertexProperty, EdgeProperty>> subGraphs;

	public: 
	    GraphID id;


	    /**
	     * @brief The graph has to be described by *edges* 
	     * (source Vertex ==> target Vertex) and
	     * the *vertices* of this graph.
	     *
	     */
	    BGL(GraphDescription graphDesc) :
		id(0){

		std::vector<unsigned> vertices     = graphDesc.first;
		std::vector<EdgeDescription> edges = graphDesc.second;

		graph = new BGLGraph(vertices.size());

		unsigned edgeCount = 0;

		for(EdgeDescription edge: edges){
		    VertexID srcVertex    = std::get<0>(edge);
		    VertexID targetVertex = std::get<1>(edge);
		    EdgeID edgeID = boost::add_edge(srcVertex, targetVertex, (*graph)).first;
		    setEdgeProperty(edgeID, std::make_pair(edgeCount++, EdgeProperty()));
		}

		// Bind vertex_descriptor and VertexProperty;
		for(unsigned vertexID = 0; vertexID < vertices.size(); ++vertexID){
		    setVertexProperty(vertexID, std::make_pair(vertexID, VertexProperty()));
		}
		
	    }	    

	    ~BGL(){

	    }

  
	    /**
	     * @brief Returns all vertices of the graph
	     * 
	     */
	    std::pair<AllVertexIter, AllVertexIter> getVertices(){
		return boost::vertices((*graph));

	    }

	    /**
	     * @brief Returns the edge between source and target vertex.
	     * 
	     */
	    std::pair<EdgeID, bool> getEdge(const VertexID source, const VertexID target){
		return boost::edge(source, target, *graph);
	    }

	    /**
	     * @brief Returns all vertices, that are adjacent (connected) to *vertex*
	     *
	     */
	    std::pair<AdjacentVertexIter, AdjacentVertexIter>  getAdjacentVertices(const VertexID id){
		return boost::adjacent_vertices(id, *graph);
	    }

	    /**
	     * @brief Returns all outgoing edges of *srcVertex* paired with its target vertex.
	     *
	     */
	    std::pair<OutEdgeIter, OutEdgeIter> getOutEdges(const VertexID id){
		return boost::out_edges((*graph).global_to_local(id), (*graph));
	    }

	    /**
	     * @brief Returns all incoming edges to *targetVertex* paired with its source vertex.
	     *
	     */
	    std::pair<InEdgeIter, InEdgeIter> getInEdges(const VertexID id){
		return boost::in_edges((*graph).global_to_local(id), (*graph));
	    }

	    /**
	     * @brief Returns the local id of *vertex* in this graph.
	     *
	     * If this graph has no supergraph (hasSuperGraph()==false) then local ids are the same as global ids.
	     */
	    VertexID getLocalID(VertexProperty vertex){
		return (*graph).global_to_local(vertex.id);
	    }
	    
	    void setVertexProperty(VertexID vertex, VertexProperty value){
		std::pair<unsigned, VertexProperty> propPair = (*graph)[vertex];
		(*graph)[vertex] = std::make_pair<propPair.first, value>;
	    }

	    void setVertexProperty(VertexID vertex, std::pair<unsigned, VertexProperty> propPair){
		(*graph)[vertex] = propPair;
	    }

	    
	    void setEdgeProperty(EdgeID edge, EdgeProperty value){
		std::pair<unsigned, EdgeProperty> propPair = (*graph)[edge];
		(*graph)[edge] = std::make_pair<propPair.first, value>;
	    }

	    void setEdgeProperty(EdgeID edge, std::pair<unsigned, EdgeProperty> propPair){
		(*graph)[edge] = propPair;
	    }	    

	    /**
	     * @brief Returns the property of *vertex*.
	     *
	     */
	    std::pair<unsigned, VertexProperty>& getVertexProperty(const VertexID vertex){
		return (*graph)[vertex];
	    }
	    
	    
	    /**
	     * @brief Return the property of *edge*.
	     *
	     */
	    std::pair<unsigned, EdgeProperty>& getEdgeProperty(const EdgeID edge){
		return (*graph)[edge];
	    }

	    /**
	     * @brief Return the vertex to which *edge* points to.
	     *
	     */
	    VertexID getEdgeTarget(EdgeID edge){
		return boost::target(edge, (*graph));
	    }

	    /**
	     * @brief Return the vertex to which *edge* points from.
	     *
	     */
	    VertexID getEdgeSource(EdgeID edge){
		return boost::source(edge, (*graph));
	    }
	      
	};

    } // namespace graphPolicy

} // namespace graybat
