#pragma once

#include <algorithm> /* std::min */
#include <vector>    /* std::vector */
#include <assert.h>  /* assert */

namespace graybat {
    
    namespace mapping {
    
	struct Consecutive {

	    template<typename T_Graph>
	    std::vector<typename T_Graph::Vertex> operator()(const unsigned processID, const unsigned processCount, T_Graph &graph){

		typedef typename T_Graph::Vertex Vertex;

		unsigned vertexCount      = graph.getVertices().size();
		unsigned vertexPerProcess = ceil((float)vertexCount / processCount);

		// More processes than vertices
		if(processID > vertexCount - 1){
		    return std::vector<Vertex>();
		}

		unsigned minVertex = processID * vertexPerProcess;
		unsigned maxVertex = minVertex + vertexPerProcess;

		// Slice maxVertex of last process
		if(minVertex > vertexCount){
		    return std::vector<Vertex>();
		}
	    
		maxVertex = std::min(maxVertex, vertexCount);
	
		assert(minVertex <= maxVertex);
	
		std::vector<Vertex> vertices = graph.getVertices();
		std::vector<Vertex> myVertices(vertices.begin() + minVertex, vertices.begin() + maxVertex);
		return myVertices;

	    }

	};

    } /* mapping */
    
} /* graybat */
