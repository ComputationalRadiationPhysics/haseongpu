#pragma once

#include <vector>    /* std::vector */
#include <stdlib.h>  /* srand, rand */

namespace graybat {

    namespace mapping {

	/**
	 * Random distribution of vertices of the *graph* to the
	 * the peers. All peers need to set the same random *seed*.
	 * Thus, all peers have the same random base. Therefore
	 * seeds that depend on varying parameters like time or
	 * pid or not applicable here.
	 *
	 * @param  seed static random seed for all peers
	 * @return random set of vertices of the *graph*
	 *
	 */
	struct Random {

	    Random(int seed)
		: seed(seed){

	    }
	  
	    template<typename T_Graph>
	    std::vector<typename T_Graph::Vertex> operator()(const unsigned processID, const unsigned processCount, T_Graph &graph){


		typedef typename T_Graph::Vertex Vertex;
		srand(seed);
		std::vector<Vertex> myVertices;
		unsigned vertexCount   = graph.getVertices().size();
		

		if(processID <= vertexCount){
		    for(Vertex v: graph.getVertices()){
			unsigned randomID = rand() % processCount;
			if(randomID == processID){
			    myVertices.push_back(v);
			}

		    }
		}
	      

		return myVertices;
	    }
	  
	private:
	    int seed;

	};

    } /* mapping */
    
} /* graybat */
