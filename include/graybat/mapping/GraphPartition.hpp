#pragma once

#include <vector>    /* std::vector */
#include <metis.h>   /* idx_t, METIS_PartGraphKway */

namespace graybat {

    namespace mapping {

	/**
	 * Partitioning of the communication graph
	 * into k parts. k is set either to the 
	 * number of peers that want to take
	 * part in communication or is given
	 * as an input parameter.
	 *
	 * @param  nParts number of parts to partition
	 * @return set of vertices that belong "together"
	 */
	struct GraphPartition {

	    GraphPartition() : GraphPartition(0){

	    }
	    
	    GraphPartition(unsigned nParts)
		: nParts(nParts){

	    }

	    /**
	     * @Brief Translates the graph into the compressed row storage
	     *        format (CSR) which can be parsed by Metis.
	     *
	     * @See http://en.wikipedia.org/wiki/Sparse_matrix#Compressed_row_Storage_.28CRS_or_CSR.29
	     *
	     */
	    template<typename T_Graph>
	    std::pair<std::vector<idx_t>, std::vector<idx_t> > toCompressedRowStorage(T_Graph &graph) {
		
		typedef typename T_Graph::Vertex Vertex;
		typedef typename T_Graph::Edge   Edge;

		unsigned i = 0;

		std::vector<idx_t> xadj(1,i);
		std::vector<idx_t> adjncy;

		
		for(Vertex v : graph.getVertices()){
		    for(Edge link : graph.getOutEdges(v)){
			adjncy.push_back(link.target.id);
			i++;
			
		    }
		    xadj.push_back(i);

		}
		
		return std::make_pair(xadj, adjncy);
	    }

	    template<typename T_Graph>
	    std::vector<typename T_Graph::Vertex> operator()(const unsigned processID, const unsigned processCount, T_Graph &graph){

		typedef typename T_Graph::Vertex Vertex;
		std::vector<Vertex> myVertices;
		auto csr = toCompressedRowStorage(graph);

		if(nParts == 0){
		    nParts = processCount;
		}

		if(nParts == 1){
		    return graph.getVertices();
		}
		
		idx_t nVertices = graph.getVertices().size();
		idx_t nWeights  = 1;
		idx_t objval;
		std::vector<idx_t> part(nVertices, 0);

		METIS_PartGraphKway(&nVertices, &nWeights,
				    csr.first.data(), csr.second.data(),
				    NULL, NULL, NULL, &nParts, NULL,
				    NULL, NULL,
				    &objval,
				    part.data());

		
		for(unsigned part_i = 0; part_i < part.size(); part_i++){
		    if(part[part_i] == (int)processID){
			myVertices.push_back(graph.getVertex(part_i));
		    }
		    
		} 
		
		return myVertices;
	    }

	private:
	    idx_t nParts;

	};

    } /* mapping */
    
} /* graybat */
