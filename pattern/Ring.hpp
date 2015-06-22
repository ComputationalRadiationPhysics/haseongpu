/*******************************************************************************
 *
 * GRAPH TOPOLOGY GENERATORS
 *
 *******************************************************************************/

namespace graybat {

    namespace pattern {

	typedef unsigned                                                        VertexID;
	typedef std::pair<VertexID, VertexID>                                   EdgeDescription;
	typedef std::pair<std::vector<VertexID>, std::vector<EdgeDescription> > GraphDescription;

    
	struct Ring {

	    const unsigned verticesCount;

	    Ring(const unsigned verticesCount) :
		verticesCount(verticesCount) {

	    }
    

	    GraphDescription operator()(){
		std::vector<VertexID> vertices(verticesCount);
		std::vector<EdgeDescription> edges;

		if(vertices.size() != 1) {
		    for(unsigned i = 0; i < vertices.size(); ++i){
			if(i == vertices.size() - 1){
			    edges.push_back(std::make_pair(i, 0));
			}
			else {
			    edges.push_back(std::make_pair(i, i + 1));
			}
		
		    }
	    
		}
		return std::make_pair(vertices,edges);
	
	    }

	};

    } /* pattern */

} /* graybat */
