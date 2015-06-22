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

    
    struct Chain {

      const unsigned verticesCount;

      Chain(const unsigned verticesCount) :
	verticesCount(verticesCount) {

      }
    

      GraphDescription operator()(){
	std::vector<VertexID> vertices(verticesCount);
    	std::vector<EdgeDescription> edges;

	if(vertices.size() != 1) {
	    for(unsigned i = 0; i < vertices.size() - 1; ++i){
		edges.push_back(std::make_pair(i, i + 1));
		
	    }
	}

	return std::make_pair(vertices,edges);
      }

    };

  } /* pattern */

} /* graybat */
