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


    struct FullyConnected {

      const unsigned verticesCount;

      FullyConnected(const unsigned verticesCount) :
	verticesCount(verticesCount){

      }
      
      GraphDescription operator()(){
	std::vector<VertexID> vertices(verticesCount);

	assert(vertices.size() == verticesCount);

	std::vector<EdgeDescription> edges;

	for(unsigned i = 0; i < vertices.size(); ++i){
	  vertices.at(i) = i;
	  for(unsigned j = 0; j < vertices.size(); ++j){
	    if(i == j){
	      continue;
	    } 
	    else {
	      edges.push_back(std::make_pair(i, j));
	    }
	  }

	}

	return std::make_pair(vertices,edges);
      }

    };
      
  } /* pattern */

} /* graybat */
