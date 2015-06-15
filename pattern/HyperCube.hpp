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

    
    struct HyperCube {

      const unsigned dimension;

      HyperCube(const unsigned dimension) :
	dimension(dimension){

      }
      
      unsigned hammingDistance(unsigned a, unsigned b){
	unsigned abXor = a xor b;
	return (unsigned) __builtin_popcount(abXor);
      }

      GraphDescription operator()(){
	assert(dimension >= 1);
	std::vector<EdgeDescription> edges;

	unsigned verticesCount = pow(2, dimension);
	std::vector<VertexID> vertices(verticesCount);

	for(unsigned i = 0; i < vertices.size(); ++i){
	  vertices.at(i) = i;
	  for(unsigned j = 0; j < vertices.size(); ++j){
	    if(hammingDistance(i, j) == 1){
	      edges.push_back(std::make_pair(i, j));
	    }

	  }
	}
    
	return std::make_pair(vertices,edges);
      }

    };

  } /* pattern */

} /* graybat */
