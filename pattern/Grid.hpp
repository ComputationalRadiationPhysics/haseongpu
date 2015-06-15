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



    struct Grid {

      const unsigned height;
      const unsigned width;
	
      Grid(const unsigned height, const unsigned width) :
	height(height),
	width(width){

      }
    
      GraphDescription operator()(){

	const unsigned verticesCount = height * width;
	std::vector<unsigned> vertices(verticesCount);
	std::vector<EdgeDescription> edges;

	for(unsigned i = 0; i < vertices.size(); ++i){
	  vertices.at(i) = i;
	    
	  if(i >= width){
	    unsigned up   = i - width;
	    edges.push_back(std::make_pair(i, up));
	  }

	  if(i < (verticesCount - width)){
	    unsigned down = i + width;
	    edges.push_back(std::make_pair(i,down));
	  }


	  if((i % width) != (width - 1)){
	    int right = i + 1;
	    edges.push_back(std::make_pair(i,right));
	  }

	  if((i % width) != 0){
	    int left = i - 1;
	    edges.push_back(std::make_pair(i, left));
	  }
	

	}

	return std::make_pair(vertices,edges);
      }

    };

  } /* pattern */

} /* graybat */
