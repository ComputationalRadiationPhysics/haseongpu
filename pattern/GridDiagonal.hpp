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

    struct GridDiagonal {

	const unsigned height;
	const unsigned width;
	
	GridDiagonal(const unsigned height, const unsigned width) :
	    height(height),
	    width(width){

	}
    
	GraphDescription operator()(){
	    const unsigned verticesCount = height * width;
	    std::vector<unsigned> vertices(verticesCount);
	    std::vector<EdgeDescription > edges;

	    for(unsigned i = 0; i < vertices.size(); ++i){
		vertices.at(i) = i;
	    
		// UP
		if(i >= width){
		    unsigned up   = i - width;
		    edges.push_back(std::make_pair(i, up));
		}

		// UP LEFT
		if(i >= width and (i % width) != 0){
		    unsigned up_left   = i - width - 1;
		    edges.push_back(std::make_pair(i, up_left));
		}

		// UP RIGHT
		if(i >= width and (i % width) != (width - 1)){
		    unsigned up_right   = i - width + 1;
		    edges.push_back(std::make_pair(i, up_right));
		}

		// DOWN
		if(i < (verticesCount - width)){
		    unsigned down = i + width;
		    edges.push_back(std::make_pair(i, down));
		}

		// DOWN LEFT
		if(i < (verticesCount - width) and (i % width) != 0){
		    unsigned down_left = i + width - 1;
		    edges.push_back(std::make_pair(i, down_left));
		}

		// DOWN RIGHT
		if(i < (verticesCount - width) and (i % width) != (width - 1)){
		    unsigned down_right = i + width + 1;
		    edges.push_back(std::make_pair(i, down_right));
		}

		// RIGHT
		if((i % width) != (width - 1)){
		    int right = i + 1;
		    edges.push_back(std::make_pair(i, right));
		}

		// LEFT
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
