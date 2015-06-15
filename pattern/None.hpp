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

    
    struct None {
      GraphDescription operator()(){
	  std::vector<VertexID> vertices;
    
	  std::vector<EdgeDescription> edges;

	  return std::make_pair(vertices,edges);
      }

    };

  } /* pattern */

} /* graybat */
