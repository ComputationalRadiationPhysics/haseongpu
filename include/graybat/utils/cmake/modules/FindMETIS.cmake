  #
  # Find the METIS includes and libraries
  #
  # ParMETIS is an MPI-based parallel library that implements a variety of algorithms for 
  # partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of 
  # sparse matrices. It can be found at:
  # 	http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
  #
  # METIS_INCLUDE_DIR - where to find autopack.h
  # METIS_LIBRARIES   - List of fully qualified libraries to link against.
  # METIS_FOUND       - Do not attempt to use if "no" or undefined.
  
  FIND_PATH(METIS_INCLUDE_DIR metis.h
    /usr/local/include
    /usr/include
  )
  
  FIND_LIBRARY(METIS_LIBRARY metis
    /usr/local/lib
    /usr/lib
  )
  
  FIND_LIBRARY(METIS_LIBRARY metis
    /usr/local/lib
    /usr/lib
  )
  
  IF(METIS_INCLUDE_DIR)
    IF(METIS_LIBRARY)
      SET( METIS_LIBRARIES ${METIS_LIBRARY} ${METIS_LIBRARY})
      SET( METIS_FOUND "YES" )
    ENDIF(METIS_LIBRARY)
  ENDIF(METIS_INCLUDE_DIR)