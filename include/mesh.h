#ifndef MESH_H
#define MESH_H 

#include <vector>
#include <curand_kernel.h> /* curand_uniform */
#include <string>

#define SMALL 1E-06
#define VERY_SMALL 0.0

struct TwoDimPoint {
  double x;
  double y;
};

typedef TwoDimPoint TwoDimDir;

struct Point {
  double x;
  double y;
  double z;

};

typedef Point Vector;

/**
 * @brief a Ray, defined by a startpoint, direction and length
 */
struct Ray {
  Point p;
  Vector dir;
  float length;
  
};

struct NormalRay {
  TwoDimPoint p;
  TwoDimDir dir;
};



/**
 * @brief Contains the structure of the crystal
 *
 * All the fixed values of how the crystal is meshed
 *
 * points The coordinates of the triangle vertices
 *        All x coordinates followed by all of the y coordinates of the triangle vertices
 *        structure: [x_1, x_2, ... x_n, y_1, y_2, ... y_n] (n == numberOfPoints)
 *
 *
 * betaValues beta values for all prisms ordered accordingly to the prismIDs:
 *            prismID = triangleID + layer * numberOfTriangles;
 *            therefore, all betaValues for a layer are grouped together
 *
 * normalVec the normal vectors for each triangle edge
 *           first half (size: 3*numberOfTriangles -> one for each side) contains
 *           the x components of each vector, second half contains the y components.
 *           the each half is ordered as follows:
 *           [ triangle1edge0, triangle2edge0, ... triangleNedge0, triangle1edge1, triangle2edge1, ... ]
 *           i.e. all first edges of each triangle, followed by all second edges of each triangle and so on.
 *
 * centers the coordinates of the center points for each triangle
 *         All x coordinates followed by all y coordinates of the triangle vertices
 *         similar to "points"
 *
 * surfaces the sizes of the surfaces of each triangle, ordered by the triangleID
 *
 * forbidden describes the relation of edge indices of adjacent triangles
 *           -1 means, there is no adjacent triangle to that edge
 *           0,1,2 describes the index of the edge as seen from the ADJACENT triangle
 *
 *           order of data is similar to normalVec:
 *           [ triangle1edge0, triangle2edge0, ... triangleNedge0, triangle1edge1, triangle2edge1, ... ]
 *           i.e. all first edges of each triangle, followed by all second edges of each triangle and so on.
 *
 * triangles contains the indices to access the "points" datastructure 
 *           (each triangle has 3 points as vertices). Each entry is an
 *           index from 0 to numberOfPoints, corresponding to the positions 
 *           of a vertex in "points".
 *           structure is similar to "forbidden":
 *           [ triangle1A, triangle2A, ... triangleNA, triangle1B, triangle2B, ... triangleNB, triangle1C, ... ]
 *           i.e. for triangles with vertices A,B,C there are all the indices
 *           of the A-vertices, followed by all the B and C vertices.
 *
 * neighbors describes the relation of triangle indices to each other.
 *           Each entry corresponds to a triangleID (see "triangles") which
 *           is adjacent to the current triangle.
 *           structure is similar to "forbidden":
 *           [ triangle1edge0, triangle2edge0, ... triangleNedge0, triangle1edge1, triangle2edge1, ... ]
 *
 * normalPoint contains indices to the x and y components of the positions where the
 *             normalVectors start (see normalVec).
 *             Indices point to locations in "points" (i.e. normal vectors start at
 *             triangle vertices!)
 *             structure is VERY similar to triangles: 
 *             [ triangle1p0, triangle2p0, ... triangleNp0, triangle1p1, triangle2p1, ... ]
 *
 */
struct Mesh {

  double *points;
  double *betaValues;
  double *normalVec;
  double *centers;
  float  *surfaces;
  int	 *forbidden;
  float  *betaCells;

  //indexstructs
  unsigned *triangles;
  int *neighbors;
  unsigned *normalPoint;

  //constants
  float surfaceTotal;
  float thickness;
  float nTot;
  float crystalFluorescence;
  unsigned numberOfTriangles;
  unsigned numberOfLevels;
  unsigned numberOfPrisms;
  unsigned numberOfPoints;
  unsigned numberOfSamples;

  ~Mesh();

  __device__ int getNeighbor(unsigned triangle, int edge);
  __device__ Point genRndPoint(unsigned triangle, unsigned level, curandStateMtgp32 *globalState);
  __device__ double getBetaValue(unsigned triangle, unsigned level);
  __device__ double getBetaValue(unsigned prism);
  __device__ NormalRay getNormal(unsigned triangle, int edge);
  __device__ Point getSamplePoint(unsigned sample);
  __device__ Point getCenterPoint(unsigned triangle, unsigned level);
  __device__ int getForbiddenEdge(unsigned triangle, int edge);

  static int parseMultiGPU(Mesh *hMesh, 
			    Mesh **dMesh, 
			    std::string root,
			    unsigned numberOfDevices,
			    unsigned *devices);
};

#endif /* MESH_H */
