/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @author Marius Melzer
 * @licence GPLv3
 *
 */

#pragma once
#include <vector>
#include <string>
#include <curand_kernel.h> /* curand_uniform */

#include <geometry.hpp>
#include <ConstHybridVector.hpp>

#define REFLECTION_SMALL 1E-3
#define SMALL 1E-5
#define VERY_SMALL 0.0


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
 * betaVolume beta values for all prisms ordered accordingly to the prismIDs:
 *            prismID = triangleID + layer * numberOfTriangles;
 *            therefore, all betaVolume for a layer are grouped together
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
 * triangleSurfaces the sizes of the surfaces of each triangle, ordered by the triangleID
 *
 * forbiddenEdge  describes the relation of edge indices of adjacent triangles
 *           -1 means, there is no adjacent triangle to that edge
 *           0,1,2 describes the index of the edge as seen from the ADJACENT triangle
 *
 *           order of data is similar to normalVec:
 *           [ triangle1edge0, triangle2edge0, ... triangleNedge0, triangle1edge1, triangle2edge1, ... ]
 *           i.e. all first edges of each triangle, followed by all second edges of each triangle and so on.
 *
 * trianglesPointIndices contains the indices to access the "points" datastructure 
 *           (each triangle has 3 points as vertices). Each entry is an
 *           index from 0 to numberOfPoints, corresponding to the positions 
 *           of a vertex in "points".
 *           structure is similar to "forbiddenEdge":
 *           [ triangle1A, triangle2A, ... triangleNA, triangle1B, triangle2B, ... triangleNB, triangle1C, ... ]
 *           i.e. for triangles with vertices A,B,C there are all the indices
 *           of the A-vertices, followed by all the B and C vertices.
 *
 * triangleNeighbors describes the relation of triangle indices to each other.
 *           Each entry corresponds to a triangleID (see "trianglePointIndices") which
 *           is adjacent to the current triangle.
 *           structure is similar to "forbiddenEdge":
 *           [ triangle1edge0, triangle2edge0, ... triangleNedge0, triangle1edge1, triangle2edge1, ... ]
 *
 * triangleNormalPoint contains indices to the x and y components of the positions where the
 *             normalVectors start (see normalVec). For each Triangle 3 points (3 sides)
 *             are stored in this list.
 *             Indices point to locations in "points" (i.e. normal vectors start at
 *             triangle vertices!)
 *             structure is VERY similar to trianglePointIndices: 
 *             [ triangle1p0, triangle2p0, ... triangleNp0, triangle1p1, triangle2p1, ... ]
 *
 * refractiveIndices [0]->bottomInside, [1]->bottomOutside, [2]->topInside, [3]->topOutside
 * 
 * reflectivities Contains the reflectivities for upper and lower surface of gain medium
 *                Structure is based on 2 layers of triangles:
 *                [refl_tri1_bott, refl_tri2_bott, ...,refl_triN_bott, refl_tri1_top, refl_tri2_top, ..., refl_triN_top]
 * 
 * totalReflectionAngles [0]-> bottomTotalReflectionAngle, [1]-> topTotalReflectionAngle
 */
class Mesh {
 public:

  Mesh(// Constants
       double claddingAbsorption,
       float surfaceTotal,
       float thickness,
       float nTot,
       float crystalTFluo,
       unsigned numberOfTriangles,
       unsigned numberOfLevels,
       unsigned numberOfPrisms,
       unsigned numberOfPoints,
       unsigned numberOfSamples,
       unsigned claddingNumber,
       // Vectors
       std::vector<double> points,
       std::vector<double> normalVec,
       std::vector<double> betaVolume,
       std::vector<double> centers,
       std::vector<float> triangleSurfaces,
       std::vector<int> forbiddenEdge,
       std::vector<double> betaCells,
       std::vector<unsigned> claddingCellTypes,
       std::vector<float> refractiveIndices,
       std::vector<float> reflectivities,
       std::vector<float> totalReflectionAngles,
       std::vector<unsigned> trianglePointIndices,
       std::vector<int> triangleNeighbors,
       std::vector<unsigned> triangleNormalPoint);


  ConstHybridVector<double>   points;
  ConstHybridVector<double>   betaVolume;
  ConstHybridVector<double>   normalVec;
  ConstHybridVector<double>   centers;
  ConstHybridVector<float>    triangleSurfaces;
  ConstHybridVector<int>      forbiddenEdge;
  ConstHybridVector<double>   betaCells;
  ConstHybridVector<unsigned> claddingCellTypes;

  ConstHybridVector<float>    refractiveIndices; 
  ConstHybridVector<float>    reflectivities;   //based on triangleIndex, with offset from bottom/top
  ConstHybridVector<float>    totalReflectionAngles;

  // Indexstructs
  ConstHybridVector<unsigned> trianglePointIndices;
  ConstHybridVector<int>      triangleNeighbors;
  ConstHybridVector<unsigned> triangleNormalPoint;

  // Constants
  double claddingAbsorption;
  float surfaceTotal;
  float thickness;
  float nTot;
  float crystalTFluo;
  unsigned numberOfTriangles;
  unsigned numberOfLevels;
  unsigned numberOfPrisms;
  unsigned numberOfPoints;
  unsigned numberOfSamples;
  unsigned claddingNumber;

  ~Mesh();

  void free();

  __device__ int getNeighbor(unsigned triangle, int edge) const;
  __device__ Point genRndPoint(unsigned triangle, unsigned level, curandStateMtgp32 *globalState) const;
  __device__ double getBetaVolume(unsigned triangle, unsigned level) const;
  __device__ double getBetaVolume(unsigned prism) const;
  __device__ NormalRay getNormal(unsigned triangle, int edge) const;
  __device__ Point getSamplePoint(unsigned sample) const;
  __device__ Point getCenterPoint(unsigned triangle, unsigned level) const;
  __device__ int getForbiddenEdge(unsigned triangle, int edge) const;
  __device__ unsigned getCellType(unsigned triangle) const;


  unsigned getMaxReflections(ReflectionPlane reflectionPlane) const;
  unsigned getMaxReflections() const;

  __device__ __host__ float getReflectivity(ReflectionPlane reflectionPlane, unsigned triangle) const;
  __device__ __host__ float getReflectionAngle(ReflectionPlane reflectionPlane) const;

  __device__ __host__ void test() const;

};

