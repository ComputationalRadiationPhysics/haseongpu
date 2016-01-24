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

//CLIB
#include <cfloat> /* DBL_MAX*/

// STL
#include <vector>

// Alpaka
#include <alpaka/alpaka.hpp>

// HASEonGPU
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

template <typename T_Acc, typename T_Host, typename T_Stream>
class Mesh {
 public:

    using Host    = T_Host;
    using Acc     = T_Acc;     
    using Stream  = T_Stream;
    using DevAcc  = alpaka::dev::Dev<Acc>;
    using DevHost = alpaka::dev::Dev<Host>;        
    using Dim =  alpaka::dim::DimInt<1u>;
    using Size = std::size_t;


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

    // Buffers
    alpaka::mem::buf::Buf<DevAcc, double,   Dim, Size> pointsBuf;
    alpaka::mem::buf::Buf<DevAcc, double,   Dim, Size> normalVecBuf;    
    alpaka::mem::buf::Buf<DevAcc, double,   Dim, Size> betaVolumeBuf;
    alpaka::mem::buf::Buf<DevAcc, double,   Dim, Size> centersBuf;
    alpaka::mem::buf::Buf<DevAcc, float,    Dim, Size> triangleSurfacesBuf;
    alpaka::mem::buf::Buf<DevAcc, int,      Dim, Size> forbiddenEdgeBuf;
    alpaka::mem::buf::Buf<DevAcc, double,   Dim, Size> betaCellsBuf;    
    alpaka::mem::buf::Buf<DevAcc, unsigned, Dim, Size> claddingCellTypesBuf;            

    alpaka::mem::buf::Buf<DevAcc, float,    Dim, Size> refractiveIndicesBuf;
    alpaka::mem::buf::Buf<DevAcc, float,    Dim, Size> reflectivitiesBuf;  //based on triangleIndex, with offset from bottom/top
    alpaka::mem::buf::Buf<DevAcc, float,    Dim, Size> totalReflectionAnglesBuf;

    alpaka::mem::buf::Buf<DevAcc, unsigned, Dim, Size> trianglePointIndicesBuf;
    alpaka::mem::buf::Buf<DevAcc, int,      Dim, Size> triangleNeighborsBuf;
    alpaka::mem::buf::Buf<DevAcc, unsigned, Dim, Size> triangleNormalPointBuf;

    // Ptr
    double* points;
    double* normalVec;    
    double* betaVolume;
    double* centers;
    float*  triangleSurfaces;
    int*    forbiddenEdge;
    double* betaCells;    
    unsigned* claddingCellTypes;            

    float*    refractiveIndices;
    float*    reflectivities;  //based on triangleIndex* with offset from bottom/top
    float*    totalReflectionAngles;

    unsigned* trianglePointIndices;
    int*      triangleNeighbors;
    unsigned* triangleNormalPoint;
    
    
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
	 // Buffers
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
	 std::vector<unsigned> triangleNormalPoint,
	 // Device
	 DevAcc &dev) :
	// Constants
	claddingAbsorption(claddingAbsorption),
	surfaceTotal(surfaceTotal),
	thickness(thickness),
	nTot(nTot),
	crystalTFluo(crystalTFluo),
	numberOfTriangles(numberOfTriangles),
	numberOfLevels(numberOfLevels),
	numberOfPrisms(numberOfPrisms),
	numberOfPoints(numberOfPoints),
	numberOfSamples(numberOfSamples),
	claddingNumber(claddingNumber),
	// Vectors
	pointsBuf(alpaka::mem::buf::alloc                <double,   Size, Size, DevAcc>(dev, points.size())),
	normalVecBuf(alpaka::mem::buf::alloc             <double,   Size, Size, DevAcc>(dev, normalVec.size())),
	betaVolumeBuf(alpaka::mem::buf::alloc            <double,   Size, Size, DevAcc>(dev, betaVolume.size())),
	centersBuf(alpaka::mem::buf::alloc               <double,   Size, Size, DevAcc>(dev, centers.size())),
	triangleSurfacesBuf(alpaka::mem::buf::alloc      <float,    Size, Size, DevAcc>(dev, triangleSurfaces.size())),
	forbiddenEdgeBuf(alpaka::mem::buf::alloc         <int,      Size, Size, DevAcc>(dev, forbiddenEdge.size())),
	betaCellsBuf(alpaka::mem::buf::alloc             <double,   Size, Size, DevAcc>(dev, betaCells.size())),
	claddingCellTypesBuf(alpaka::mem::buf::alloc     <unsigned, Size, Size, DevAcc>(dev, claddingCellTypes.size())),
	refractiveIndicesBuf(alpaka::mem::buf::alloc     <float,    Size, Size, DevAcc>(dev, refractiveIndices.size())),
	reflectivitiesBuf(alpaka::mem::buf::alloc        <float,    Size, Size, DevAcc>(dev, reflectivities.size())),
	totalReflectionAnglesBuf(alpaka::mem::buf::alloc <float,    Size, Size, DevAcc>(dev, totalReflectionAngles.size())),
	trianglePointIndicesBuf(alpaka::mem::buf::alloc  <unsigned, Size, Size, DevAcc>(dev, trianglePointIndices.size())),
	triangleNeighborsBuf(alpaka::mem::buf::alloc     <int,      Size, Size, DevAcc>(dev, triangleNeighbors.size())),
	triangleNormalPointBuf(alpaka::mem::buf::alloc   <unsigned, Size, Size, DevAcc>(dev, triangleNormalPoint.size()))
    {

	// FIXIT: Is this the most general
        DevHost devHost (alpaka::dev::DevMan<Host>::getDevByIdx(0));
	Stream  stream  (dev);
	
	alpaka::mem::view::ViewPlainPtr<DevHost, double,   Dim, Size> hPoints(points.data(), devHost,  alpaka::Vec<Dim, Size>(points.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, double,   Dim, Size> hNormalVec(normalVec.data(), devHost, alpaka::Vec<Dim, Size>(normalVec.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, double,   Dim, Size> hBetaVolume(betaVolume.data(), devHost, alpaka::Vec<Dim, Size>(betaVolume.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, double,   Dim, Size> hCenters(centers.data(), devHost, alpaka::Vec<Dim, Size>(centers.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, float,    Dim, Size> hTriangleSurfaces(triangleSurfaces.data(), devHost, alpaka::Vec<Dim, Size>(triangleSurfaces.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, int,      Dim, Size> hForbiddenEdge(forbiddenEdge.data(), devHost, alpaka::Vec<Dim, Size>(forbiddenEdge.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, double,   Dim, Size> hBetaCells(betaCells.data(), devHost, alpaka::Vec<Dim, Size>(betaCells.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, unsigned, Dim, Size> hCladdingCellTypes(claddingCellTypes.data(), devHost, alpaka::Vec<Dim, Size>(claddingCellTypes.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, float,    Dim, Size> hRefractiveIndices(refractiveIndices.data(), devHost, alpaka::Vec<Dim, Size>(refractiveIndices.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, float,    Dim, Size> hReflectivities(reflectivities.data(), devHost, alpaka::Vec<Dim, Size>(reflectivities.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, float,    Dim, Size> hTotalReflectionsAngles(totalReflectionAngles.data(), devHost, alpaka::Vec<Dim, Size>(totalReflectionAngles.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, unsigned, Dim, Size> hTrianglePointIndices(trianglePointIndices.data(), devHost, alpaka::Vec<Dim, Size>(trianglePointIndices.size()));
	alpaka::mem::view::ViewPlainPtr<DevHost, int     , Dim, Size> hTriangleNeighbors(triangleNeighbors.data(), devHost, alpaka::Vec<Dim, Size>(triangleNeighbors.size()));
        alpaka::mem::view::ViewPlainPtr<DevHost, unsigned, Dim, Size> hTriangleNormalPoint(triangleNormalPoint.data(), devHost, alpaka::Vec<Dim, Size>(triangleNormalPoint.size()));		

	alpaka::mem::view::copy(stream, this->pointsBuf, hPoints, points.size());
	alpaka::mem::view::copy(stream, this->normalVecBuf, hNormalVec, normalVec.size());
	alpaka::mem::view::copy(stream, this->betaVolumeBuf, hBetaVolume, betaVolume.size());
	alpaka::mem::view::copy(stream, this->centersBuf, hCenters, centers.size());
	alpaka::mem::view::copy(stream, this->triangleSurfacesBuf, hTriangleSurfaces, triangleSurfaces.size());
	alpaka::mem::view::copy(stream, this->forbiddenEdgeBuf, hForbiddenEdge, forbiddenEdge.size());
	alpaka::mem::view::copy(stream, this->betaCellsBuf, hBetaCells, betaCells.size());
	alpaka::mem::view::copy(stream, this->claddingCellTypesBuf, hCladdingCellTypes, claddingCellTypes.size());
	alpaka::mem::view::copy(stream, this->refractiveIndicesBuf, hRefractiveIndices, refractiveIndices.size());
	alpaka::mem::view::copy(stream, this->reflectivitiesBuf, hReflectivities, reflectivities.size());

        alpaka::mem::view::copy(stream, this->totalReflectionAnglesBuf, hTotalReflectionsAngles, totalReflectionAngles.size());
	alpaka::mem::view::copy(stream, this->trianglePointIndicesBuf, hTrianglePointIndices, trianglePointIndices.size());
	alpaka::mem::view::copy(stream, this->triangleNeighborsBuf, hTriangleNeighbors, triangleNeighbors.size());
	alpaka::mem::view::copy(stream, this->triangleNormalPointBuf, hTriangleNormalPoint, triangleNormalPoint.size());

        this->points = alpaka::mem::view::getPtrNative(pointsBuf);
        this->normalVec = alpaka::mem::view::getPtrNative(normalVecBuf);    
        this->betaVolume = alpaka::mem::view::getPtrNative(betaVolumeBuf);
        this->centers = alpaka::mem::view::getPtrNative(centersBuf);
        this->triangleSurfaces = alpaka::mem::view::getPtrNative(triangleSurfacesBuf);
        this->forbiddenEdge = alpaka::mem::view::getPtrNative(forbiddenEdgeBuf);
        this->betaCells = alpaka::mem::view::getPtrNative(betaCellsBuf);    
        this->claddingCellTypes = alpaka::mem::view::getPtrNative(claddingCellTypesBuf);            
        this->refractiveIndices = alpaka::mem::view::getPtrNative(refractiveIndicesBuf);
        this->reflectivities = alpaka::mem::view::getPtrNative(reflectivitiesBuf);
        this->totalReflectionAngles = alpaka::mem::view::getPtrNative(totalReflectionAnglesBuf);
        this->trianglePointIndices = alpaka::mem::view::getPtrNative(trianglePointIndicesBuf);
        this->triangleNeighbors = alpaka::mem::view::getPtrNative(triangleNeighborsBuf);
        this->triangleNormalPoint = alpaka::mem::view::getPtrNative(triangleNormalPointBuf);
        
    }

    template <class T, class B, class E>
    void assertRange(const std::vector<T> &v, const B minElement,const E maxElement, const bool equals){
	if(equals){
	    assert(*std::min_element(v.begin(),v.end()) == minElement);
	    assert(*std::max_element(v.begin(),v.end()) == maxElement);
	}else{
	    assert(*std::min_element(v.begin(),v.end()) >= minElement);
	    assert(*std::max_element(v.begin(),v.end()) <= maxElement);
	}
    }

    template <class T, class B>
    void assertMin(const std::vector<T>& v,const  B minElement,const bool equals){
	if(equals){
	    assert(*std::min_element(v.begin(),v.end()) == minElement);
	}else{
	    assert(*std::min_element(v.begin(),v.end()) >= minElement);
	}
    }


    /**
     * @brief fetch the id of an adjacent triangle
     *
     * @param triangle the index of the triangle of which you want the neighbor
     * @param edge the side of the triangle, for whih you want the neighbor
     *
     * @return the index of the neighbor triangle
     */
    ALPAKA_FN_ACC int getNeighbor(unsigned triangle, int edge) const{
	return triangleNeighbors[triangle + edge*numberOfTriangles];
    }

    /**
     * @brief generate a random point within a prism
     *
     * @param triangle the triangle to describe the desired prism
     * @param the level of the desired prism
     * @param *globalState a global state for the Mersenne Twister PRNG
     *
     * @return random 3D point inside the desired prism
     *
     * Uses a Mersenne Twister PRNG and Barycentric coordinates to generate a
     * random position inside a given triangle in a specific depth
     */
    


    //FIXIT: use random number generator of alpaka (picongpu: src/libPMACC/startposition/RandImpl)

    template <typename T_Rand>
    ALPAKA_FN_ACC Point genRndPoint(T_Rand &rand, unsigned triangle, unsigned level) const{
	// Random number generator
	// FIXIT: No need to initialize this again and again ?
	
	Point startPoint = {0,0,0};
	double u = rand();
	double v = rand();

	if((u+v)>1) {
		u = 1-u;
		v = 1-v;
	}
	
	double w = 1-u-v;
	int t1 = trianglePointIndices[triangle];
	int t2 = trianglePointIndices[triangle + numberOfTriangles];
	int t3 = trianglePointIndices[triangle + 2 * numberOfTriangles];

	// convert the random startpoint into coordinates
	startPoint.z = (level + rand()) * thickness;
	startPoint.x =
	    (points[t1] * u) +
	    (points[t2] * v) +
	    (points[t3] * w);
	startPoint.y =
	    (points[t1+numberOfPoints] * u) +
	    (points[t2+numberOfPoints] * v) +
	    (points[t3+numberOfPoints] * w);

	return startPoint;
    }



    /**
     * @brief get a betaVolume for a specific triangle and level
     *
     * @param triangle the id of the desired triangle
     * @param level the level of the desired prism
     *
     * @return a beta value
     */
    ALPAKA_FN_ACC double getBetaVolume(unsigned triangle, unsigned level) const{
	return betaVolume[triangle + level * numberOfTriangles];
    }

    /**
     * @brief get a betaVolume for a specific prism
     *
     * @param the id of the desired prism
     *
     * @return a beta value
     */
    ALPAKA_FN_ACC double getBetaVolume(unsigned prism) const{
	return betaVolume[prism];
    }

    /**
     * @brief generates a normal vector for a given side of a triangle
     *
     * @param triangle the desired triangle
     * @param the edge (0,1,2) of the triangle
     *
     * @return a normal vector with length 1
     */
    ALPAKA_FN_ACC NormalRay getNormal(unsigned triangle, int edge) const{
	NormalRay ray = { {0,0},{0,0}};
	int offset =  edge*numberOfTriangles + triangle;
	ray.p.x = points[triangleNormalPoint[offset] ];
	ray.p.y = points[triangleNormalPoint[offset] + numberOfPoints ];

	ray.dir.x = normalVec[offset];
	ray.dir.y = normalVec[offset + 3*numberOfTriangles];

	return ray;
    }	

    /**
     * @brief genenerates a point with the coordinates of a given vertex
     *
     * @param sample the id of the desired samplepoint
     *
     * @return the Point with correct 3D coordinates
     */
    ALPAKA_FN_ACC Point getSamplePoint(unsigned sample_i) const{
	Point p = {0,0,0};
	unsigned level = sample_i/numberOfPoints;
	p.z = level*thickness;
	unsigned pos = sample_i - (numberOfPoints * level);
	p.x = points[pos];
	p.y = points[pos + numberOfPoints];
	return p;
    }

    /**
     * @brief get a Point in the center of a prism
     *
     * @param triangle the id of the desired triangle
     * @param level the level of the desired prism
     *
     * @return a point with the coordinates (3D) of the prism center
     */
    ALPAKA_FN_ACC Point getCenterPoint(unsigned triangle,unsigned level) const{
	Point p = {0,0,(level+0.5)*thickness};
	p.x = centers[triangle];
	p.y = centers[triangle + numberOfTriangles];
	return p;
    }

    /**
     * @brief gets the edge-id which will be forbidden
     *
     * @param trianle the index of the triangle you are currently in
     * @param edge the index of the edge through which you are leaving the triangle
     *
     * @return the id of the edge, which will be forbidden in the new triangle
     *
     * The forbidden edge corresponds to the edge through which you left the
     * previous triangle (has a different index in the new triangle)
     */
    ALPAKA_FN_ACC int getForbiddenEdge(unsigned triangle,int edge) const{
	return forbiddenEdge[edge * numberOfTriangles + triangle];
    }


    ALPAKA_FN_ACC unsigned getCellType(unsigned triangle) const{
	return claddingCellTypes[triangle];
    }


    ALPAKA_FN_HOST double distance2D(const TwoDimPoint p1, const TwoDimPoint p2) const{
	return fabs(sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)));
    }

    ALPAKA_FN_HOST double getMaxDistance(std::vector<TwoDimPoint> points) const{
	double maxDistance = -1;

	for(unsigned p1=0 ; p1 < points.size() ; ++p1)
	    for(unsigned p2 = p1; p2 < points.size() ; ++p2)
		maxDistance = std::max(maxDistance,distance2D(points[p1],points[p2]));

	return maxDistance;
    }

    ALPAKA_FN_HOST double calculateMaxDiameter(const double* points, const unsigned offset) const{
	TwoDimPoint minX = {DBL_MAX,0};
	TwoDimPoint minY = {0,DBL_MAX};
	TwoDimPoint maxX = {DBL_MIN,0};
	TwoDimPoint maxY = {0,DBL_MIN};

	for(unsigned p=0; p<offset; ++p){
	    TwoDimPoint np = {points[p], points[p+offset]};
	    minX = (points[p] < minX.x) ? np : minX;
	    maxX = (points[p] > maxX.x) ? np : maxX;
	}
	for(unsigned p=offset;p<2*offset;++p){
	    TwoDimPoint np = {points[p-offset],points[p]};
	    minY = points[p] < minY.y ? np : minY;
	    maxY = points[p] > maxY.y ? np : maxY;
	}

	std::vector<TwoDimPoint> extrema;
	extrema.push_back(minX);
	extrema.push_back(minY);
	extrema.push_back(maxX);
	extrema.push_back(maxY);
	

	return getMaxDistance(extrema);
    }

    ALPAKA_FN_HOST unsigned getMaxReflections (ReflectionPlane reflectionPlane) const{
	double d    = calculateMaxDiameter(points, numberOfPoints);
	float alpha = getReflectionAngle(reflectionPlane) * M_PI / 180.;
	double h    = numberOfLevels * thickness; 
	double z    = d/tan(alpha);
	return ceil(z/h);
    }

    ALPAKA_FN_HOST unsigned getMaxReflections() const{
	unsigned top = getMaxReflections(TOP_REFLECTION);
	unsigned bottom = getMaxReflections(BOTTOM_REFLECTION);
	return std::max(top, bottom);
    }

    ALPAKA_FN_HOST_ACC float getReflectivity(ReflectionPlane reflectionPlane, unsigned triangle) const{
	switch(reflectionPlane){
	case BOTTOM_REFLECTION:
	    return reflectivities[triangle];
	case TOP_REFLECTION:
	    return reflectivities[triangle + numberOfTriangles];
	}
	return 0;
    }

    ALPAKA_FN_HOST_ACC float getReflectionAngle(ReflectionPlane reflectionPlane) const{
	switch(reflectionPlane){
	case BOTTOM_REFLECTION:
	    return totalReflectionAngles[0];
	    //return asin(refractiveIndices[1]/refractiveIndices[0]);
	case TOP_REFLECTION:
	    return totalReflectionAngles[1];
	    //return asin(refractiveIndices[3]/refractiveIndices[2]);
	}
	return  0;
    }

    ALPAKA_FN_HOST_ACC float getTriangleSurface(const unsigned triangle) const{
	return triangleSurfaces[triangle];
	    
    }

    ALPAKA_FN_HOST_ACC  void test() const{
	printf("Constants:\n");
	printf("claddingAbsorption: %f\n", claddingAbsorption);
	printf("surfaceTotal: %f\n", surfaceTotal);
	printf("thickness: %f\n", thickness);
	printf("nTot: %f\n", nTot);
	printf("crystalTFluo: %f\n", crystalTFluo);
	printf("numberOfTriangles: %u\n", numberOfTriangles);
	printf("numberOfLevels: %u\n", numberOfLevels);
	printf("numberOfPrisms: %u\n", numberOfPrisms);
	printf("numberOfPoints: %u\n", numberOfPoints);
	printf("numberOfSamples: %u\n", numberOfSamples);
	printf("claddingNumber: %u\n", claddingNumber);
  
    }

};

