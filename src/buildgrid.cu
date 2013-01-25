#include <cstdio>
#include <cassert>
#include <cmath>
#include <list>
#include <vector>
#include <algorithm>

#include "buildgrid.h"

//set numbers of cells in each direction of the grid
//Input:
//numbers of triangles in a layer and numbers of layers
//numbers of prisms per cell if they are in a homogeneous range
void setDimGrid(Grid *grid, int numTrianglesPerLayer, int numLayers, int numPrismPerCell) {
    float invQuotient;
    invQuotient = pow(numPrismPerCell, -0.3333);

    grid->dim.z = numLayers * invQuotient;
    if ( grid->dim.z < 1 ) grid->dim.z = 1;

    grid->dim.x = grid->dim.y = sqrtf(numTrianglesPerLayer * invQuotient);
    if ( grid->dim.x < 1 ) grid->dim.x = 1;
    if ( grid->dim.y < 1 ) grid->dim.y = 1;

}

//determines the boudning box of the grid with the boudning box array of the prisms
void calcGridAabb(Grid *grid, const Aabb *objectAabbs, int numAabbs){
    
    int i;
    //BoundingBox of the Grid
    //set the start aabb of the grid with the aabb-parameters of the first aabb element
    grid->aabb.max.x = objectAabbs->max.x;
    grid->aabb.min.x = objectAabbs->min.x;
    grid->aabb.max.y = objectAabbs->max.y;
    grid->aabb.min.y = objectAabbs->min.y;
    grid->aabb.max.z = objectAabbs->max.z;
    grid->aabb.min.z = objectAabbs->min.z;

    for (i = 1; i < numAabbs; ++i) {
        if (grid->aabb.max.x < objectAabbs[i].max.x)
            grid->aabb.max.x = objectAabbs[i].max.x;
        if (grid->aabb.min.x > objectAabbs[i].min.x)
            grid->aabb.min.x = objectAabbs[i].min.x;
        if (grid->aabb.max.y < objectAabbs[i].max.y)
            grid->aabb.max.y = objectAabbs[i].max.y;
        if (grid->aabb.min.y > objectAabbs[i].min.y)
            grid->aabb.min.y = objectAabbs[i].min.y;
        if (grid->aabb.max.z < objectAabbs[i].max.z)
            grid->aabb.max.z = objectAabbs[i].max.z;
        if (grid->aabb.min.z > objectAabbs[i].min.z)
            grid->aabb.min.z = objectAabbs[i].min.z;
    }
}

//buildGrid
//fillGrid
//build and fill grid
//before using this function set the dimension of the grid !!
//   therefore you can use setDimGrid
void buildGrid(Grid *grid, const Aabb *aabbs, int numAabbs) {

    int numCells = grid->dim.x * grid->dim.y * grid->dim.z;
    grid->cellList = new GridCell[numCells];

    std::vector<std::list<int> > cellLists(numCells);

    int i, j, k, l;
    //BoundingBox of the Grid
    calcGridAabb(grid, aabbs, numAabbs);

    float deltaX,deltaY,deltaZ;
    deltaX=grid->aabb.max.x - grid->aabb.min.x;
    deltaY=grid->aabb.max.y - grid->aabb.min.y;
    deltaZ=grid->aabb.max.z - grid->aabb.min.z;

//    //increase Boundingbox 0.1%
//    deltaX*=1.001;
//    deltaY*=1.001;
//    deltaZ*=1.001;
//    grid->aabb.max.x+=deltaX*0.001;
//    grid->aabb.max.y+=deltaY*0.001;
//    grid->aabb.max.z+=deltaZ*0.001;
    
    int imin, imax, jmin, jmax, kmin, kmax;

    float tmpX = grid->dim.x / deltaX;
    float tmpY = grid->dim.y / deltaY;
    float tmpZ = grid->dim.z / deltaZ;

    int cellIndex;
    int dimXdimY=grid->dim.x*grid->dim.y;
    
    //*
    //fill Grid
    //loop over prism bounding boxes
    for (l = 0; l < numAabbs; ++l) {
        //determine cell indexes min and max for each dimension
        imin = (aabbs[l].min.x - grid->aabb.min.x) * tmpX;
        imax = (aabbs[l].max.x - grid->aabb.min.x) * tmpX;
        if(imax==grid->dim.x) --imax;
        
        jmin = (aabbs[l].min.y - grid->aabb.min.y) * tmpY;
        jmax = (aabbs[l].max.y - grid->aabb.min.y) * tmpY;
        if(jmax==grid->dim.y) --jmax;
        
        kmin = (aabbs[l].min.z - grid->aabb.min.z) * tmpZ;
        kmax = (aabbs[l].max.z - grid->aabb.min.z) * tmpZ;
        if(kmax==grid->dim.z) --kmax;
        
        //add prism index to cellList if prism is inside cell
        for (i = imin; i <= imax; ++i)
            for (j = jmin; j <= jmax; ++j)
                for (k = kmin; k <= kmax; ++k){
                    cellIndex = k * dimXdimY + j * grid->dim.x + i;
                    cellLists[cellIndex].push_back(l);
                }
    }
    
    //*
    //copy lists to grid structure
    
    //pointer
    int* p;
    //std::list iterator
    std::list<int>::const_iterator it;
    
    //loop over grid cells
    for(l=0;l<numCells;++l){
        grid->cellList[l].length = cellLists[l].size();
        grid->cellList[l].prismIdxList=new int[grid->cellList->length];
        //set pointer to array begin
        p=grid->cellList[l].prismIdxList;
        //copy
        for(it=cellLists[l].begin();it!=cellLists[l].end();++it)
            *p++=*it;
    }
}
