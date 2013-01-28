/* 
 * File:   main.cpp
 * Author: s7740034
 *
 * Created on 26. Januar 2013, 15:18
 */

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;

#include "grid.h"
#include "buildgrid.h"
#include "filtergrid.h"

#define CUDA_CHECK(cmd) {cudaError_t error = cmd; if(error!=cudaSuccess){printf("<%s>:%i ",__FILE__,__LINE__); printf("[CUDA] Error: %s\n", cudaGetErrorString(error));}}
/*start kernel, wait for finish and check errors*/
#define CUDA_CHECK_KERNEL_SYNC(...) __VA_ARGS__;CUDA_CHECK(cudaDeviceSynchronize())
/*only check if kernel start is valid*/
#define CUDA_CHECK_KERNEL(...) __VA_ARGS__;CUDA_CHECK(cudaGetLastError())

void cellIdx2ijk(const Grid *grid, int cellIdx, int *i, int *j, int *k){
	
	int dimXdimY=grid->dim.x*grid->dim.y;
	int tmp;
	
	*k = cellIdx / dimXdimY;
	
	tmp=cellIdx - *k * dimXdimY;
	//*j = (cellIdx - *k * dimXdimY) / grid->dimX;
	*j = tmp / grid->dim.x;
	
	//*i = cellIdx - *k * dimXdimY - *j* grid->dimX;
	*i = tmp - *j * grid->dim.x;
}

__global__ void testkernel(Grid grid, int *results, int *resultsize) {
  RayCu ray;
  PointCu start;
  start.x = 1.1;
  start.y = 1.1;
  start.z = 1.1;

  VectorCu dir;
  dir.x = 3.5;
  dir.y = 3.5;
  dir.z = 24;

  ray.P = start;
  ray.direction = dir;

  results = filter(&grid, &ray, resultsize);

  for(int i=0; i<*resultsize; ++i) {
    printf("%i, ", results[i]);
  }
}

/*
 * 
 */
int main(int argc, char** argv) 
{

    //height of a layer
    float hLayer=3.0f;
    
    //number of points per layer
    int numPperLayer=18;
    
    //number of triangles per layer
    int numTperLayer=21;
    
    int numLayers=9;
    
    //number of points
    int numP=numPperLayer*numLayers;
    
    //number of triangles
    int numT=numTperLayer*numLayers;
    
    PointCu pList[numP];
    
    pList[0].x=1.0f;  pList[0].y=1.0f;  pList[0].z=0.0f;  pList[0].w=hLayer;
    pList[1].x=1.0f;  pList[1].y=2.0f;  pList[1].z=0.0f;  pList[1].w=hLayer;
    pList[2].x=1.0f; pList[2].y=3.0f;  pList[2].z=0.0f; pList[2].w=hLayer;
    pList[3].x=1.0f; pList[3].y=4.0f;  pList[3].z=0.0f; pList[3].w=hLayer;
    
    pList[4].x=2.0f; pList[4].y=0.5f;  pList[4].z=0.0f; pList[4].w=hLayer;
    pList[5].x=2.0f; pList[5].y=1.5f;  pList[5].z=0.0f; pList[5].w=hLayer;
    pList[6].x=2.0f; pList[6].y=2.5f;  pList[6].z=0.0f; pList[6].w=hLayer;
    pList[7].x=2.0f; pList[7].y=3.5f;  pList[7].z=0.0f; pList[7].w=hLayer;
    pList[8].x=2.0f; pList[8].y=4.5f;  pList[8].z=0.0f; pList[8].w=hLayer;
    
    pList[9] .x=3.0f; pList[9].y=1.0f;  pList[9].z=0.0f; pList[9].w=hLayer;
    pList[10].x=3.0f; pList[10].y=2.0f; pList[10].z=0.0f; pList[10].w=hLayer;
    pList[11].x=3.0f; pList[11].y=3.0f; pList[11].z=0.0f; pList[11].w=hLayer;
    pList[12].x=3.0f; pList[12].y=4.0f; pList[12].z=0.0f; pList[12].w=hLayer;
    
    pList[13].x=4.0f; pList[13].y=0.5f; pList[13].z=0.0f; pList[13].w=hLayer;
    pList[14].x=4.0f; pList[14].y=1.5f; pList[14].z=0.0f; pList[14].w=hLayer;
    pList[15].x=4.0f; pList[15].y=2.5f; pList[15].z=0.0f; pList[15].w=hLayer;
    pList[16].x=4.0f; pList[16].y=3.5f; pList[16].z=0.0f; pList[16].w=hLayer;
    pList[17].x=4.0f; pList[17].y=4.5f; pList[17].z=0.0f; pList[17].w=hLayer;    
    

    //Test Zuweisungsoperator
    //PointCu ff;
    //ff=pList[0];
    //cout<<"ff="<<ff<<endl;
    
    //#################################
    //prisms in the first layer:
    PrismCu prismList[numT];
    
    prismList[0].t1.A=pList[0];  prismList[0].t1.B=pList[5];  prismList[0].t1.C=pList[4];
    
    prismList[1].t1.A=pList[0];  prismList[1].t1.B=pList[1];  prismList[1].t1.C=pList[5];
    
    prismList[2].t1.A=pList[1];  prismList[2].t1.B=pList[6];  prismList[2].t1.C=pList[5];
    
    prismList[3].t1.A=pList[1];  prismList[3].t1.B=pList[2];  prismList[3].t1.C=pList[6];
    
    prismList[4].t1.A=pList[2];  prismList[4].t1.B=pList[7];  prismList[4].t1.C=pList[6];
    
    prismList[5].t1.A=pList[2];  prismList[5].t1.B=pList[3];  prismList[5].t1.C=pList[7];
    
    prismList[6].t1.A=pList[3];  prismList[6].t1.B=pList[8];  prismList[6].t1.C=pList[7];

    
    prismList[7].t1.A=pList[4];  prismList[7].t1.B=pList[5];  prismList[7].t1.C=pList[9];
    
    prismList[8].t1.A=pList[5];  prismList[8].t1.B=pList[10];  prismList[8].t1.C=pList[9];
   
    prismList[9].t1.A=pList[5];  prismList[9].t1.B=pList[6];  prismList[9].t1.C=pList[10];
    
    prismList[10].t1.A=pList[6];  prismList[10].t1.B=pList[11];  prismList[10].t1.C=pList[10];
    
    prismList[11].t1.A=pList[6];  prismList[11].t1.B=pList[7];  prismList[11].t1.C=pList[11];
    
    prismList[12].t1.A=pList[7];  prismList[12].t1.B=pList[12];  prismList[12].t1.C=pList[11];
    
    prismList[13].t1.A=pList[7];  prismList[13].t1.B=pList[8];  prismList[13].t1.C=pList[12];
    
    
    prismList[14].t1.A=pList[9];  prismList[14].t1.B=pList[14];  prismList[14].t1.C=pList[13];
    
    prismList[15].t1.A=pList[9];  prismList[15].t1.B=pList[10];  prismList[15].t1.C=pList[14];
    
    prismList[16].t1.A=pList[10];  prismList[16].t1.B=pList[15];  prismList[16].t1.C=pList[14];
    
    prismList[17].t1.A=pList[10];  prismList[17].t1.B=pList[11];  prismList[17].t1.C=pList[15];
    
    prismList[18].t1.A=pList[11];  prismList[18].t1.B=pList[16];  prismList[18].t1.C=pList[15];
    
    prismList[19].t1.A=pList[11];  prismList[19].t1.B=pList[12];  prismList[19].t1.C=pList[16];
    
    prismList[20].t1.A=pList[12];  prismList[20].t1.B=pList[17];  prismList[20].t1.C=pList[16];
    
    
    //###################################
    //points in other Layers
    int i,j;
    int offsetIdxP;
    
    int offsetIdxPrism;
    
    float zValue;
    for(i=0;i<numLayers;++i){
        offsetIdxP=i*numPperLayer;
        zValue=i*hLayer;
        
        //sinnvoll??
        for(j=0;j<numPperLayer;++j){
            pList[offsetIdxP+j].x=pList[j].x;
            pList[offsetIdxP+j].y=pList[j].y;
            pList[offsetIdxP+j].z=zValue;
            pList[offsetIdxP+j].w=hLayer;
        }
        
        //Prisms
        offsetIdxPrism=i*numTperLayer;
        for(j=0;j<numTperLayer;++j){
            prismList[offsetIdxPrism+j].t1.A=prismList[j].t1.A;
            prismList[offsetIdxPrism+j].t1.A.z+=zValue;
            prismList[offsetIdxPrism+j].t1.B=prismList[j].t1.B;
            prismList[offsetIdxPrism+j].t1.B.z+=zValue;
            prismList[offsetIdxPrism+j].t1.C=prismList[j].t1.C;
            prismList[offsetIdxPrism+j].t1.C.z+=zValue;
        }
    }
    
    
    //calc AABBs
    Aabb aabbs[numT];
    calcAabbs(prismList,aabbs,numT);
    
    
    //for(i=0;i<numT;++i){
    //    cout<<"i="<<i<<" \tmin="<<aabbs[i].min<<" \tmax="<<aabbs[i].max<<endl;
    //}
    
    
    //Grid on Host
    Grid grid;

    //set dimensions of the grid
    grid.dim.x=2;
    grid.dim.y=3;
    grid.dim.z=4;
    
    //build  and fill grid on Host
    buildGrid(&grid, aabbs, numT);
    
    
    
    cout<<"\n\nBoundingBox Grid:\nxmin,ymin,zmin: "
        <<grid.aabb.min.x<<" ,"<<grid.aabb.min.y<<" ,"<<grid.aabb.min.z
        <<"\nxmax,ymax,zmax: "
        <<grid.aabb.max.x<<" ,"<<grid.aabb.max.y<<" ,"<<grid.aabb.max.z<<"\n";
    
    cout<<"Zellenlaenge in x: "<<(grid.aabb.max.x-grid.aabb.min.x)/grid.dim.x<<endl;
    cout<<"Zellenlaenge in y: "<<(grid.aabb.max.y-grid.aabb.min.y)/grid.dim.y<<endl;
    cout<<"Zellenlaenge in z: "<<(grid.aabb.max.z-grid.aabb.min.z)/grid.dim.z<<endl;
    
    int numCells = grid.dim.x * grid.dim.y * grid.dim.z;
    //Ausgabe der Zellen und der zugehoerigen Listen
    for(int cellIdx=0;cellIdx<numCells;++cellIdx){
        int i,j,k;
        cellIdx2ijk(&grid,cellIdx,&i,&j,&k);
        cout<<"cellIdx="<<cellIdx<<" ->i="<<i<<" ->j="<<j<<" ->k="<<k<<"\n"
            <<"Anzahl Elem: "<<grid.cellList[cellIdx].length<<" Indizes: ";
        for(int i2=0;i2<grid.cellList[cellIdx].length;++i2)
            cout<<grid.cellList[cellIdx].prismIdxList[i2]<<" ";
        cout<<endl<<endl;
    }

    Grid d_grid;
    copyGrid2GPU(&grid, &d_grid);

    int *results;
    int *d_results;
    int size;
    int *d_size;
    CUDA_CHECK(cudaMalloc((void**) &d_size, sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**) &d_results, sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_size, &size, sizeof(int), cudaMemcpyHostToDevice));

    CUDA_CHECK_KERNEL(testkernel<<<1,1>>>(d_grid, d_results, d_size));

    CUDA_CHECK(cudaMemcpy(&size, d_size, sizeof(int), cudaMemcpyDeviceToHost));
    //results = (int*) malloc(size*sizeof(int));
    results = new int[size];
    CUDA_CHECK(cudaMemcpy(results, d_results, size*sizeof(int), cudaMemcpyDeviceToHost));

    printf("Got %i results: ", size);
    for(int i=0; i<size; ++i) {
      printf("%i, ", results[i]);
    }
    puts("");

    return 0;
}

