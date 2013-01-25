
#include "buildgrid.h"

void prepareGrid(Grid *d_grid, int4 dimGrid, const PrismCu *prisms, int numPrisms)
{
	
	//Grid on Host
	Grid grid;
	
	//set dimensions of the grid
	grid.dim.x=dimGrid.x;
	grid.dim.y=dimGrid.y;
	grid.dim.z=dimGrid.z;
	
	
	//array with BoudningBoxes for prisms
	Aabb* aabbs=new Aabb[numPrisms];
	
	//calc AABBs of the prisms
	calcAabbs(prisms, aabbs, numPrisms);
	
	//build  and fill grid on Host
	buildGrid(&grid, aabbs, numPrisms);
	
	//copy grid to GPU
	copyGrid2GPU(&grid,d_grid);
	
	delete [] aabbs;
}
