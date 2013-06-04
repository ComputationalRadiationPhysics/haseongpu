//#include "ray_propagation_gpu.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include <vector>
#include "curand_kernel.h"
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>

#define TEST_VALUES true
#define USE_IMPORTANCE false
#define DIVIDE_PI true
#define SMALL 1E-06
#define VERY_SMALL 1E-14

#define CUDA_CHECK_RETURN(value) {				\
  cudaError_t _m_cudaStat = value;				\
  if (_m_cudaStat != cudaSuccess) {				\
    fprintf(stderr, "Error %s at line %d in file %s\n",	\
    cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);	\
    exit(1);							\
  }								\
}
#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x) != CURAND_STATUS_SUCCESS) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

__device__ double cladAbsorption;
__device__ double nTot;
__device__ double sigmaE;
__device__ double sigmaA;
__device__ double thicknessOfPrism;
__device__ int numberOfLevels;
__device__ int cladNumber;
__device__ int numberOfPoints;
__device__ int numberOfTriangles;

__device__ double propagationOld(
		double x_pos, 
		double y_pos, 
		double z_pos, 
		double x_dest, 
		double y_dest, 
		double z_dest, 
		int t_start, 
		int mesh_start, 
		double *p_in,
		double *n_x,
		double *n_y,
		int *n_p,
		int *neighbors,
		int *forbidden,
		int *cell_type,
		double *beta_v){
//    in first try no reflections
//    calculate the vector and make the the calculation, which surface would be the shortest to reach
//    then get the length, make the integration, get the information about the next cell out of the array
//    set the point to the surface (this surface is "forbidden" in the calculations)
//    proceed until you hit a the point or the surface
//    if you are closer then "small" stop and return the value
    double vec_x, vec_y,vec_z, norm;
    double distance, length, length_help, distance_total;
    double gain=1;
    double nominator, denominator;
    int tri, cell_z; // the current triangle number and position concerning the z's
    int decider; // which one is the shortest - info
    int tri_next, cell_z_next, forb, forb_dump;
	int ct; // used to read out cell_type 
	int size_p = numberOfPoints;
	int N_cells = numberOfTriangles;
	double z_mesh = thicknessOfPrism;
	double sigma_a = sigmaA;
	double sigma_e = sigmaE;
	double clad_num = cladNumber;
	double clad_abs = cladAbsorption;
	double N_tot = nTot;
	
    
    
//    initial positions
    tri = t_start;
    cell_z = mesh_start;
    
    
//    definition of the vectors without reflections
    vec_x = (x_dest - x_pos);
    vec_y = (y_dest - y_pos);
    vec_z = (z_dest - z_pos);
    
    norm = sqrt(vec_x*vec_x+vec_y*vec_y+vec_z*vec_z);
    
    vec_x = vec_x/norm;
    vec_y = vec_y/norm;
    vec_z = vec_z/norm;
    
//    now calculate the length to travel
    distance = sqrt((x_dest - x_pos)*(x_dest - x_pos)+(y_dest - y_pos)*(y_dest - y_pos)+(z_dest - z_pos)*(z_dest - z_pos));
    distance_total = distance;
	// does this make sense?
    length = distance;
    
    forb = -1;
	
//	mexPrintf("Propagation called");
//    mexEvalString("drawnow;");
    
//    the ray has to be set to be ALIVE before!
//    now do the unlimited for loop - break!!!
    for(;;)
    {
	
//	  mexPrintf("Propagation for part called\n\n");
//    mexEvalString("drawnow;");
//        definition for decider
//        0,1,2: int for the neighbors
//        3: hor plane up
//        4: hor plane down
        
//        at first set the decider = -1;
        decider = -1;
        length = distance;
		
		
//		  read, which type of cell it is you are propagation in
		ct = cell_type[tri];
        
//        mexPrintf("forb: %i\n",forb);
//        mexEvalString("drawnow;");
        
//        try the triangle faces
//        remember the correlation between the normals and the points
//        n1: p1-2, n2: p1-3, n3:p2-3
//        the third coordinate (z) of the particpating points for the surfaces can be set to be z=0, 
//        as everything uses triangular "tubes/prisms", as well as n_z=0 in this case!
        if (forb != 0){
            nominator = (n_x[tri]*p_in[n_p[tri]] + n_y[tri]*p_in[n_p[tri]+size_p]) - (n_x[tri]*x_pos + n_y[tri]*y_pos);
            denominator = n_x[tri]*vec_x + n_y[tri]*vec_y;
            if (denominator != 0.0)
            {
                length_help = nominator/denominator;
                if (length_help < length && length_help > VERY_SMALL)
                {
                    length = length_help;
                    decider = 0;
                    forb_dump = (forbidden[tri]);
                }
            }
        }
        
        if (forb != 1){
            nominator = (n_x[tri+N_cells]*p_in[n_p[tri+N_cells]] + n_y[tri+N_cells]*p_in[n_p[tri+N_cells]+size_p]) - (n_x[tri+N_cells]*x_pos + n_y[tri+N_cells]*y_pos);
            denominator = n_x[tri+N_cells]*vec_x + n_y[tri+N_cells]*vec_y;
            if (denominator != 0.0)
            {
                length_help = nominator/denominator;
                if (length_help < length && length_help > VERY_SMALL)
                {
                    length = length_help;
                    decider = 1;
                    forb_dump = (forbidden[tri+N_cells]);
                }
            }
        }
        
        if (forb !=2){
            nominator = (n_x[tri+2*N_cells]*p_in[n_p[tri+2*N_cells]] + n_y[tri+2*N_cells]*p_in[n_p[tri+2*N_cells]+size_p]) - (n_x[tri+2*N_cells]*x_pos + n_y[tri+2*N_cells]*y_pos);
            denominator = n_x[tri+2*N_cells]*vec_x + n_y[tri+2*N_cells]*vec_y;
            if (denominator != 0.0)
            {
                length_help = nominator/denominator;
                if (length_help < length && length_help > VERY_SMALL)
                {
                    length = length_help;
                    decider = 2;
                    forb_dump = (forbidden[tri+2*N_cells]);
                }
            }
        }
        
//        try the horizontal planes, which one is the shortest, n_x and n_y are zero!, n_z =1!
//        at first the upper plane
        if (forb != 3){
            nominator = (cell_z+1)*z_mesh - z_pos;
            denominator = z_pos*vec_z;
            if (denominator != 0.0)
            {
                length_help = nominator/denominator;
                if (length_help < length && length_help > VERY_SMALL)
                {
                    length = length_help;
                    decider = 3;
                    forb_dump = 4; // you are not allowed to go down in the next step
                }
            }
        }
        
//        next is the lower plane
        if (forb != 4){
            nominator = (cell_z)*z_mesh - z_pos;
            denominator = z_pos*vec_z;
            
            if (denominator != 0.0)
            {
                length_help = nominator/denominator;
                if (length_help < length && length_help > VERY_SMALL)
                {
                    length = length_help;
                    decider = 4;
                    forb_dump = 3; // you are not allowed to go up in the next step
                }
            }
        }
        
        forb = forb_dump;
		
        
//        now make a switch to differ the different cases
        switch(decider){
                
            case 0:
//                this is the case for the intersection with the first choice triangle-surface
                tri_next = neighbors[tri];
                cell_z_next = cell_z;
                break;
                
            case 1:
//                second triangle surface
                tri_next = neighbors[tri+N_cells];
                cell_z_next = cell_z;
                break;
                
            case 2:
//                third triangle surface
                tri_next = neighbors[tri+2*N_cells];
                cell_z_next = cell_z;
                break;
                
            case 3:
//                go one plane up
                tri_next = tri;
                cell_z_next = cell_z + 1;
                break;
                
            case 4:
//                go one plane down
                tri_next = tri;
                cell_z_next = cell_z - 1;
                break;
                
            default:
//                make an error statement
                break;
        }
        
//        now we know where to go, let's make the integration
//        take the beta_v[tri+cell_z*N_cells] 

//		  at this position do the decision whether it is a gain part or cladding
//		  it might be absorbing or amplifying, for the cladding only absorbing
//		  a simple "if then"

		if (ct == clad_num){
			gain = gain * exp(-clad_abs * length);
		}
		else {
			gain = gain * exp(N_tot*(beta_v[tri+cell_z*N_cells]*(sigma_e + sigma_a)-sigma_a)*length);
		}
//        gain = LineIntegralMCRK4_S(3, tri, cell_z, gain, length);
        
//        after integration make the propagation
        
//        mexPrintf("Distance: %f, Length: %f\n",distance, length);
//        mexPrintf("decider: %i, forbidden: %i\n",decider, forb);
//        mexPrintf("vec_x: %f, vec_y: %f, vec_z: %f\n", vec_x, vec_y, vec_z);
//        mexPrintf("current_x: %f current_y: %f current_z: %f\n", x_pos, y_pos, z_pos);
//        mexPrintf("tri: %i, tri_next: %i, cell_z: %i, cell_next: %i\n", tri, tri_next, cell_z, cell_z_next);
//        mexEvalString("drawnow;");
//        str=mxCreateString("Press a key");
//        mexCallMATLAB(1,&dump,1,&str,"input"); 
//        str and dump should be defined to be a *mxArray and don't forget to kill them at the end
        
        distance -= length;
        
//        return 1;
//        
        
        x_pos = x_pos + length*vec_x;
        y_pos = y_pos + length*vec_y;
        z_pos = z_pos + length*vec_z;
        
        if (abs(distance)< SMALL)
        {
            break;
        }
        
        
//        now set the next cell
        tri = tri_next;
        cell_z = cell_z_next;      
        
//        break;
//        now we should make the integration routine
    }
    
    gain /= (distance_total*distance_total);

    return gain;
}

/**
 * @brief Propagate a ray between 2 points and calculate the resulting ASE-Flux at the Destination
 *
 * @params x_pos		the x-coordinate where the ray starts
 *         y_pos		the y-coordinate where the ray starts
 *         z_pos		the z-coordinate where the ray starts
 *         x_dest		the destination of the ray (x-coordinate)
 *         y_dest		the destination of the ray (y-coordinate)
 *         z_dest		the destination of the ray (z-coordinate)
 *         t_start		the index of the triangle, in which the ray starts
 *         mesh_start	the level of the mesh (the slice) in which the ray starts
 *		   p_in			coordinates of the sample-points of one layer (first all x-coordinates, then all y-coordinates)
 *		   n_x			x-coordinates for the normal-vectors for the 3 rectangular sides of each prism
 *		   n_y			y-coordinates for the normal-vectors for the 3 rectangular sides of each prism
 *		   n_p			indices of the points where the normal-vectors start	
 *		   neighbors	indices of the adjacent triangles	
 *		   forbidden	sides of the new triangles which are forbidden, after coming from an adjacent triangle
 *		   cell_type	contains the material-constant for each cell/prism
 *		   beta_v		contains the beta-values for each cell/prism
 *
 */
__device__ double rayPropagationGpu(
		double xPos,
		double yPos,
		double zPos,
		double xDestination,
		double yDestination,
		double zDestination,
		int firstTriangle,
		int firstLevel,
		double *points,
		double *xOfNormals,
		double *yOfNormals,
		int *positionsOfNormalVectors,
		int *neighbors,
		int *forbidden,
		int* cellTypes,
		double* betaValues){
	//    no reflections
	//
	//    create the vector between both points and calculate which surface would be the shortest to reach.
	//    then get the length, make the integration, get the information about the next cell out of the array
	//    set the point to the surface (this surface is "forbidden" in the calculations)
	//    proceed until you hit a the point or the surface
	
	double vec_x, vec_y,vec_z;
	double distanceRemaining, length, lengthHelp, distanceTotal;
	double nominator, denominator;
	double gain=1.;
	int tri, cell_z; // the current triangle number and position concerning the z's
	int tri_next, cell_z_next, forb, forb_dump;
	int offset;
#if TEST_VALUES==true
	double testDistance = 0;
#endif


	//    initial positions
	tri = firstTriangle;
	cell_z = firstLevel;

	// direction-vector (without reflections)
	vec_x = (xDestination - xPos);
	vec_y = (yDestination - yPos);
	vec_z = (zDestination - zPos);

	// total distance to travel
	distanceTotal = sqrt(vec_x*vec_x+vec_y*vec_y+vec_z*vec_z);
	// normalized direction-vector
	vec_x = vec_x/distanceTotal;
	vec_y = vec_y/distanceTotal;
	vec_z = vec_z/distanceTotal;

	// remaining distance to travel
	distanceRemaining = distanceTotal;

	// at the beginning, all surfaces are possible
	forb = -1;
	int loopbreaker=0;
	for(;;)
	{
		// the length of the ray-part inside the current prism. We try to minimize this value
		length = distanceRemaining;
		//        definition for decider
		//        0,1,2: int for the neighbors
		//        3: hor plane up
		//        4: hor plane down
		//        try the triangle faces
		//        remember the correlation between the normals and the points
		//        n1: p1-2, n2: p1-3, n3:p2-3
		//        the third coordinate (z) of the particpating points for the surfaces can be set to be z=0, 
		//        as everything uses triangular "prisms", as well as n_z=0 in this case!
		
		// forb describes the surface, from which the ray enters the prism.
		// this surface is no suitable candidate, since the length would be 0!
		if (forb != 0){
			denominator = xOfNormals[tri]*vec_x + yOfNormals[tri]*vec_y;
			// see if we intersect at all
			if (denominator != 0.0)
			{
				nominator = (xOfNormals[tri]*points[positionsOfNormalVectors[tri]] + yOfNormals[tri]*points[positionsOfNormalVectors[tri]+ numberOfPoints]) - (xOfNormals[tri]*xPos + yOfNormals[tri]*yPos);
				lengthHelp = nominator/denominator;
				// if we found a new smallest length, use it
				if (lengthHelp < length && lengthHelp > VERY_SMALL)
				{
					length = lengthHelp;
					forb_dump = (forbidden[tri]);	
					tri_next = neighbors[tri];
					cell_z_next = cell_z;

				}
			}
		}

		// see forb !=0 case
		if (forb != 1){
			//offset, since the 3 rectangular surfaces are stored at different positions in the array
			offset = tri+numberOfTriangles;
			denominator = xOfNormals[offset]*vec_x + yOfNormals[offset]*vec_y;
			if (denominator != 0.0)
			{
				nominator = (xOfNormals[offset]*points[positionsOfNormalVectors[offset]] + yOfNormals[offset]*points[positionsOfNormalVectors[offset]+ numberOfPoints]) - (xOfNormals[offset]*xPos + yOfNormals[offset]*yPos);
				lengthHelp = nominator/denominator;
				if (lengthHelp < length && lengthHelp > VERY_SMALL)
				{
					length = lengthHelp;
					forb_dump = (forbidden[offset]);
					tri_next = neighbors[offset];
					cell_z_next = cell_z;
				}
			}
		}

		// see forb !=0 case
		if (forb !=2){
			offset = tri+2*numberOfTriangles;
		denominator = xOfNormals[offset]*vec_x + yOfNormals[offset]*vec_y;
			if (denominator != 0.0)
			{
				nominator = (xOfNormals[offset]*points[positionsOfNormalVectors[offset]] + yOfNormals[offset]*points[positionsOfNormalVectors[offset]+ numberOfPoints]) - (xOfNormals[offset]*xPos + yOfNormals[offset]*yPos);
				lengthHelp = nominator/denominator;
				if (lengthHelp < length && lengthHelp > VERY_SMALL)
				{
					length = lengthHelp;
					forb_dump = (forbidden[offset]);
					tri_next = neighbors[offset];
					cell_z_next = cell_z;
				}
			}
		}

		// if-structure "optimized"
		denominator = zPos*vec_z;
		if (denominator != 0.0){
			if (forb != 3){
				{
					nominator = (cell_z+1)* thicknessOfPrism - zPos;
					lengthHelp = nominator/denominator;
					if (lengthHelp < length && lengthHelp > 0.0)
					{
						length = lengthHelp;
						//decider = 3;
						forb_dump = 4; // you are not allowed to go down in the next step
						tri_next = tri;
						cell_z_next = cell_z + 1;
					}
				}
			}

			// next is the lower plane
			if (forb != 4){
				nominator = (cell_z)* thicknessOfPrism - zPos;
				lengthHelp = nominator/denominator;
				if (lengthHelp < length && lengthHelp > 0.0)
				{
					length = lengthHelp;
					//decider = 4;
					forb_dump = 3; // you are not allowed to go up in the next step
					tri_next = tri;
					cell_z_next = cell_z - 1;
				}
			}
		}

		// set the new forbidden surface
		forb = forb_dump;

		// switch is now directly included into the if-statements
		//
		// now we know where to go, let's make the integration
		// take the betaValues[tri+cell_z*N_cells]
		//
		// at this position do the decision whether it is a gain part or cladding
		// it might be absorbing or amplifying, for the cladding only absorbing
		// a simple "if then"

		if (cellTypes[tri] == cladNumber){
			gain *= exp((-1)*(cladAbsorption * length));
		}
		else {
			gain *= (double) exp(length* nTot*(betaValues[tri+cell_z*numberOfTriangles]*(sigmaE + sigmaA)-sigmaA));
		}


		// the remaining distance is decreased by the length we travelled through the prism
		distanceRemaining -= length;
		

#if TEST_VALUES==true
		testDistance += length;
		if(loopbreaker>500){
			printf("Loopbreaker reached. firstTriangle: %d, level: %d, length: %f, distanceTotal:%f, testDistance%f, distanceRemaining:%f\n",firstTriangle,firstLevel,length,distanceTotal,testDistance,distanceRemaining);
			return 0.;
		}else{
			loopbreaker++;
		}
#endif
		// if the distance between the destination and our current position is small enough, we are done
		if (fabs(distanceRemaining) < SMALL){
			break;
		}

		// now set the next cell and position
		xPos = xPos + length*vec_x;
		yPos = yPos + length*vec_y;
		zPos = zPos + length*vec_z;

		tri = tri_next;
		cell_z = cell_z_next;      

	}

#if TEST_VALUES==true
	if(fabs(distanceTotal-testDistance) > SMALL)
		printf("Distance too big! firstTriangle: %d, level: %d, length: %f, distanceTotal:%f, testDistance%f, distanceRemaining:%f\n",firstTriangle,firstLevel,length,distanceTotal,testDistance,distanceRemaining);
#endif
	
	gain /=(distanceTotal*distanceTotal);
	return gain;
}

#if USE_IMPORTANCE==true
__device__ void importf(curandState localstate, int point, int mesh_start, double *importance, int *N_rays, double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int *forbidden, int *cell_type, double *beta_v, double *center_x, double *center_y, int *surface, int NumRays)
{
    int i_t, i_z, Rays_dump=0, rays_left, i_r, rand_t, rand_z;
    double sum_phi=0.0, surf_tot=0.0;
    double distance, x_pos, y_pos, z_pos;
    double prop;

//    calculate the gain from the centers of each of the boxes to the observed point
//    calculate the gain and make a "mapping"
//    receipt: pick the point in the center of one cell, 
//    calculate the gain from this point to the observed point,
//    estimate the inner part of the Phi_ASE - Integral,
//    scale the amount of rays proportionally with it
//    sum the amount of rays and scale it to Int=1, which gives the inverse weights
//    the number of rays is determined via floor(), with ceil(), zero-redions could be added

//    use the routine "propagation"!, test: no reflections, just exponential
    x_pos = p_in[point];
    y_pos = p_in[point+numberOfPoints];
    z_pos = mesh_start * thicknessOfPrism;

    for (i_t=0;i_t<numberOfTriangles;i_t++)
    {

        for (i_z=0;i_z<(numberOfLevels-1);i_z++) //remember the definition differences MatLab/C for indices
        {
//            at this point replace the following routine with propagation(...)
//            later expand this with the beta/tau values...
			prop = naivePropagation(center_x[i_t], center_y[i_t], thicknessOfPrism*(i_z+0.5),  x_pos, y_pos, z_pos, i_t, i_z, p_in, n_x, n_y, n_p, neighbors, forbidden , cell_type, beta_v);
			// Propagation vom Zentrum jedes Prismas zu jedem Samplepunkt 
			//
//			prop = propagation(center_x[i_t], center_y[i_t], z_mesh*(i_z+0.5), x_pos, y_pos, z_pos, i_t, i_z);
            importance[i_t + i_z*numberOfTriangles] = beta_v[i_t+i_z*numberOfTriangles]*(prop);
            sum_phi += importance[i_t + i_z*numberOfTriangles];

        }
        surf_tot += surface[i_t];
    }

//    now calculate the number of rays
    for (i_t=0;i_t<numberOfTriangles;i_t++)
    {
        for (i_z=0;i_z<(numberOfLevels-1);i_z++) //remember the definition differences MatLab/C for indices
        {
//            this is the amount of the sampled rays out of the cells
            N_rays[i_t + i_z*numberOfTriangles] = (int)(floor(importance[i_t + i_z*numberOfTriangles]/sum_phi*NumRays));

            Rays_dump +=  N_rays[i_t + i_z*numberOfTriangles];
        }
    }

    rays_left = NumRays-Rays_dump;
//    distribute the remaining not distributed rays randomly
    if ((rays_left)>0)
    {
        for (i_r=0;i_r<rays_left;i_r++)
        {
            rand_t = (int )(curand_uniform(&localstate)*numberOfTriangles);
            rand_z = (int )(curand_uniform(&localstate)*(numberOfLevels-1));
            N_rays[rand_t + rand_z*N_cells]++;
        }
    }

//    now think about the mount of rays which would come out of this volume(surface)
//    dividing this number with the new amount of rays gives the final importance weight for this area!
    for (i_t=0;i_t<N_cells;i_t++)
    {
        for (i_z=0;i_z<(numberOfLevels-1);i_z++) //remember the definition differences MatLab/C for indices
        {
//            this is the amount of the sampled rays out of the cells
            if (N_rays[i_t + i_z*numberOfTriangles]>0)
            {
                importance[i_t + i_z*numberOfTriangles] = NumRays*surface[i_t]/surf_tot/N_rays[i_t + i_z*N_cells];
//                importance[i_t + i_z*numberOfTriangles] = NumRays*surface[i_t]/surf_tot;
            }
            else
            {
                importance[i_t + i_z*numberOfTriangles] = 0; // case of beta of this point == 0 e.g.
            }
        }
    }
}
#endif


/**
 * Initializes the global variables of the GPU with the correct values.
 * All those values are from the original propagation-function which we ported.
 */
__global__ void setupGlobalVariablesKernel ( 
		double host_sigma_e,
		double host_sigma_a, 
		int host_clad_num, 
		double host_clad_abs, 
		double host_N_tot, 
		int host_N_cells, 
		double host_z_mesh, 
		int host_mesh_z, 
		int hostNumberOfPoints )
{
	sigmaE = host_sigma_e;	
	sigmaA = host_sigma_a;
	cladNumber = host_clad_num;
	cladAbsorption = host_clad_abs;
	nTot = host_N_tot;
	numberOfTriangles = host_N_cells;
	thicknessOfPrism = host_z_mesh;
	numberOfLevels = host_mesh_z;
	numberOfPoints = hostNumberOfPoints;
	//printf("Sigma_e in setup=%f\tSigma_eHost=%f\n",sigma_e,host_sigma_e);
} 

#if USE_IMPORTANCE==true
__global__ void importanceKernel( curandState *globalState, double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int *forbidden, int* cell_type, int hostNumberOfTriangles, double* beta_v, double *importance, int *N_rays, double *center_x, double *center_y, int *surface, int NumRays) {
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	for(int i=0; i< hostNumberOfTriangles; ++i){
		for(int j=0; j< numberOfLevels; ++j){
			importf(globalState[id], i,j, importance, N_rays, p_in, n_x, n_y, n_p, neighbors, forbidden, cell_type, beta_v, center_x, center_y,surface, NumRays);

		}
	}
	
}
#endif

/**
 * Prints some of the global device variables.
 * Is only used for testing
 */
__global__ void testKernel (  ){
	printf("\nSigmaE=%.6e",sigmaE);
	printf("\nSigmaA=%.6e",sigmaA);
	printf("\nNumberOfLevels=%d",numberOfLevels);
	printf("\nNumberOfPoints=%d",numberOfPoints);
	printf("\nthicknessOfPrism_=%.6e",thicknessOfPrism);
	printf("\nnumberOfTriangles=%d",numberOfTriangles);
	printf("\nnTot=%.6e",nTot);
} 


/**
 * Does the raytracing for a single Sample point (in a defined level).
 * This Kernel has to be started for each sample point with the same value for iterations
 * and the same number of blocks/threads.
 *
 * \var globalState the state of the mersenneTwister PRNG
 * 		(has a maximum of 200 positions!)
 * \var phi points to a memory region which is initialized with 0
 * 		(can hold one value for each sample point)
 * \var point2D the index of the current sample point (points to p_in)
 * \var level the level of the current sample point (how deep we are through the material)
 * \var iterations the number rays which are computed by this thread
 * 		(always for the same combination of startprism+samplepoint
 */
#if USE_IMPORTANCE==true
__global__ void raytraceStep( curandStateMtgp32* globalState, float* phi, int point2D, int level, int iterations, double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int *forbidden, int* t_in, int* cell_type, int hostNumberOfTriangles, double* beta_v, double *importance) {
#else
__global__ void raytraceStep(
		curandStateMtgp32* globalState, 
		float* phi, 
		const int point2D, 
		const int level, 
		const int iterations, 
		double *p_in, 
		double *n_x, 
		double *n_y, 
		int *n_p, 
		int *neighbors, 
		int *forbidden, 
		int* t_in, 
		int* cell_type, 
		double* beta_v) {
#endif
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	const unsigned numberOfPrisms = (numberOfTriangles * (numberOfLevels-1));
	const unsigned threadsPerPrism = blockDim.x * gridDim.x/numberOfPrisms;
	// break, if we have more threads than we need
	if(id >= threadsPerPrism*numberOfPrisms)
		return;

	double gain = 0.;
	const int endPointX = p_in[point2D];
	const int endPointY = p_in[ numberOfPoints + point2D];
	const int endPointZ = level* thicknessOfPrism;


	// this should give the same start values multiple times (so that every thread uses the same prism, which yields
	// big benefits for the memory access (and caching!)
	unsigned startPrism = id/threadsPerPrism;
	int startLevel = (startPrism)/numberOfTriangles;
	int startTriangle=(startPrism-(numberOfTriangles*startLevel));

	// the indices of the vertices of the starttriangle
	int t1 = t_in[startTriangle];
	int t2 = t_in[startTriangle+ numberOfTriangles];
	int t3 = t_in[startTriangle+2*numberOfTriangles];

	// do all this multiple times (we can't have more than 200 blocks due to restrictions of the Mersenne Twister)
	for (int i=0; i<iterations ; ++i){
		// random startpoint generation
		double u = curand_uniform(&globalState[blockIdx.x]);
		double v = curand_uniform(&globalState[blockIdx.x]);

		if((u+v)>1)
		{
			u = 1-u;
			v = 1-v;
		}
		double w = 1-u-v;

		// convert the random startpoint into coordinates
		double zRand = (startLevel + curand_uniform(&globalState[blockIdx.x]))* thicknessOfPrism;
		double xRand = p_in[t1]*u + p_in[t2]*v + p_in[t3]*w;
		double yRand = p_in[ numberOfPoints + t1]*u + p_in[ numberOfPoints + t2]*v + p_in[ numberOfPoints + t3]*w;

		__syncthreads();
		gain += double(rayPropagationGpu(xRand, yRand, zRand, endPointX, endPointY, endPointZ, startTriangle, startLevel ,p_in, n_x, n_y, n_p, neighbors, forbidden , cell_type, beta_v));
	//	gain += double(propagationOld(xRand, yRand, zRand, endPointX, endPointY, endPointZ, startTriangle, startLevel ,p_in, n_x, n_y, n_p, neighbors, forbidden , cell_type, beta_v));
	}
	

	// do the multiplication just at the end of all iterations
	// (gives better numeric behaviour)
	gain *= beta_v[startTriangle + numberOfTriangles*startLevel];
#if USE_IMPORTANCE==true
	atomicAdd(&(phi[point2D + level*numberOfPoints]),float(gain*importance[startTriangle + numberOfTriangles*startLevel]));
#else
	atomicAdd(&(phi[point2D + level*numberOfPoints]),float(gain)); 
#endif
	return;
}

double* doubleVectorToArray(std::vector<double> *input){
	double* output;
	output = (double*) malloc(sizeof(double) * input->size());
	for(int i=0; i< input->size(); ++i){
		output[i] = input->at(i);	
	}
	return output;
}
int* intVectorToArray(std::vector<int> *input){
	int* output;
	output = (int*) malloc(sizeof(int) * input->size());
	for(int i=0; i< input->size(); ++i){
		output[i] = input->at(i);	
	}
	return output;
}
unsigned* unsignedVectorToArray(std::vector<unsigned> *input){
	unsigned* output;
	output = (unsigned*) malloc(sizeof(unsigned) * input->size());
	for(int i=0; i< input->size(); ++i){
		output[i] = input->at(i);	
	}
	return output;
}

//----------------------------------------------------
// Host Code
//----------------------------------------------------
float runRayPropagationGpu(
		std::vector<double> *ase, 
		unsigned &threads, 
		unsigned &blocks, 
		unsigned &totalNumberOfRays,
		std::vector<double> *betaValuesVector,
		std::vector<double> *xOfNormalsVector,
		std::vector<double> *yOfNormalsVector,
		std::vector<unsigned> *cellTypesVector,
		std::vector<unsigned> *triangleIndicesVector,
		std::vector<int> *forbiddenVector,
		std::vector<int> *neighborsVector,
		std::vector<int> *positionsOfNormalVectorsVector,
		std::vector<double> *pointsVector,
		std::vector<double> *betaCellsVector,
		float hostCladAbsorption,
		int hostCladNumber,
		float hostNTot,
		float hostSigmaA,
		float hostSigmaE,
		unsigned hostNumberOfPoints,
		unsigned hostNumberOfTriangles,
		unsigned hostNumberOfLevels,
		float hostThicknessOfPrism,
		float hostCrystalFluorescence)
{
	/** GPU Kernel Variables
	 * The idea is, that the number of threads is fixed (to maximize GPU occupancy)
	 * and the number of blocks as well (200 is the maximum for the standard
	 * Mersenne Twister implementaion). Therefore, the number of rays per sample
	 * are fixed to be k*200*256.
	 * That means, sometimes we have to increase the number of rays a little.
	 *
	 * \var iterations is used to give every thread k iterations (to simulate k rays)
	 *
	 * note that every samplepoint receives the exact same number of rays.
	 */
	
	//
	
	fprintf(stderr, "\nConverting Vectors to Arrays\n");

	double* hostBetaValues = doubleVectorToArray(betaValuesVector);
	double* hostXOfNormals = doubleVectorToArray(xOfNormalsVector);
	double* hostYOfNormals = doubleVectorToArray(xOfNormalsVector);
	unsigned* hostCellTypes = unsignedVectorToArray(cellTypesVector);
	unsigned* hostTriangleIndices = unsignedVectorToArray(triangleIndicesVector);
	int* hostForbidden = intVectorToArray(forbiddenVector);
	int* hostNeighbors = intVectorToArray(neighborsVector);
	int* hostPositionsOfNormalVectors = intVectorToArray(positionsOfNormalVectorsVector);
	double* hostPoints = doubleVectorToArray(pointsVector);

	
	fprintf(stderr, "\nCalculating optimal number of Rays\n");
	threads = 256;
	blocks = 200;
	unsigned raysPerSample = ceil(totalNumberOfRays/double(hostNumberOfPoints * (hostNumberOfLevels)));
	unsigned hostNumberOfPrisms = (hostNumberOfTriangles * (hostNumberOfLevels-1));
	unsigned hostThreadsPerPrism = threads*blocks/hostNumberOfPrisms; //9
	unsigned raysPerIterationPerSample = hostThreadsPerPrism*hostNumberOfPrisms; //48600
	int iterations = ceil(double(raysPerSample) / raysPerIterationPerSample);
	//fprintf(stderr, "raysPerSample=%d\n",raysPerSample);
	//fprintf(stderr, "interations=%d\n",iterations);
	raysPerSample = raysPerIterationPerSample * iterations;
	totalNumberOfRays = unsigned(raysPerSample * (hostNumberOfPoints * (hostNumberOfLevels)));
	//fprintf(stderr, "After Normalization:\n");
	//fprintf(stderr, "raysPerSample=%d\n",raysPerSample);
	//fprintf(stderr, "totalNumberOfRays=%d\n",totalNumberOfRays);

	fprintf(stderr, "(%d per Sample)\n",raysPerSample);
	

	/** Variables for the device
	 * These are on-GPU representations of the input parameters
	 * of variable size.
	 *
	 * \var p_in: coordinates of the sample-points of one layer (first all x-coordinates, then all y-coordinates)
	 * \var n_*: values of the normal-vectors for the 3 rectangular sides of each prism (described in 2D)
	 * \var beta_v: the beta values of the prisms
	 * \var phi: the accumulated ASE-Flux for each sample point
	 * \var forbidden: the side of the prism through which the ray "entered" the prism
	 * \var n_p: the points where the normals (n_x,n_y) start
	 * \var neighbors: indices to the adjacent triangles in t_in
	 * \var t_in: indices of the points which are considered to be a triangle (A points start from 0, B points from size_t, C points from size_t*2)
	 * \var cell_type: determines which cell type we are looking at.
	 * other input parameters are put to the GPU by the setupGlobalVariablesKernel
	 */
	double  *points, *xOfNormals, *yOfNormals, *betaValues;
	float *phi;
	int *forbidden, *positionsOfNormalVectors, *neighbors, *triangleIndices, *cellTypes;


	/** Variables for benchmarking and results
	 * \var hostPhi will contain the ASE-Flux integral (not normalized)
	 */
	float runtimeGpu = 0.0;
	fprintf(stderr, "Initializing Timer\n");
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	fprintf(stderr, "Initializing ASE-Values with 0\n");
	float hostPhi[hostNumberOfPoints * (hostNumberOfLevels)];
	
	for(int i=0;i<hostNumberOfPoints*(hostNumberOfLevels);++i){
			hostPhi[i] = 0.;
	}


	/** The Mersenne Twister PRNG
	 * These calls are more or less directly from the curand
	 * \var devMTGPStates the current randomness states of each block (may not exceed 200!)
	 *
	 */
	curandStateMtgp32 *devMTGPStates;
	{
		fprintf(stderr, "\nInitializing Mersenne Twister PRNG\n");
		/**Allocate space for PRNG states on device */
		CUDA_CALL(cudaMalloc((void **)&devMTGPStates, blocks * sizeof(curandStateMtgp32)));

		/** Allocate space for MTGP kernel parameters */
		mtgp32_kernel_params *devKernelParams;
		CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));

		/**Reformat from predefined parameter sets to kernel format,
		 * and copy kernel parameters to device memory */
		CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));

		/** Initialize one state per thread block */
		/** \TODO initialize with time */
		CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, blocks, 1234));
	}

/** Parameters for Importance Sampling
 * These are conditionally compiled
 * \var devStates is for the random number generator (not the Mersenne Twister!)
 * \var center_* describes the center of a triangle
 */
#if USE_IMPORTANCE==true
	curandState *devStates;
	int *N_rays, *surface;
	double *center_x, *center_y, *importance;
	double host_importance[hostNumberOfPoints * (hostNumberOfLevels)];
	int host_N_rays[hostNumberOfTriangles * hostNumberOfLevels];
	for(int i=0;i<hostNumberOfPoints*(hostNumberOfLevels);++i){
		host_importance[i] = 0.;
	}
	for(int i=0;i<hostNumberOfTriangles*hostNumberOfLevels;++i){
		host_N_rays[i] = 0;
	}
#endif



	/** Allocation of memory on the GPU and setting of global GPU-variables
	 *
	 */
	{
		fprintf(stderr, "\nFilling the device Variables\n");
		//Create constant values on GPU
		setupGlobalVariablesKernel<<<1,1>>>(
				double(hostSigmaE), 
				double(hostSigmaA),
				hostCladNumber,
				double(hostCladAbsorption),
				double(hostNTot), 
				hostNumberOfTriangles, 
				double(hostThicknessOfPrism),
				hostNumberOfLevels, 
				hostNumberOfPoints); //@OPTIMIZE: initialize the constants as constants...

		cudaThreadSynchronize();
		testKernel<<<1,1>>>();

		// Memory allocation on device
		CUDA_CHECK_RETURN(cudaMalloc(&points, 2 * hostNumberOfPoints * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&xOfNormals, 3 * hostNumberOfTriangles * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&yOfNormals, 3 * hostNumberOfTriangles * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&neighbors, 3* hostNumberOfTriangles * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&forbidden, 3* hostNumberOfTriangles * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&positionsOfNormalVectors, 3* hostNumberOfTriangles * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&triangleIndices, 3* hostNumberOfTriangles * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&cellTypes,hostNumberOfTriangles * hostNumberOfLevels * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&betaValues,hostNumberOfTriangles * hostNumberOfLevels * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&phi,hostNumberOfPoints * (hostNumberOfLevels) * sizeof(float)));

#if USE_IMPORTANCE==true /// This part only appears if we compile with the importance sampling
		CUDA_CHECK_RETURN(cudaMalloc(&importance,hostNumberOfPoints * (hostNumberOfLevels) * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&center_x,hostNumberOfTriangles * (hostNumberOfLevels) * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&center_y,hostNumberOfTriangles * (hostNumberOfLevels) * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&N_rays,hostNumberOfPoints * (hostNumberOfLevels) * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&surface,hostNumberOfPoints * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&devStates, iterations * threads * blocks * sizeof(curandState)));

#endif

		/// Copy data from host to device
		CUDA_CHECK_RETURN(cudaMemcpy(points, hostPoints, 2 * hostNumberOfPoints * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(xOfNormals, hostXOfNormals, 3 * hostNumberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(yOfNormals, hostYOfNormals, 3 * hostNumberOfTriangles * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(neighbors, hostNeighbors, 3* hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(forbidden,hostForbidden, 3* hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(positionsOfNormalVectors ,hostPositionsOfNormalVectors, 3* hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(triangleIndices ,hostTriangleIndices, 3* hostNumberOfTriangles * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(cellTypes,hostCellTypes, hostNumberOfTriangles *  hostNumberOfLevels * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(betaValues, hostBetaValues, hostNumberOfTriangles * (hostNumberOfLevels-1) * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(phi, hostPhi, hostNumberOfPoints * (hostNumberOfLevels) * sizeof(float), cudaMemcpyHostToDevice));

#if USE_IMPORTANCE==true /// This part only appears if we compile with the importance sampling
		CUDA_CHECK_RETURN(cudaMemcpy(importance, host_importance, hostNumberOfPoints * (hostNumberOfLevels) * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(center_x, host_center_x, hostNumberOfTriangles * hostNumberOfLevels * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(center_y, host_center_y, hostNumberOfTriangles * hostNumberOfLevels * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(N_rays, host_N_rays, hostNumberOfPoints * (hostNumberOfLevels) * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(surface, host_surface, hostNumberOfPoints * sizeof(int), cudaMemcpyHostToDevice));
#endif
	}


#if USE_IMPORTANCE==true
		fprintf(stderr, "\nStarting the Importance Sampling\n");
		random_setup_kernel <<< blocks, threads >>> ( devStates, time(NULL) );
		cudaThreadSynchronize();

		importanceKernel<<< blocks * iterations, threads>>>(devStates, points, xOfNormals, yOfNormals, positionsOfNormalVectors, neighbors, forbidden, cellTypes, hostNumberOfTriangles, betaValues, importance, N_rays, center_x, center_y, surface, blocks * threads * iterations * hostNumberOfPoints * (hostNumberOfLevels+1) );
#endif
	// start the Kernels
	{
		fprintf(stderr, "\nStarting the Naive Propagation\n");
		cudaEventRecord(start, 0);

		unsigned kernelcount=0;
		// start a new kernel for every(!) sample point of our mesh
		for(int point2D = 0; point2D < hostNumberOfPoints ; ++point2D){
			for(int level = 0; level < hostNumberOfLevels; ++level){
				cudaThreadSynchronize();
#if USE_IMPORTANCE==true
				raytraceStep<<< blocks, threads >>> ( devMTGPStates, phi, point2D, level, iterations, points, xOfNormals, yOfNormals, positionsOfNormalVectors, neighbors, forbidden, triangleIndices, cellTypes, betaValues, importance);
#else
				raytraceStep<<< blocks, threads >>> ( devMTGPStates, phi, point2D, level, iterations, points, xOfNormals, yOfNormals, positionsOfNormalVectors, neighbors, forbidden, triangleIndices, cellTypes, betaValues);

				if(kernelcount % 200 == 0)
					fprintf(stderr, "Sampling point %d done\n",kernelcount);
				kernelcount++;
#endif
			}
		}

		fprintf(stderr, "\nNaive Propagation done\n");
		cudaThreadSynchronize();
	}


	// Evaluate device data
	{
		fprintf(stderr, "\nStopping the Timer\n");
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&runtimeGpu, start, stop);
		fprintf(stderr, "\nRetrieving ASE Data\n");

		CUDA_CHECK_RETURN(cudaMemcpy(hostPhi, phi, hostNumberOfPoints * (hostNumberOfLevels) * sizeof(float), cudaMemcpyDeviceToHost));
#if USE_IMPORTANCE==true
		CUDA_CHECK_RETURN(cudaMemcpy(host_importance, importance, hostNumberOfPoints * (hostNumberOfLevels) * sizeof(double), cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(host_N_rays, N_rays, hostNumberOfPoints * (hostNumberOfLevels) * sizeof(int), cudaMemcpyDeviceToHost));
#endif
		fprintf(stderr, "done\n");
		for(int i=0; i< hostNumberOfPoints;++i){
			for(int j=0 ; j<hostNumberOfLevels ; ++j)
			{
				int pos = i*hostNumberOfLevels+j;
				hostPhi[pos] =float( (double(hostPhi[pos]) / (raysPerSample * 4.0f * 3.14159)));
				double gain_local = double(hostNTot)*(betaCellsVector->at(pos))*double(hostSigmaE+hostSigmaA)-double(hostNTot*hostSigmaA);
				ase->at(pos) = gain_local*hostPhi[pos] / double(hostCrystalFluorescence);
			}
		}
		fprintf(stderr, "\nValues in Vector\n");
	}
		fprintf(stderr, "\nFreeing Memory\n");

	// Free memory on device
	{
		cudaFree(points);
		cudaFree(xOfNormals);
		cudaFree(yOfNormals);
		cudaFree(neighbors);
		cudaFree(forbidden);
		cudaFree(positionsOfNormalVectors);
		cudaFree(betaValues);
#if USE_IMPORTANCE==true
		cudaFree(importance);
		cudaFree(N_rays);
		cudaFree(center_x);
		cudaFree(center_y);
#endif
	}

	// print statistics
	   {
	/*
	   fprintf(stderr, "\n\n");
	   fprintf(stderr, "Vertices       : %d\n", host_size_p * (host_mesh_z+1));
	   fprintf(stderr, "Levels         : %d\n", host_mesh_z);
	   fprintf(stderr, "Prisms         : %d\n", host_N_cells * host_mesh_z);
	   fprintf(stderr, "Rays per Vertex: %d\n", raysPerSample);
	   fprintf(stderr, "Rays Total     : %d\n", raysPerSample * host_size_p * (host_mesh_z+1));
	   fprintf(stderr, "GPU Blocks     : %d\n", blocks);
	   fprintf(stderr, "iterations     : %d\n", iterations);
	   fprintf(stderr, "GPU Threads    : %d\n", threads * blocks);
	   fprintf(stderr, "Runtime_GPU    : %f s\n", runtimeGpu / 1000.0);
	   fprintf(stderr, "\n");
	*/
	   }
	{
		cudaFreeHost(hostBetaValues);
		cudaFreeHost(hostXOfNormals);
		cudaFreeHost(hostYOfNormals);
		cudaFreeHost(hostCellTypes);
		cudaFreeHost(hostTriangleIndices);
		cudaFreeHost(hostForbidden);
		cudaFreeHost(hostNeighbors);
		cudaFreeHost(hostPositionsOfNormalVectors);
		cudaFreeHost(hostPoints);
	}
		fprintf(stderr, "done\n");
	cudaDeviceReset();
	return runtimeGpu;
}

