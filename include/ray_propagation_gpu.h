#ifndef RAY_PROPAGATION_GPU_KERNEL_H
#define RAY_PROPAGATION_GPU_KERNEL_H
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include <vector>
#include "curand_kernel.h"
#include "datatypes.h"
/* include MTGP host helper functions */
#include <curand_mtgp32_host.h>
/* include MTGP pre-computed parameter sets */
#include <curand_mtgp32dc_p_11213.h>
#include <cuda_runtime_api.h>
#include "testdata_transposed.h"

#define TEST_VALUES false 
#define USE_IMPORTANCE false
#define DIVIDE_PI true
#define SMALL 1E-06

#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x) != CURAND_STATUS_SUCCESS) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

//----------------------------------------------------
// Device Code
//----------------------------------------------------

__device__ double clad_abs;
__device__ double N_tot;
__device__ double sigma_e;
__device__ double sigma_a;
__device__ double z_mesh;
__device__ int mesh_z;
__device__ int clad_num;
__device__ int size_p;
__device__ int N_cells;




__device__ int selectTriangle(int id, int trianglesInOneLevel){
	int totalNumberOfThreads = blockDim.x * gridDim.x;
	int threadsPer2DTriangle = ceil( float(totalNumberOfThreads) / float(trianglesInOneLevel) );
	return id / threadsPer2DTriangle;
}

__device__ int selectLevel(int id, int totalNumberOfLevels){
	int totalNumberOfThreads = blockDim.x * gridDim.x;
	int threadsPerLevel = ceil( float(totalNumberOfThreads) / float(totalNumberOfLevels) );
	return id / threadsPerLevel;
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
__device__ float rayPropagationGpu(double x_pos, double y_pos, double z_pos, double x_dest, double y_dest, double z_dest, int t_start, int mesh_start,  double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int *forbidden, int* cell_type, double* beta_v){
	//    in first try no reflections
	//
	//    create the vector between both points and calculate which surface would be the shortest to reach.
	//    then get the length, make the integration, get the information about the next cell out of the array
	//    set the point to the surface (this surface is "forbidden" in the calculations)
	//    proceed until you hit a the point or the surface
	
	double vec_x, vec_y,vec_z;
	double distance, length, length_help;
	double nominator, denominator;
	double gain=1.;
	int tri, cell_z; // the current triangle number and position concerning the z's
	int tri_next, cell_z_next, forb, forb_dump;
	int offset;
#if TEST_VALUES==true
	double testDistance = 0;
#endif


	//    initial positions
	tri = t_start;
	cell_z = mesh_start;

	// direction-vector (without reflections)
	vec_x = (x_dest - x_pos);
	vec_y = (y_dest - y_pos);
	vec_z = (z_dest - z_pos);

	// total distance to travel
	distance_total = sqrt(vec_x*vec_x+vec_y*vec_y+vec_z*vec_z);

	// normalized direction-vector
	vec_x = vec_x/distance_total;
	vec_y = vec_y/distance_total;
	vec_z = vec_z/distance_total;

	// remaining distance to travel
	double distance_remaining = distance_total;

	// at the beginning, all surfaces are possible
	forb = -1;

	for(;;)
	{
		// the length of the ray-part inside the current prism. We try to minimize this value
		length = distance_remaining;

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
			denominator = n_x[tri]*vec_x + n_y[tri]*vec_y;
			// see if we intersect at all
			if (denominator != 0.0)
			{
				nominator = (n_x[tri]*p_in[n_p[tri]] + n_y[tri]*p_in[n_p[tri]+ size_p]) - (n_x[tri]*x_pos + n_y[tri]*y_pos);
				length_help = nominator/denominator;
				// if we found a new smallest length, use it
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
					forb_dump = (forbidden[tri]);	
					tri_next = neighbors[tri];
					cell_z_next = cell_z;

				}
			}
		}

		// see forb !=0 case
		if (forb != 1){
			//offset, since the 3 rectangular surfaces are stored at different positions in the array
			offset = tri+N_cells;
			denominator = n_x[offset]*vec_x + n_y[offset]*vec_y;
			if (denominator != 0.0)
			{
				nominator = (n_x[offset]*p_in[n_p[offset]] + n_y[offset]*p_in[n_p[offset]+ size_p]) - (n_x[offset]*x_pos + n_y[offset]*y_pos);
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
					forb_dump = (forbidden[offset]);
					tri_next = neighbors[offset];
					cell_z_next = cell_z;
				}
			}
		}

		// see forb !=0 case
		if (forb !=2){
			offset = tri+2*N_cells;
		denominator = n_x[offset]*vec_x + n_y[offset]*vec_y;
			if (denominator != 0.0)
			{
				nominator = (n_x[offset]*p_in[n_p[offset]] + n_y[offset]*p_in[n_p[offset]+ size_p]) - (n_x[offset]*x_pos + n_y[offset]*y_pos);
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
					forb_dump = (forbidden[offset]);
					tri_next = neighbors[offset];
					cell_z_next = cell_z;
				}
			}
		}

		// if-structure "optimized"
		denominator = z_pos*vec_z;
		if (denominator != 0.0){
			if (forb != 3){
				{
					nominator = (cell_z+1)* z_mesh - z_pos;
					length_help = nominator/denominator;
					if (length_help < length && length_help > 0.0)
					{
						length = length_help;
						//decider = 3;
						forb_dump = 4; // you are not allowed to go down in the next step
						tri_next = tri;
						cell_z_next = cell_z + 1;
					}
				}
			}

			// next is the lower plane
			if (forb != 4){
				nominator = (cell_z)* z_mesh - z_pos;
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
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
		// take the beta_v[tri+cell_z*N_cells] 
		//
		// at this position do the decision whether it is a gain part or cladding
		// it might be absorbing or amplifying, for the cladding only absorbing
		// a simple "if then"

		if (cell_type[tri] == clad_num){
			gain *= exp((-1)*(clad_abs * length));
		}
		else {
			gain *= exp(length* N_tot*(beta_v[tri+cell_z*N_cells]*(sigma_e + sigma_a)-sigma_a));
		}
#if TEST_VALUES==true
		testDistance += length;
#endif

		// the remaining distance is decreased by the length we travelled through the prism
		distance_remaining -= length;


		// if the distance between the destination and our current position is small enough, we are done
		if (fabs(distance_remaining)< SMALL)
			break;

		// now set the next cell and position
		x_pos = x_pos + length*vec_x;
		y_pos = y_pos + length*vec_y;
		z_pos = z_pos + length*vec_z;

		tri = tri_next;
		cell_z = cell_z_next;      

	}

#if TEST_VALUES==true
	assert(fabs(distance_total-testDistance) < 0.000001);
#endif

	return gain/(distance_total*distance_total);
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
    y_pos = p_in[point+size_p];
    z_pos = mesh_start * z_mesh;

    for (i_t=0;i_t<N_cells;i_t++)
    {

        for (i_z=0;i_z<(mesh_z-1);i_z++) //remember the definition differences MatLab/C for indices
        {
//            at this point replace the following routine with propagation(...)
//            later expand this with the beta/tau values...
			prop = naivePropagation(center_x[i_t], center_y[i_t], z_mesh*(i_z+0.5),  x_pos, y_pos, z_pos, i_t, i_z, p_in, n_x, n_y, n_p, neighbors, forbidden , cell_type, beta_v);
			// Propagation vom Zentrum jedes Prismas zu jedem Samplepunkt 
			//
//			prop = propagation(center_x[i_t], center_y[i_t], z_mesh*(i_z+0.5), x_pos, y_pos, z_pos, i_t, i_z);
            importance[i_t + i_z*N_cells] = beta_v[i_t+i_z*N_cells]*(prop);
            sum_phi += importance[i_t + i_z*N_cells];

        }
        surf_tot += surface[i_t];
    }

//    now calculate the number of rays
    for (i_t=0;i_t<N_cells;i_t++)
    {
        for (i_z=0;i_z<(mesh_z-1);i_z++) //remember the definition differences MatLab/C for indices
        {
//            this is the amount of the sampled rays out of the cells
            N_rays[i_t + i_z*N_cells] = (int)(floor(importance[i_t + i_z*N_cells]/sum_phi*NumRays));

            Rays_dump +=  N_rays[i_t + i_z*N_cells];
        }
    }

    rays_left = NumRays-Rays_dump;
//    distribute the remaining not distributed rays randomly
    if ((rays_left)>0)
    {
        for (i_r=0;i_r<rays_left;i_r++)
        {
            rand_t = (int )(curand_uniform(&localstate)*N_cells);
            rand_z = (int )(curand_uniform(&localstate)*(mesh_z-1));
            N_rays[rand_t + rand_z*N_cells]++;
        }
    }

//    now think about the mount of rays which would come out of this volume(surface)
//    dividing this number with the new amount of rays gives the final importance weight for this area!
    for (i_t=0;i_t<N_cells;i_t++)
    {
        for (i_z=0;i_z<(mesh_z-1);i_z++) //remember the definition differences MatLab/C for indices
        {
//            this is the amount of the sampled rays out of the cells
            if (N_rays[i_t + i_z*N_cells]>0)
            {
                importance[i_t + i_z*N_cells] = NumRays*surface[i_t]/surf_tot/N_rays[i_t + i_z*N_cells];
//                importance[i_t + i_z*N_cells] = NumRays*surface[i_t]/surf_tot;
            }
            else
            {
                importance[i_t + i_z*N_cells] = 0; // case of beta of this point == 0 e.g.
            }
        }
    }
}
#endif


/**
 * Initializes the global variables of the GPU with the correct values.
 * All those values are from the original propagation-function which we ported.
 */
__global__ void setupGlobalVariablesKernel ( double host_sigma_e, double host_sigma_a, int host_clad_num, int host_clad_abs, double host_N_tot, int host_N_cells, double host_z_mesh, int host_mesh_z, int host_size_p )
{
	sigma_e = host_sigma_e;	
	sigma_a = host_sigma_a;
	clad_num = host_clad_num;
	clad_abs = host_clad_abs;
	N_tot = host_N_tot;
	N_cells = host_N_cells;
	z_mesh = host_z_mesh;
	mesh_z = host_mesh_z;
	size_p = host_size_p;
	//printf("Sigma_e in setup=%f\tSigma_eHost=%f\n",sigma_e,host_sigma_e);
} 

#if USE_IMPORTANCE==true
__global__ void importanceKernel( curandState *globalState, double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int *forbidden, int* cell_type, int host_size_t, double* beta_v, double *importance, int *N_rays, double *center_x, double *center_y, int *surface, int NumRays) {
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	for(int i=0; i< host_size_t; ++i){
		for(int j=0; j< mesh_z; ++j){
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
	printf("\nSigma_e=%.30f",sigma_e);
	printf("\nSigma_a=%.30f",sigma_a);
	printf("\nmesh_z=%d",mesh_z);
	printf("\nz_mesh_=%.30f",z_mesh);
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
__global__ void raytraceStep( curandStateMtgp32* globalState, float* phi, int point2D, int level, int iterations, double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int *forbidden, int* t_in, int* cell_type, int host_size_t, double* beta_v, double *importance) {
#else
__global__ void raytraceStep( curandStateMtgp32* globalState, float* phi, int point2D, int level, int iterations, double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int *forbidden, int* t_in, int* cell_type, int host_size_t, double* beta_v) {
#endif
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	double gain = 0.;
	const int endPointX = p_in[point2D];
	const int endPointY = p_in[ size_p + point2D];
	const int endPointZ = level* z_mesh;


	// this should give the same start values multiple times (so that every thread uses the same prism, which yields
	// big benefits for the memory access (and caching!)
	const int startTriangle = selectTriangle(id,host_size_t); 
	const int startLevel = selectLevel(id, mesh_z); 

	// the indices of the vertices of the starttriangle
	const int t1 = t_in[startTriangle];
	const int t2 = t_in[startTriangle+ N_cells];
	const int t3 = t_in[startTriangle+2*N_cells];

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
		double zRand = (startLevel + curand_uniform(&globalState[blockIdx.x]))* z_mesh;
		double xRand = p_in[t1]*u + p_in[t2]*v + p_in[t3]*w;
		double yRand = p_in[ size_p + t1]*u + p_in[ size_p + t2]*v + p_in[ size_p + t3]*w;
		__syncthreads();
		gain += double(rayPropagationGpu(xRand, yRand, zRand, endPointX, endPointY, endPointZ, startTriangle, startLevel ,p_in, n_x, n_y, n_p, neighbors, forbidden , cell_type, beta_v));
	}

	// do the multiplication just at the end of all iterations
	// (gives better numeric behaviour)
	gain *=  beta_v[startTriangle + N_cells*startLevel];
#if USE_IMPORTANCE==true
	atomicAdd(&(phi[point2D + level*size_p]),float(gain*importance[startTriangle + N_cells*startLevel]));  
#else
	atomicAdd(&(phi[point2D + level*size_p]),float(gain)); 
#endif
}

//----------------------------------------------------
// Host Code
//----------------------------------------------------
float runRayPropagationGpu(std::vector<double> *ase, unsigned &threads, unsigned &blocks, unsigned &totalNumberOfRays)
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
	unsigned raysPerSample = ceil(totalNumberOfRays/float(host_size_p * (host_mesh_z+1)));
	threads = 256;
	blocks = 200;
	int iterations = ceil(float(raysPerSample) / float(blocks * threads));
	//fprintf(stderr, "raysPerSample=%d\n",raysPerSample);
	//fprintf(stderr, "interations=%d\n",iterations);
	raysPerSample = threads * blocks * iterations;
	totalNumberOfRays = unsigned(raysPerSample * (host_size_p * (host_mesh_z+1)));
	//fprintf(stderr, "After Normalization:\n");
	//fprintf(stderr, "raysPerSample=%d\n",raysPerSample);
	//fprintf(stderr, "totalNumberOfRays=%d\n",totalNumberOfRays);
#if TEST_VALUES==true
	assert(raysPerSample == threads * blocks * iterations);
#endif


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
	double  *p_in, *n_x, *n_y, *beta_v;
	float *phi;
	int *forbidden, *n_p, *neighbors, *t_in, *cell_type;


	/** Variables for benchmarking and results
	 * \var hostPhi will contain the ASE-Flux integral (not normalized)
	 */
	float runtimeGpu = 0.0;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float hostPhi[host_size_p * (host_mesh_z +1)];
	for(int i=0;i<host_size_p*(host_mesh_z+1);++i){
			hostPhi[i] = 0.;
	}


	/** The Mersenne Twister PRNG
	 * These calls are more or less directly from the curand
	 * \var devMTGPStates the current randomness states of each block (may not exceed 200!)
	 *
	 */
	curandStateMtgp32 *devMTGPStates;
	{
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
	double host_importance[host_size_p * (host_mesh_z+1)];
	int host_N_rays[host_size_t * host_mesh_z];
	for(int i=0;i<host_size_p*(host_mesh_z+1);++i){
		host_importance[i] = 0.;
	}
	for(int i=0;i<host_size_t*host_mesh_z;++i){
		host_N_rays[i] = 0;
	}
#endif



	/** Allocation of memory on the GPU and setting of global GPU-variables
	 *
	 */
	{
		//Create constant values on GPU
		setupGlobalVariablesKernel<<<1,1>>>(host_sigma_e, host_sigma_a, host_clad_num, host_clad_abs, host_N_tot, host_size_t, host_z_mesh, host_mesh_z, host_size_p); //@OPTIMIZE: initialize the constants as constants...

		cudaThreadSynchronize();

		// Memory allocation on device
		CUDA_CHECK_RETURN(cudaMalloc(&p_in, 2 * host_size_p * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&n_x, host_size_t * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&n_y, host_size_t * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&neighbors, 3* host_size_t * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&forbidden, 3* host_size_t * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&n_p, 3* host_size_t * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&t_in, 3* host_size_t * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&cell_type,host_size_t * host_mesh_z * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&beta_v,host_size_t * host_mesh_z * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&phi,host_size_p * (host_mesh_z +1) * sizeof(float)));

#if USE_IMPORTANCE==true /// This part only appears if we compile with the importance sampling
		CUDA_CHECK_RETURN(cudaMalloc(&importance,host_size_p * (host_mesh_z +1) * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&center_x,host_size_t * (host_mesh_z) * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&center_y,host_size_t * (host_mesh_z) * sizeof(double)));
		CUDA_CHECK_RETURN(cudaMalloc(&N_rays,host_size_p * (host_mesh_z +1) * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&surface,host_size_p * sizeof(int)));
		CUDA_CHECK_RETURN(cudaMalloc(&devStates, iterations * threads * blocks * sizeof(curandState)));

#endif

		/// Copy data from host to device
		CUDA_CHECK_RETURN(cudaMemcpy(p_in, host_p_in, 2 * host_size_p * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(n_x, host_n_x, host_size_t * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(n_y, host_n_y, host_size_t * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(neighbors, host_neighbors, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(forbidden,host_forbidden, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(n_p ,host_n_p, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(t_in ,host_t_in, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(cell_type,host_cell_type, host_size_t *  host_mesh_z * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(beta_v, host_beta_v, host_size_t * host_mesh_z * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(phi, hostPhi, host_size_p * (host_mesh_z+1) * sizeof(float), cudaMemcpyHostToDevice));

#if USE_IMPORTANCE==true /// This part only appears if we compile with the importance sampling
		CUDA_CHECK_RETURN(cudaMemcpy(importance, host_importance, host_size_p * (host_mesh_z+1) * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(center_x, host_center_x, host_size_t * host_mesh_z * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(center_y, host_center_y, host_size_t * host_mesh_z * sizeof(double), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(N_rays, host_N_rays, host_size_p * (host_mesh_z+1) * sizeof(int), cudaMemcpyHostToDevice));
		CUDA_CHECK_RETURN(cudaMemcpy(surface, host_surface, host_size_p * sizeof(int), cudaMemcpyHostToDevice));
#endif
	}


#if USE_IMPORTANCE==true
		fprintf(stderr, "\nStarting the Importance Sampling\n");
		random_setup_kernel <<< blocks, threads >>> ( devStates, time(NULL) );
		cudaThreadSynchronize();

		importanceKernel<<< blocks * iterations, threads>>>(devStates, p_in, n_x, n_y, n_p, neighbors, forbidden, cell_type, host_size_t, beta_v, importance, N_rays, center_x, center_y, surface, blocks * threads * iterations * host_size_p * (host_mesh_z+1) );
#endif
	// start the Kernels
	{
		fprintf(stderr, "\nStarting the Naive Propagation\n");
		cudaEventRecord(start, 0);

		for(int point2D = 0; point2D < host_size_p ; ++point2D){
			for(int level = 0; level <= host_mesh_z; ++ level){
				cudaThreadSynchronize();
#if USE_IMPORTANCE==true
				raytraceStep<<< blocks, threads >>> ( devMTGPStates, phi, point2D, level, iterations, p_in, n_x, n_y, n_p, neighbors, forbidden, t_in, cell_type, host_size_t, beta_v, importance);
#else
				raytraceStep<<< blocks, threads >>> ( devMTGPStates, phi, point2D, level, iterations, p_in, n_x, n_y, n_p, neighbors, forbidden, t_in, cell_type, host_size_t, beta_v);

#endif
			}
		}

		cudaThreadSynchronize();
	}

	// Evaluate device data
	{
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&runtimeGpu, start, stop);

		CUDA_CHECK_RETURN(cudaMemcpy(hostPhi, phi, host_size_p * (host_mesh_z+1) * sizeof(int), cudaMemcpyDeviceToHost));
#if USE_IMPORTANCE==true
		CUDA_CHECK_RETURN(cudaMemcpy(host_importance, importance, host_size_p * (host_mesh_z+1) * sizeof(double), cudaMemcpyDeviceToHost));
		CUDA_CHECK_RETURN(cudaMemcpy(host_N_rays, N_rays, host_size_p * (host_mesh_z+1) * sizeof(int), cudaMemcpyDeviceToHost));
#endif
		for(int i=0; i< host_size_p*(host_mesh_z+1); ++i){
			//fprintf(stderr, "\nPhi_ase[%d]= %.10f",i, hostPhi[i] / raysPerSample);
#if DIVIDE_PI==true
			ase->at(i) = (double(double(hostPhi[i]) / (raysPerSample * 4.0f * 3.14159)));
#else
			ase->at(i) = (double(double(hostPhi[i]) / raysPerSample));
#endif
		}
	}

	// Free memory on device
	{
		cudaFree(p_in);
		cudaFree(n_x);
		cudaFree(n_y);
		cudaFree(neighbors);
		cudaFree(forbidden);
		cudaFree(n_p);
		cudaFree(beta_v);
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
	cudaDeviceReset();
	return runtimeGpu;
}

#endif
