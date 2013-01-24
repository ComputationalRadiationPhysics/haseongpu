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
#include "testdata_1.h"

#define REAL_VALUES false
#define SMALL 1E-06
#define CUDA_CHECK_RETURN(value) {				\
	cudaError_t _mCudaStat = value;				\
	if (_mCudaStat != cudaSuccess) {				\
		fprintf(stderr, "Error %s at line %d in file %s\n",	\
				cudaGetErrorString(_mCudaStat), __LINE__, __FILE__);	\
		exit(1);							\
	}								\
}
#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if((x) != CURAND_STATUS_SUCCESS) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__); \
	return EXIT_FAILURE;}} while(0)

//----------------------------------------------------
// Device Code
//----------------------------------------------------

/**
  @brief Calculates A-B for 2 float4-based inputs
 **/
__device__ PointCu subtractPoints(PointCu A, PointCu B){
	PointCu C;
	C.x = A.x - B.x;
	C.y = A.y - B.y;
	C.z = A.z - B.z;
	C.w = A.w - B.w;
	return C;
}

__device__ double clad_abs;
__device__ double N_tot;
__device__ double sigma_e;
__device__ double sigma_a;
__device__ double z_mesh;
__device__ int mesh_z;
__device__ int clad_num;
__device__ int size_p;
__device__ int N_cells;


__device__ RayCu generateRayGpu(PointCu vertexPoint, PrismCu startPrism, curandStateMtgp32 *randomstate, int bid){
	float u = curand_uniform(&randomstate[bid]);
	float v = curand_uniform(&randomstate[bid]);
	if((u+v) > 1){ //OPTIMIZE: remove if
		u = 1-u;
		v = 1-v;
	}
	const float w = 1-(u+v);

	PointCu A = startPrism.t1.A;
	PointCu B = startPrism.t1.B;
	PointCu C = startPrism.t1.C;

	// Get x and y coordinates from the random barycentric values
	const float xRand = u*A.x + v*B.x + w*C.x ;
	const float yRand = u*A.y + v*B.y + w*C.y ;

	// Take one of the given z-coordinates and add a random part of the prism height
	const float zRand = A.z + curand_uniform(&randomstate[bid]) * startPrism.t1.A.w;

	float ase=0.f;

	// Take the values to assemble a ray
	RayCu r = {
		{xRand, yRand, zRand, ase},
		vertexPoint};
	return r;
}

__device__ float distance(PointCu a, PointCu b){
	float d = sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y),2) + pow((b.z - a.z),2));
	return fabs(d);
}

__device__ PrismCu selectPrism(int id, PrismCu prisms[], int totalNumberOfPrisms){
	int totalNumberOfThreads = blockDim.x * gridDim.x;
	int threadsPerPrism = ceil( float(totalNumberOfThreads) / float(totalNumberOfPrisms) );
	int prism = id / threadsPerPrism;

	return prisms[prism];
}

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

__device__ float propagate(RayCu ray, PrismCu prisms[], PrismCu startprism){
	//	float gain = 1.f;
	//	float vecX = ray.direction.x - ray.P.x;
	//	float vecY = ray.direction.y - ray.P.y;
	//	float vecZ = ray.direction.z - ray.P.z;
	//
	//	float distanceTotal = sqrt(vecX*vecX+vecY*vecY+vecZ*vecZ);
	//	float distance = distanceTotal;
	//	float length = distanceTotal;
	//	vecX /= distanceTotal;
	//	vecY /= distanceTotal;
	//	vecZ /= distanceTotal;
	//
	//	PrismCu current = startprism;
	//
	//	int i=0;
	//	float lengthHelp = 0.f;
	//	float d = 0.f;
	//	float nominator = 0.f;
	//	float denominator = 0.f;
	//	TriangleCu t1;
	//	VectorCu AB,AC;
	//	PlaneCu pl;
	//
	//	for(;;){
	//		length = distance;
	//		//generate the triangle surfaces of the prism
	//		t1 = current.t1;
	//		TriangleCu t2 = { 
	//			{t1.A.x, t1.A.y, t1.A.z + t1.A.w, 1},
	//			{t1.B.x, t1.B.y, t1.B.z + t1.B.w, 1},
	//			{t1.C.x, t1.C.y, t1.C.z + t1.C.w, 1}
	//		};
	//
	//		// OPTIMIZE: make use of the rectangles!
	//		TriangleCu surfaces[8] = {
	//			t1,
	//			t2,
	//			{t1.A, t1.B, t2.A},
	//			{t1.B, t2.B, t2.A},
	//			{t1.B, t1.C, t2.C},
	//			{t1.B, t2.B, t2.C},
	//			{t1.A, t1.C, t2.C},
	//			{t1.A, t2.A, t2.C}
	//		};
	//
	//		for(i=0; i<8 ; ++i){ //OPTIMIZE: unroll, so that every surface can be optimized differently
	//			// get the generating vectors for the plane
	//			AB = subtractPoints(surfaces[i].B, surfaces[i].A);
	//			AC = subtractPoints(surfaces[i].C, surfaces[i].A);
	//
	//			pl.P = surfaces[i].A;
	//			// cross product of the vectors
	//			pl.normal.x = AB.y*AC.z - AB.z*AC.y;
	//			pl.normal.y = AB.z*AC.x - AB.x*AC.z;
	//			pl.normal.z = AB.x*AC.y - AB.y*AC.x;
	//
	//			// direction * pl.normal
	//			denominator = (ray.direction.x * pl.normal.x) + (ray.direction.y * pl.normal.y) + (ray.direction.z * pl.normal.z);
	//			if(denominator != 0.f) //OPTIMIZE: check if we have a lot of branch diversion, or if all threads behave the same
	//			{
	//				// A * pl.normal
	//				d = (surfaces[i].A.x * pl.normal.x) + (surfaces[i].A.y * pl.normal.y) + (surfaces[i].A.z * pl.normal.z);
	//				// d - (P * pl.normal)
	//				nominator = d - ((ray.P.x * pl.normal.x) + (ray.P.y * pl.normal.y) + (ray.P.z * pl.normal.y)); 
	//				lengthHelp = nominator/denominator;
	//				if(lengthHelp < length && lengthHelp > 0.f) //OPTIMIZE: most threads should do the same?
	//				{
	//					length = lengthHelp;
	//				}
	//			}
	//		}
	//
	//
	//		//with the new length, get the gain and add it
	//		// @TODO
	//		gain += length;
	//
	//		// calculate values for next iteration
	//		distance -= length;
	//
	//		
	//		if(abs(distance) < SMALL)
	//		{
	//			break;
	//		}
	//
	//		ray.P.x += length*vecX;
	//		ray.P.y += length*vecY;
	//		ray.P.z += length*vecZ;
	//
	//		//@TODO:
	//		// calculate the next PRISM (maybe with help of some neighbor-datastructure?
	//
	//	}
	//
	//
	//	return gain;
	return 0;
}

__device__ float naive_propagation(double x_pos, double y_pos, double z_pos, double x_dest, double y_dest, double z_dest, int t_start, int mesh_start,  double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int *forbidden, int* cell_type, double* beta_v){
	//    in first try no reflections
	//    calculate the vector and make the the calculation, which surface would be the shortest to reach
	//    then get the length, make the integration, get the information about the next cell out of the array
	//    set the point to the surface (this surface is "forbidden" in the calculations)
	//    proceed until you hit a the point or the surface
	//    if you are closer then "small" stop and return the value
	double vec_x, vec_y,vec_z, norm;
	double distance, length, length_help, distance_total;
	double nominator, denominator;
	double gain=1;
	int tri, cell_z; // the current triangle number and position concerning the z's
	int decider; // which one is the shortest - info
	int tri_next, cell_z_next, forb, forb_dump;
	int ct;

#if REAL_VALUES==true
	gain=1;
#else
	gain=0;
#endif


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

		//        mexPrintf("forb: %i\n",forb);
		//        mexEvalString("drawnow;");

		//        try the triangle faces
		//        remember the correlation between the normals and the points
		//        n1: p1-2, n2: p1-3, n3:p2-3
		//        the third coordinate (z) of the particpating points for the surfaces can be set to be z=0, 
		//        as everything uses triangular "tubes/prisms", as well as n_z=0 in this case!
		if (forb != 0){
			denominator = n_x[tri]*vec_x + n_y[tri]*vec_y;
			if (denominator != 0.0)
			{
				nominator = (n_x[tri]*p_in[n_p[tri]] + n_y[tri]*p_in[n_p[tri]+ size_p]) - (n_x[tri]*x_pos + n_y[tri]*y_pos);
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
					decider = 0;
					forb_dump = (forbidden[tri]);

				}
			}
		}

		if (forb != 1){
			denominator = n_x[tri+ N_cells]*vec_x + n_y[tri+ N_cells]*vec_y;
			if (denominator != 0.0)
			{
				nominator = (n_x[tri+ N_cells]*p_in[n_p[tri+ N_cells]] + n_y[tri+ N_cells]*p_in[n_p[tri+ N_cells]+ size_p]) - (n_x[tri+ N_cells]*x_pos + n_y[tri+ N_cells]*y_pos);
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
					decider = 1;
					forb_dump = (forbidden[tri+ N_cells]);
				}
			}
		}

		if (forb !=2){
			denominator = n_x[tri+2* N_cells]*vec_x + n_y[tri+2* N_cells]*vec_y;
			if (denominator != 0.0)
			{
				nominator = (n_x[tri+2* N_cells]*p_in[n_p[tri+2*  N_cells]] + n_y[tri+2* N_cells]*p_in[n_p[tri+2* N_cells]+ size_p]) - (n_x[tri+2* N_cells]*x_pos + n_y[tri+2* N_cells]*y_pos);
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
					decider = 2;
					forb_dump = (forbidden[tri+2* N_cells]);
				}
			}
		}

		//        try the horizontal planes, which one is the shortest, n_x and n_y are zero!, n_z =1!
		//        at first the upper plane
		if (forb != 3){
			denominator = z_pos*vec_z;
			if (denominator != 0.0)
			{
				nominator = (cell_z+1)* z_mesh - z_pos;
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
					decider = 3;
					forb_dump = 4; // you are not allowed to go down in the next step
				}
			}
		}

		//        next is the lower plane
		if (forb != 4){
			denominator = z_pos*vec_z;

			if (denominator != 0.0)
			{
				nominator = (cell_z)* z_mesh - z_pos;
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
					decider = 4;
					forb_dump = 3; // you are not allowed to go up in the next step
				}
			}
		}

		forb = forb_dump;


		//        now make a switch to differ the different cases
		//@TODO: include this into the if statements?
		switch(decider){

			case 0:
				//                this is the case for the intersection with the first choice triangle-surface
				tri_next = neighbors[tri];
				cell_z_next = cell_z;
				break;

			case 1:
				//                second triangle surface
				tri_next = neighbors[tri+ N_cells];
				cell_z_next = cell_z;
				break;

			case 2:
				//                third triangle surface
				tri_next = neighbors[tri+2* N_cells];
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

#if REAL_VALUES==true

		ct = cell_type[tri]; 
		if (ct == clad_num){
			gain = gain * exp((-1)*(clad_abs * length));
		}
		else {
			gain = gain * exp(N_tot*(beta_v[tri+cell_z*N_cells]*(sigma_e + sigma_a)-sigma_a)*length);
		}
#else

		gain += length; 

#endif
		distance -= length;


		if (fabs(distance)< SMALL)
			break;

		//        now set the next cell
		x_pos = x_pos + length*vec_x;
		y_pos = y_pos + length*vec_y;
		z_pos = z_pos + length*vec_z;

		tri = tri_next;
		cell_z = cell_z_next;      

	}

#if REAL_VALUES==true
	gain /= (distance_total*distance_total); //@OPTIMIZE
#endif
	return gain;
}


__global__ void setupKernel ( double host_sigma_e, double host_sigma_a, int host_clad_num, int host_clad_abs, double host_N_tot, int host_N_cells, double host_z_mesh, int host_mesh_z, int host_size_p ){
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
	//int id = threadIdx.x + blockDim.x*blockIdx.x;
	//	curand_init ( seed, id, 0, &state[id] );
	//	// OPTIMIZE: Use MersenneTwister or even a better PRNG
} 
__global__ void testKernel (  ){
	if(threadIdx.x == 0)
		printf("\nSigma_e=%.30f",sigma_e);
	printf("\nSigma_a=%.30f",sigma_a);
	printf("\nmesh_z=%d",mesh_z);
	printf("\nz_mesh_=%.30f",z_mesh);
} 

// does the raytracing for a single ray (randomly generated) and a single (given) Vertex
__global__ void raytraceStep( curandStateMtgp32* globalState, float* phi, int point2D, int level, int iterations, double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int *forbidden, int* t_in, int* cell_type, int host_size_t, double* beta_v) {
	int id = threadIdx.x + blockIdx.x * blockDim.x;

	//OPTIMIZE: the Octree should/could produce a subset of the prism-array!

	int endpoint_x = p_in[point2D];
	int endpoint_y = p_in[ size_p + point2D];
	int endpoint_z = level* z_mesh;
	PointCu endpoint = {endpoint_x, endpoint_y, endpoint_z};
	double gain = 0.;
	float initial_distance = 50.;

	for (int i=0; i<iterations ; ++i){
		// this should give the same prism multiple times (so that every thread uses the same prism, which yields
		// big benefits for the memory access (and caching!)
		//PrismCu startprism = selectPrism(id, prisms, prismCount);	
		int starttriangle = selectTriangle(id,host_size_t); //@TODO: second parameter is number of 2D Triangles
		int startlevel = selectLevel(id, mesh_z); //@TODO: second parameter is number of Levels 

		// the indices of the vertices of the starttriangle
		int t_1 = t_in[starttriangle];
		int t_2 = t_in[starttriangle+ N_cells];
		int t_3 = t_in[starttriangle+2* N_cells];

		// random startpoint generation
		float  u = curand_uniform(&globalState[blockIdx.x]);
		float  v = curand_uniform(&globalState[blockIdx.x]);

		if((u+v)>1)
		{
			u = 1-u;
			v = 1-v;
		}

		float w = 1-u-v;

		// convert the random startpoint into coordinates
		double z_rand = (startlevel + curand_uniform(&globalState[blockIdx.x]))* z_mesh;
		double x_rand = p_in[t_1]*u + p_in[t_2]*v + p_in[t_3]*w;
		double y_rand = p_in[ size_p + t_1]*u + p_in[ size_p + t_2]*v + p_in[ size_p + t_3]*w;
		PointCu startpoint = {x_rand, y_rand, z_rand};

		//RayCu ray = generateRayGpu(vertices[vertex_index].P,startprism, globalState,blockIdx.x);
		initial_distance = distance(startpoint, endpoint);

		gain = naive_propagation(x_rand, y_rand, z_rand, endpoint_x, endpoint_y, endpoint_z, starttriangle, startlevel ,p_in, n_x, n_y, n_p, neighbors, forbidden , cell_type, beta_v); //@TODO: why is there a beta_v reference in the sequential-code?

		atomicAdd(&(phi[point2D + level*size_p]),float(gain)); //@TODO: importance-value 

#if REAL_VALUES!=true
		assert(fabs(gain-initial_distance) < 0.0001);
#endif
	}
}

//----------------------------------------------------
// Host Code
//----------------------------------------------------
int main(){

	// Variables from the mexFunction 
	double  *p_in, *n_x, *n_y, *beta_v;
	float *phi;
	int *forbidden, *n_p, *neighbors, *t_in, *cell_type;

	//int size_t=10, N_cells, size_p=10, size_z=2;

	int host_N_cells = host_size_t;

	//double z_mesh = 0.0356;
	//double  *host_p_in, *host_n_x, *host_n_y;
	//int *host_forbidden, *host_n_p, *host_neighbors, *host_t_in;
	// double host_clad_abs;

	/*
	   host_p_in = (double *)mxGetData(prhs[0]); //point coordinates in 2D . at first size_p x-values, then size_p y-values
	   host_n_x = (double *)mxGetData(prhs[4]); //normals of the facets, x-components // probably size_t values?? //@TODO check this
	   host_n_y = (double *)mxGetData(prhs[5]); //normals of the facets, y-components
	   host_forbidden = (int *)mxGetData(prhs[11]);//are the correspondance of the face used in the previous triangle, which must not be tested (rounding errors)
	   host_t_in = (int *)mxGetData(prhs[1]);  //association triangle-points - c-indexing-sytle!
	   host_n_p = (int *)mxGetData(prhs[10]); // gives the Index to one of the points in p_in which is in the plane of the normals - c-indexing-sytle!
	   host_neighbors = (int *)mxGetData(prhs[6]); //for each cell in t_in, its neighboring cell in plane geometry  - c-indexing-sytle!
	   host_size_t = (int )mxGetM(prhs[1]); //number of triangles per sheet
	   N_cells = host_size_t;
	   size_p = (int )mxGetM(prhs[0]); //number of points
	   z_mesh = (double)((double *)mxGetData(prhs[14]))[0];
	   size_z = (int )mxGetN(prhs[2]); //number of meshing in z-direction
	   host_clad_abs = (double)(cladabs[0]);

	 */
	//Variable definitions
	int threads = 256;
	float runtimeGpu = 0.0;
	cudaEvent_t start, stop;
	bool useGpu = true;
	curandStateMtgp32 *devMTGPStates;
	mtgp32_kernel_params *devKernelParams;
	float host_phi[host_size_p * (host_mesh_z +1)];

	// Generate testdata
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// GPU Raytracing
	unsigned rays_per_sample = 100; //10^5
	int rays_per_thread = ceil(rays_per_sample / float(threads));
	int blocks = 200;
	int iterations = ceil(rays_per_thread / float(blocks));

	if(useGpu){

		{

			//Create constant values on GPU
			setupKernel<<<1,1>>>(host_sigma_e, host_sigma_a, host_clad_num, host_clad_abs, host_N_tot, host_N_cells, host_z_mesh, host_mesh_z, host_size_p);


			cudaThreadSynchronize();
			//testKernel<<<1,1>>>();
			/*
			   cudaMemcpyToSymbol("sigma_a", &host_sigma_a, sizeof(double), cudaMemcpyHostToDevice);
			   cudaMemcpyToSymbol("clad_abs", &host_clad_abs, sizeof(double), cudaMemcpyHostToDevice);
			   cudaMemcpyToSymbol("clad_num", &host_clad_num, sizeof(int), cudaMemcpyHostToDevice);
			   cudaMemcpyToSymbol("N_tot", &host_N_tot, sizeof(double), cudaMemcpyHostToDevice);
			   cudaMemcpyToSymbol("N_cells", &host_N_cells, sizeof(int), cudaMemcpyHostToDevice);
			   cudaMemcpyToSymbol("z_mesh", &host_z_mesh, sizeof(double), cudaMemcpyHostToDevice);
			   cudaMemcpyToSymbol("mesh_z", &host_mesh_z, sizeof(int), cudaMemcpyHostToDevice);
			   cudaMemcpyToSymbol("size_p", &host_size_p, sizeof(int), cudaMemcpyHostToDevice);

			 */

			// Memory allocation on device (mexFunction Variables)
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

			// Copy data from host to device (mex Function Variables)
			CUDA_CHECK_RETURN(cudaMemcpy(p_in, host_p_in, 2 * host_size_p * sizeof(double), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(n_x, host_n_x, host_size_t * sizeof(double), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(n_y, host_n_y, host_size_t * sizeof(double), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(neighbors, host_neighbors, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(forbidden,host_forbidden, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(n_p ,host_n_p, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(t_in ,host_t_in, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(cell_type,host_cell_type, host_size_t *  host_mesh_z * sizeof(int), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(beta_v, host_beta_v, host_size_t * host_mesh_z * sizeof(double), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(phi, host_phi, host_size_p * (host_mesh_z+1) * sizeof(float), cudaMemcpyHostToDevice));
		}

		// Generating Random Numbers
		// Allocate space for prng states on device 
		CUDA_CALL(cudaMalloc((void **)&devMTGPStates, blocks * sizeof(curandStateMtgp32)));

		// Setup MTGP prng states 

		// Allocate space for MTGP kernel parameters 
		CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));

		// Reformat from predefined parameter sets to kernel format, 
		// and copy kernel parameters to device memory               
		CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));

		// Initialize one state per thread block
		CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, blocks, 1234));

		// State setup is complete 


		cudaEventRecord(start, 0);
		// start the Kernels

		for(int point2D = 0; point2D < host_size_p ; ++point2D){
			for(int level = 0; level <= host_mesh_z; ++ level){
				cudaThreadSynchronize();
				raytraceStep<<< blocks, threads >>> ( devMTGPStates, phi, point2D, level, iterations, p_in, n_x, n_y, n_p, neighbors, forbidden, t_in, cell_type, host_size_t, beta_v);
			}
		}

		cudaThreadSynchronize();
		CUDA_CHECK_RETURN(cudaMemcpy(host_phi, phi, host_size_p * (host_mesh_z+1) * sizeof(int), cudaMemcpyDeviceToHost));
		
		for(int i=0; i< host_size_p*(host_mesh_z+1); ++i){
			printf("\nPhi_ase= %.20f", host_phi[i] / rays_per_sample);
		
		}


		fprintf(stderr, "\n");

		// Evaluate device data
		{
			cudaEventRecord(stop, 0);
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&runtimeGpu, start, stop);
		}
		// Free memory on device
	}

	// print statistics
	//{

	fprintf(stderr, "\n");
	fprintf(stderr, "Vertices       : %d\n", host_size_p * (host_mesh_z+1));
	fprintf(stderr, "Levels         : %d\n", host_mesh_z);
	fprintf(stderr, "Prisms         : %d\n", host_N_cells * host_mesh_z);
	fprintf(stderr, "Rays per Vertex: %d\n", rays_per_sample);
	fprintf(stderr, "Rays Total     : %d\n", rays_per_sample * host_size_p * (host_mesh_z+1));
	fprintf(stderr, "GPU Blocks     : %d\n", blocks);
	fprintf(stderr, "iterations     : %d\n", iterations);
	fprintf(stderr, "GPU Threads    : %d\n", threads*blocks);
	fprintf(stderr, "Runtime_GPU    : %f s\n", runtimeGpu / 1000.0);
	fprintf(stderr, "\n");
	//	}
	// Cleanup;
	{
		// Cleanup mex Variables
		cudaFree(p_in);
		cudaFree(n_x);
		cudaFree(n_y);
		cudaFree(neighbors);
		cudaFree(forbidden);
		cudaFree(n_p);
		cudaFree(beta_v);

	}

	cudaDeviceReset();
	return 0;
}


