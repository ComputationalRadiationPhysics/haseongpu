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
// Structures
//----------------------------------------------------
typedef struct point {
	float x;
	float y;
	float z;
} POINT;

typedef struct vector {
	float x;
	float y;
	float z;
} VECTOR;

typedef struct ray {
	point start;
	vector direction;
} RAY;

typedef struct triangle {
	point a;
	point b;
	point c;
} TRIANGLE;

typedef struct plane {
	point start;
	vector normal;

} PLANE;

//------------------------------------------

//----------------------------------------------------
// Auxillary function declaration
//----------------------------------------------------

float distance(point a, point b);
void  printPoint(point p);

// New functions
bool  collide(TriangleCu t, PointCu p);
bool  collide(TriangleCu t, RayCu r);
bool  collide(PrismCu pr, RayCu r);
float4 toBarycentric(TriangleCu t, PointCu p);
PointCu intersection(PlaneCu p, RayCu r);
std::vector<TriangleCu> generateTriangles(int height, int width, float level);
std::vector<PrismCu> generatePrisms(int height, int width, float level);
std::vector<RayCu> generateRays(int height, int width, int level, unsigned maxRays);
RayCu   generateRay(int height, int weight, int level);
std::vector<VertexCu> generateSamples(int height, int width, int level);

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

__device__ float naive_propagation(double x_pos, double y_pos, double z_pos, double x_dest, double y_dest, double z_dest, int t_start, int mesh_start,  double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int N_cells, int size_p, int *forbidden, int z_mesh){
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
	//	int ct; // used to read out cell_type 


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
		//		ct = cell_type[tri]; //@TODO reimplement

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
				nominator = (n_x[tri]*p_in[n_p[tri]] + n_y[tri]*p_in[n_p[tri]+size_p]) - (n_x[tri]*x_pos + n_y[tri]*y_pos);
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
			denominator = n_x[tri+N_cells]*vec_x + n_y[tri+N_cells]*vec_y;
			if (denominator != 0.0)
			{
				nominator = (n_x[tri+N_cells]*p_in[n_p[tri+N_cells]] + n_y[tri+N_cells]*p_in[n_p[tri+N_cells]+size_p]) - (n_x[tri+N_cells]*x_pos + n_y[tri+N_cells]*y_pos);
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
				{
					length = length_help;
					decider = 1;
					forb_dump = (forbidden[tri+N_cells]);
				}
			}
		}

		if (forb !=2){
			denominator = n_x[tri+2*N_cells]*vec_x + n_y[tri+2*N_cells]*vec_y;
			if (denominator != 0.0)
			{
				nominator = (n_x[tri+2*N_cells]*p_in[n_p[tri+2*N_cells]] + n_y[tri+2*N_cells]*p_in[n_p[tri+2*N_cells]+size_p]) - (n_x[tri+2*N_cells]*x_pos + n_y[tri+2*N_cells]*y_pos);
				length_help = nominator/denominator;
				if (length_help < length && length_help > 0.0)
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
			denominator = z_pos*vec_z;
			if (denominator != 0.0)
			{
				nominator = (cell_z+1)*z_mesh - z_pos;
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
				nominator = (cell_z)*z_mesh - z_pos;
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

		//		if (ct == clad_num){
		//			gain = gain * exp(-clad_abs * length);
		//		}
		//		else {
		//			gain = gain * exp(N_tot*(beta_v[tri+cell_z*N_cells]*(sigma_e + sigma_a)-sigma_a)*length);
		//		}

		gain += length; //@TODO replace with actual gain computation
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

	gain /= (distance_total*distance_total);

	return gain;
}

__global__ void setupKernel ( curandState * state, unsigned long seed ){
	int id = threadIdx.x + blockDim.x*blockIdx.x;
	curand_init ( seed, id, 0, &state[id] );
	// OPTIMIZE: Use MersenneTwister or even a better PRNG
} 

// does the raytracing for a single ray (randomly generated) and a single (given) Vertex
__global__ void raytraceStep( curandStateMtgp32* globalState, VertexCu* vertices, int vertex_index, PrismCu* prisms, int prismCount, double *p_in, double *n_x, double *n_y, int *n_p, int *neighbors, int N_cells, int size_p, int host_size_t, int size_z, int *forbidden, double z_mesh, int* t_in) {
	int id = threadIdx.x + blockIdx.x * blockDim.x;

	//OPTIMIZE: the Octree should/could produce a subset of the prism-array!


	// this should give the same prism multiple times (so that every thread uses the same prism, which yields
	// big benefits for the memory access (and caching!)
	//PrismCu startprism = selectPrism(id, prisms, prismCount);	
	int starttriangle = selectTriangle(id,host_size_t); //@TODO: second parameter is number of 2D Triangles
	int startlevel = selectLevel(id,size_z); //@TODO: second parameter is number of Levels 

	// the indices of the vertices of the starttriangle
	int t_1 = t_in[starttriangle];
	int t_2 = t_in[starttriangle+N_cells];
	int t_3 = t_in[starttriangle+2*N_cells];

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
	double z_rand = (startlevel + curand_uniform(&globalState[blockIdx.x]))*z_mesh;
	double x_rand = p_in[t_1]*u + p_in[t_2]*v + p_in[t_3]*w;
	double y_rand = p_in[size_p + t_1]*u + p_in[size_p + t_2]*v + p_in[size_p + t_3]*w;

	//RayCu ray = generateRayGpu(vertices[vertex_index].P,startprism, globalState,blockIdx.x);
	//float initial_distance = distance(ray.P, ray.direction);

	float gain = 0.;
	//float gain = naive_propagation(x_rand, y_rand, z_rand, ray.direction.x, ray.direction.y, ray.direction.z, starttriangle, startlevel ,p_in, n_x, n_y, n_p, neighbors, N_cells, size_p, forbidden, z_mesh);
	//float gain = naive_propagation(ray.P.x, ray.P.y, ray.P.z, ray.direction.x, ray.direction.y, ray.direction.z, startprism.t1, startprism.t1.A.w ,p_in, n_x, n_y, n_p, neighbors, N_cells, size_p, forbidden, z_mesh);

	//printf("Thread: %d\t G=%.5f\t real_distance=%.5f\n",id,gain,initial_distance);

	printf("Thread: %d\t RAND=%f\n",id,curand_uniform(&globalState[0]));

	//@TODO: improve gain calculation (beta_v value, ImportanceSampling)
	//assert(fabs(gain-initial_distance) < 0.001);


	atomicAdd(&(vertices[vertex_index].P.w),gain);
}

//----------------------------------------------------
// Host Code
//----------------------------------------------------
int main(){

	// Variables from the mexFunction 
	double  *p_in, *n_x, *n_y;
	int *forbidden, *n_p, *neighbors, *t_in;
	//int size_t=10, N_cells, size_p=10, size_z=2;
	int size_z = mesh_z;

	int N_cells = host_size_t;
	//double z_mesh = 0.0356;
	//double  *host_p_in, *host_n_x, *host_n_y;
	//int *host_forbidden, *host_n_p, *host_neighbors, *host_t_in;
	
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
	 */
	//Variable definitions
	int threads = 16;
	unsigned prism_i, vertex_i;
	float runtimeGpu = 0.0;
	float runtimeCpu = 0.0;
	cudaEvent_t start, stop;
	bool useGpu = true;
	curandStateMtgp32 *devMTGPStates;
	mtgp32_kernel_params *devKernelParams;

	// Generate testdata
	std::vector<VertexCu> vertices = generateSamples(2, 2, 2);
	std::vector<PrismCu> prisms = generatePrisms(2, 2, 2);
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// GPU Raytracing
	PrismCu* hPrisms, *dPrisms;
	VertexCu* hVertices, *dVertices;
	int blocks = ceil(prisms.size() / float(threads));
	if(useGpu){

		//initialize memory
		{
			// Memory allocation on host
			CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&hPrisms, prisms.size() * sizeof(PrismCu), cudaHostAllocDefault));
			CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&hVertices, vertices.size() * sizeof(VertexCu), cudaHostAllocDefault));

			// Memory initialisation on host
			for(prism_i = 0; prism_i < prisms.size(); ++prism_i){
				hPrisms[prism_i] = prisms[prism_i];
			}
			for(prism_i = 0; prism_i < vertices.size() ; ++prism_i){
				hVertices[prism_i] = vertices[prism_i];
			}


			// Memory allocation on device
			CUDA_CHECK_RETURN(cudaMalloc(&dPrisms, prisms.size() * sizeof(PrismCu)));
			CUDA_CHECK_RETURN(cudaMalloc(&dVertices, vertices.size() * sizeof(PrismCu)));


			// Memory allocation on device (mexFunction Variables)
			CUDA_CHECK_RETURN(cudaMalloc(&p_in, 2 * size_p * sizeof(double)));
			CUDA_CHECK_RETURN(cudaMalloc(&n_x, host_size_t * sizeof(double)));
			CUDA_CHECK_RETURN(cudaMalloc(&n_y, host_size_t * sizeof(double)));
			CUDA_CHECK_RETURN(cudaMalloc(&neighbors, 3* host_size_t * sizeof(int)));
			CUDA_CHECK_RETURN(cudaMalloc(&forbidden, 3* host_size_t * sizeof(int)));
			CUDA_CHECK_RETURN(cudaMalloc(&n_p, 3* host_size_t * sizeof(int)));
			CUDA_CHECK_RETURN(cudaMalloc(&t_in, 3* host_size_t * sizeof(int)));



			// Copy data from host to device
			CUDA_CHECK_RETURN(cudaMemcpy(dPrisms, hPrisms, prisms.size() * sizeof(PrismCu), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(dVertices, hVertices, vertices.size() * sizeof(VertexCu), cudaMemcpyHostToDevice));


			// Copy data from host to device (mex Function Variables)
			CUDA_CHECK_RETURN(cudaMemcpy(p_in, host_p_in, 2 * size_p * sizeof(double), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(n_x, host_n_x, host_size_t * sizeof(double), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(n_y, host_n_y, host_size_t * sizeof(double), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(neighbors, host_neighbors, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(forbidden,host_forbidden, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(n_p ,host_n_p, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(t_in ,host_t_in, 3* host_size_t * sizeof(int), cudaMemcpyHostToDevice));
		}

		// Generating Random Numbers
		//setupKernel<<< blocks, threads >>> ( devStates, time(NULL) );
		/* Allocate space for prng states on device */
		CUDA_CALL(cudaMalloc((void **)&devMTGPStates, blocks * sizeof(curandStateMtgp32)));

		/* Setup MTGP prng states */

		/* Allocate space for MTGP kernel parameters */
		CUDA_CALL(cudaMalloc((void**)&devKernelParams, sizeof(mtgp32_kernel_params)));

		/* Reformat from predefined parameter sets to kernel format, */
		/* and copy kernel parameters to device memory               */
		CURAND_CALL(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213, devKernelParams));

		/* Initialize one state per thread block */
		CURAND_CALL(curandMakeMTGP32KernelState(devMTGPStates, mtgp32dc_params_fast_11213, devKernelParams, blocks, 1234));

		/* State setup is complete */


		fprintf(stderr, "\nbetween the kernel");
		cudaEventRecord(start, 0);
		// start the Kernels
		for(vertex_i = 0; vertex_i < vertices.size(); ++vertex_i){
			raytraceStep<<< blocks, threads >>> ( devMTGPStates, dVertices, vertex_i,  dPrisms, prisms.size(),p_in, n_x, n_y, n_p, neighbors, N_cells, size_p, host_size_t, size_z, forbidden, z_mesh, t_in);
		}

		fprintf(stderr, "\nafter the kernel");

		CUDA_CHECK_RETURN(cudaMemcpy(hVertices, dVertices, vertices.size() * sizeof(VertexCu), cudaMemcpyDeviceToHost));

		fprintf(stderr, "\n");
		for(vertex_i = 0; vertex_i < vertices.size(); ++vertex_i){
			fprintf(stderr, "Vertex %d:\t G=%.5f\n", vertex_i, hVertices[vertex_i].P.w);
		}


		// Evaluate device data
		{
			cudaEventRecord(stop, 0);
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&runtimeGpu, start, stop);
		}
		// Free memory on device
	}

	// print statistics
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "Prism       : %d\n", prisms.size());
		fprintf(stderr, "Triangles   : %d\n", prisms.size() * 8);
		fprintf(stderr, "GPU Blocks  : %d\n", blocks);
		fprintf(stderr, "GPU Threads : %d\n", threads);
		fprintf(stderr, "Runtime_GPU : %f s\n", runtimeGpu / 1000.0);
		fprintf(stderr, "Runtime_CPU : %f s\n", runtimeCpu / 1000.0);
		fprintf(stderr, "\n");
	}
	// Cleanup;
	{
		cudaFreeHost(hPrisms);
		cudaFreeHost(hVertices);
		cudaFree(dPrisms);
		cudaFree(dVertices);

		// Cleanup mex Variables
		cudaFree(p_in);
		cudaFree(n_x);
		cudaFree(n_y);
		cudaFree(neighbors);
		cudaFree(forbidden);
		cudaFree(n_p);
		cudaDeviceReset();
	}

	return 0;
}


//----------------------------------------------------
// Auxillary function definition
//----------------------------------------------------

float4 toBarycentric(TriangleCu t, PointCu p){
	float x1,x2,x3, y1,y2,y3, x,y;
	float4 b;

	x1 = t.A.x;
	x2 = t.B.x;
	x3 = t.C.x;

	y1 = t.A.y;
	y2 = t.B.y;
	y3 = t.C.y;

	x = p.x;
	y = p.y;

	b.x = ((y2-y3)*(x-x3)+(x3-x2)*(y-y3)) / ((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
	b.y = ((y3-y1)*(x-x3)+(x1-x3)*(y-y3)) / ((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
	b.z = 1 - b.x - b.y;
	b.w = 0;

	// In case of division by 0 --> nan
	if((fabs((b.x + b.y + b.z) - 1)) != (fabs((b.x + b.y + b.z) - 1)))
		b.z = 2;
	return b;
}

/**
  @brief Detects collisions of triangle and point with
  precondition, that the point is on the same 
  plane as the point.
 **/
bool collide(TriangleCu t, PointCu p){
	float4 b = toBarycentric(t, p);
	return (b.x > 0) && (b.x < 1) && (b.y > 0) && (b.y < 1) && (b.z > 0) && (b.z < 1) && (b.z == b.z);
}


/**
  @brief Detects collisions of a triangle and a ray without
  a precondition.
 **/
bool collide(TriangleCu t, RayCu r){
	PlaneCu pl;
	float b1, b2, b3, c1, c2, c3;

	b1 = t.B.x;
	b2 = t.B.y;
	b3 = t.B.z;

	c1 = t.C.x;
	c2 = t.C.y;
	c3 = t.C.z;

	pl.P = t.A;
	pl.normal.x = (b2*c3 - b3*c2);
	pl.normal.y = (b3*c1 - b1*c3);
	pl.normal.z = (b1*c2 - b2*c1);

	return collide(t, intersection(pl, r));
}

bool collide(PrismCu pr, RayCu r){
	bool hasCollide;
	PointCu A1 = pr.t1.A;
	PointCu B1 = pr.t1.B;
	PointCu C1 = pr.t1.C;
	PointCu A2 = {pr.t1.A.x, pr.t1.A.y, pr.t1.A.w, 1};
	PointCu B2 = {pr.t1.B.x, pr.t1.B.y, pr.t1.B.w, 1};
	PointCu C2 = {pr.t1.C.x, pr.t1.C.y, pr.t1.C.w, 1};

	TriangleCu triangles[8] = {
		pr.t1,
		{A2, B2, C2},
		{A1, B1, A2},
		{B1, B2, A2},
		{B1, C1, C2},
		{B1, B2, C2},
		{A1, C1, C2},
		{A1, A2, C2}};

	hasCollide = 
		collide(triangles[0], r)
		|| collide(triangles[1], r)
		|| collide(triangles[2], r) 
		|| collide(triangles[3], r)
		|| collide(triangles[4], r) 
		|| collide(triangles[5], r) 
		|| collide(triangles[6], r) 
		|| collide(triangles[7], r);

	return hasCollide;
}

/**
  @brief Intersection calculates the intersection between a plane p
  and a ray r. There is no detection for rays in the plane
  or for parallel plane. 

  It uses the normal of the plane to derive the coordinate form 
  of the plane. With the help of a coordinate form it is very
  easy to get the intersection point between a ray and a plane.

  ray   g: y~ = x~ + t*p~
  plane E: y~ = a~ + r*b~ + s*c~
  d  = n1*(x1+t*p1) + n2*(x2+t*p2) + n3*(x3+t*p3)
  d  = n~ * a~
 **/
PointCu intersection(PlaneCu pl, RayCu r){
	PointCu intersectionPoint = {0.0,0.0,0.0};

	float t, d;

	// vector coordinates
	float n1, n2, n3, x1, x2, x3, p1, p2, p3, a1, a2, a3;

	// just get the coordinates from the structs
	n1 = pl.normal.x;
	n2 = pl.normal.y;
	n3 = pl.normal.z;

	a1 = pl.P.x;
	a2 = pl.P.y;
	a3 = pl.P.z;

	x1 = r.P.x;
	x2 = r.P.y;
	x3 = r.P.z;

	p1 = r.direction.x;
	p2 = r.direction.y;
	p3 = r.direction.z;

	// calculation of intersection
	d = n1*a1 + n2*a2 + n3*a3;
	t = (d - n1*x1 - n2*x2 - n3*x3) / (n1*p1 + n2*p2 + n3*p3);

	intersectionPoint.x = x1 + t * p1;
	intersectionPoint.y = x2 + t * p2;
	intersectionPoint.z = x3 + t * p3;

	return intersectionPoint;

}

float distance(point a, point b){
	float d = sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y),2) + pow((b.z - a.z),2));
	return fabs(d);
}

std::vector<TriangleCu> generateTriangles(int height, int weight, float level){
	int h,w;
	std::vector<TriangleCu> triangles;
	for(h = 0; h < height; ++h){
		for(w = 0; w < weight; ++w){
			TriangleCu t1 = {
				{float(h), float(w), level, 1},
				{float(h), float(w+1), level, 1},
				{float(h+1), float(w), level, 1}};
			TriangleCu t2 = {
				{float(h), float(w+1), level, 1},
				{float(h+1), float(w+1), level, 1},
				{float(h+1), float(w), level, 1}};
			triangles.push_back(t1);
			triangles.push_back(t2);

		}

	}

	return triangles;
}

std::vector<PrismCu> generatePrisms(int height, int weight, float level){
	int h,w,l;
	std::vector<PrismCu> prisms;
	for(l = 0; l < level; ++l){
		for(h = 0; h < height; ++h){
			for(w = 0; w < weight; ++w){
				TriangleCu a1 = {
					{float(h), float(w), l, l+1},
					{float(h), float(w+1), l, l+1},
					{float(h+1), float(w), l, l+1}};
				TriangleCu b1 = {
					{float(h), float(w+1), l, 1+1},
					{float(h+1), float(w+1), l, 1+1},
					{float(h+1), float(w), l, 1+1}};

				PrismCu pr1 = {a1};
				PrismCu pr2 = {b1};

				prisms.push_back(pr1);
				prisms.push_back(pr2);

			}

		}

	}

	return prisms;
}

RayCu generateRay(const int heigth, const int width, const int level){
	float randHeigth = float(rand() % heigth) + (rand() / (float) RAND_MAX);
	float randWidth  = float(rand() % width ) + (rand() / (float) RAND_MAX);
	float rand_level  = float(rand() % level ) + (rand() / (float) RAND_MAX);

	float dirX = (rand() / (float) RAND_MAX);
	float dirY = (rand() / (float) RAND_MAX);
	float dirZ = (rand() / (float) RAND_MAX);

	RayCu r = {
		{randHeigth, randWidth, rand_level, 1},
		{dirX, dirY, dirZ, 0}};
	return r;
}


std::vector<RayCu> generateRays(const int height, const int width, const int level, const unsigned maxRays){
	std::vector<RayCu> rays;
	unsigned ray_i;
	for(ray_i = 0; ray_i < maxRays; ++ray_i){
		RayCu ray = generateRay(height, width, level);
		rays.push_back(ray);
	}
	return rays;
}

void printPoint(point p){
	fprintf(stdout, "Point\n");
	fprintf(stdout, "x: %f\n", p.x);
	fprintf(stdout, "y: %f\n", p.y);
	fprintf(stdout, "z: %f\n", p.z);

}
std::vector<VertexCu> generateSamples(int height, int width, int level){
	std::vector<VertexCu> samplePoints;
	int h,w,l;
	for(l = 0; l <= level; ++l){
		for(h = 0; h <= height; ++h){
			for(w = 0; w <= width; ++w){

				VertexCu p = {{float(h), float(w), float(l)}};

				samplePoints.push_back(p);
			}
		}
	}
	return samplePoints;
}
