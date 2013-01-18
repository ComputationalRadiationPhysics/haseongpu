#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector_types.h"
#include "assert.h"
#include <vector>
#include "curand_kernel.h"
#include "datatypes.h"

#define SMALL 1E-06
#define CUDA_CHECK_RETURN(value) {				\
	cudaError_t _mCudaStat = value;				\
	if (_mCudaStat != cudaSuccess) {				\
		fprintf(stderr, "Error %s at line %d in file %s\n",	\
				cudaGetErrorString(_mCudaStat), __LINE__, __FILE__);	\
		exit(1);							\
	}								\
}


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

__device__ RayCu generateRayGpu(PointCu vertexPoint, PrismCu startPrism, curandState randomstate){
	float u = curand_uniform(&randomstate);
	float v = curand_uniform(&randomstate);
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
	const float zRand = A.z + curand_uniform(&randomstate) * startPrism.t1.A.w;

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

__device__ float propagate(RayCu ray, PrismCu prisms[], PrismCu startprism){
	float gain = 1.f;
	float vecX = ray.direction.x - ray.P.x;
	float vecY = ray.direction.y - ray.P.y;
	float vecZ = ray.direction.z - ray.P.z;

	float distanceTotal = sqrt(vecX*vecX+vecY*vecY+vecZ*vecZ);
	float distance = distanceTotal;
	float length = distanceTotal;
	vecX /= distanceTotal;
	vecY /= distanceTotal;
	vecZ /= distanceTotal;

	PrismCu current = startprism;


	for(;;){
		length = distance;
		//generate the triangle surfaces of the prism
		const TriangleCu t1 = current.t1;
		const TriangleCu t2 = { 
			{t1.A.x, t1.A.y, t1.A.z + t1.A.w, 1},
			{t1.B.x, t1.B.y, t1.B.z + t1.B.w, 1},
			{t1.C.x, t1.C.y, t1.C.z + t1.C.w, 1}
		};

		// OPTIMIZE: make use of the rectangles!
		const TriangleCu surfaces[8] = {
			t1,
			t2,
			{t1.A, t1.B, t2.A},
			{t1.B, t2.B, t2.A},
			{t1.B, t1.C, t2.C},
			{t1.B, t2.B, t2.C},
			{t1.A, t1.C, t2.C},
			{t1.A, t2.A, t2.C}
		};

		int i=0;
		float lengthHelp = 0.f;
		for(i=0; i<8 ; ++i){ //OPTIMIZE: unroll, so that every surface can be optimized differently
			// get the generating vectors for the plane
			VectorCu AB = subtractPoints(surfaces[i].B, surfaces[i].A);
			VectorCu AC = subtractPoints(surfaces[i].C, surfaces[i].A);

			PlaneCu pl;
			pl.P = surfaces[i].A;
			// cross product of the vectors
			pl.normal.x = AB.y*AC.z - AB.z*AC.y;
			pl.normal.y = AB.z*AC.x - AB.x*AC.z;
			pl.normal.z = AB.x*AC.y - AB.y*AC.x;

			// direction * pl.normal
			float denominator = (ray.direction.x * pl.normal.x) + (ray.direction.y * pl.normal.y) + (ray.direction.z * pl.normal.z);
			float d = 0.f;
			float nominator = 0.f;
			if(denominator != 0.f) //OPTIMIZE: check if we have a lot of branch diversion, or if all threads behave the same
			{
				// A * pl.normal
				d = (surfaces[i].A.x * pl.normal.x) + (surfaces[i].A.y * pl.normal.y) + (surfaces[i].A.z * pl.normal.z);
				// d - (P * pl.normal)
				nominator = d - ((ray.P.x * pl.normal.x) + (ray.P.y * pl.normal.y) + (ray.P.z * pl.normal.y)); 
				lengthHelp = nominator/denominator;
				if(lengthHelp < length && lengthHelp > 0.f) //OPTIMIZE: most threads should do the same?
				{
					length = lengthHelp;
				}
			}
		}


		//with the new length, get the gain and add it
		// @TODO
		gain += length;

		// calculate values for next iteration
		distance -= length;
		if(abs(distance) < SMALL)
		{
			break;
		}

		ray.P.x += length*vecX;
		ray.P.y += length*vecY;
		ray.P.z += length*vecZ;

		//@TODO:
		// calculate the next PRISM (maybe with help of some neighbor-datastructure?

	}


	return gain;
}

__global__ void setupKernel ( curandState * state, unsigned long seed ){
	int id = threadIdx.x + blockDim.x*blockIdx.x;
	curand_init ( seed, id, 0, &state[id] );
	// OPTIMIZE: Use MersenneTwister or even a better PRNG
} 

// does the raytracing for a single ray (randomly generated) and a single (given) Vertex
//__global__ void raytraceStep( curandState* globalState, VertexCu vertex, PrismCu prisms[], int prismCount) {
__global__ void raytraceStep( VertexCu vertex, PrismCu prisms[], int prismCount) {
	//int id = threadIdx.x + blockDim.x*blockIdx.x;
	//curandState localState = globalState[id];

	//OPTIMIZE: the Octree should/could produce a subset of the prism-array!


	// this should give the same prism multiple times (so that every thread uses the same prism, which yields
	// big benefits for the memory access (and caching!)
	//PrismCu startprism = selectPrism(id, prisms, prismCount);	

	//RayCu ray = generateRayGpu(vertex.P,startprism, localState); //@TODO:verify
	//float initial_distance = distance(ray.P, ray.direction);

	//float gain = propagate(ray,prisms,startprism);

	//assert(fabs(gain-initial_distance) > 0.001);
		
	//atomicAdd(&(vertex.P.w),gain);
	//@TODO: find out, why atomic add won't compile
	//globalState[id] = localState;
}


//----------------------------------------------------
// Host Code
//----------------------------------------------------
int main(){

	//Variable definitions
	const unsigned maxRays = 100;
	const unsigned maxTriangles = 100;
	const unsigned maxVertices = 5;
	const unsigned length = ceil(sqrt(maxTriangles / 2));
	const unsigned depth  = 10;
	const unsigned maxPrisms = length * length * depth * 2;
	unsigned prism_i, vertex_i;
	float runtimeGpu = 0.0;
	float runtimeCpu = 0.0;
	cudaEvent_t start, stop;
	bool useGpu = true;
	curandState* devStates;

	// Generate testdata
	std::vector<VertexCu> vertices = generateSamples(length, length, depth);
	std::vector<PrismCu> prisms = generatePrisms(length, length, depth);
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// GPU Raytracing
	PrismCu* hPrisms, *dPrisms;
	VertexCu* hVertices, *dVertices;
	int threads = 256;
	int blocks = ceil(maxPrisms / threads);
	if(useGpu){

		//initialize memory
		{
			// Memory allocation on host
			CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&hPrisms, prisms.size() * sizeof(PrismCu), cudaHostAllocDefault));
			CUDA_CHECK_RETURN(cudaHostAlloc( (void**)&hVertices, vertices.size() * sizeof(VertexCu), cudaHostAllocDefault));

			// Memory initialisation on host
			for(prism_i = 0; prism_i < maxPrisms; ++prism_i){
				hPrisms[prism_i] = prisms[prism_i];
			}
			for(prism_i = 0; prism_i < vertices.size() ; ++prism_i){
				hVertices[prism_i] = vertices[prism_i];
			}


			// Memory allocation on device
			CUDA_CHECK_RETURN(cudaMalloc(&dPrisms, prisms.size() * sizeof(PrismCu)));
			CUDA_CHECK_RETURN(cudaMalloc(&dVertices, vertices.size() * sizeof(PrismCu)));
			CUDA_CHECK_RETURN(cudaMalloc(&devStates, threads*blocks*sizeof( curandState )));

			// Copy data from host to device
			CUDA_CHECK_RETURN(cudaMemcpy(dPrisms, hPrisms, maxPrisms * sizeof(PrismCu), cudaMemcpyHostToDevice));
			CUDA_CHECK_RETURN(cudaMemcpy(dVertices, hVertices, vertices.size() * sizeof(VertexCu), cudaMemcpyHostToDevice));
		}


		// Generating Random Numbers
	//	setupKernel<<< threads, blocks >>> ( devStates, time(NULL) );


		cudaEventRecord(start, 0);
		// start the Kernels
//		for(vertex_i = 0; vertex_i < maxVertices; ++vertex_i){
//			raytraceStep<<< threads, blocks >>> ( devStates , vertices[vertex_i] , dPrisms, maxPrisms);
//		}
		fprintf(stderr, "\nbetween the kernel");
		raytraceStep<<< threads, blocks >>> ( dVertices[0] , dPrisms, prisms.size() );

		fprintf(stderr, "\nafter the kernel");
		// Free memory on device
		cudaFree(devStates);

		// Evaluate device data
		{
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&runtimeGpu, start, stop);
		
		
		}
	}

	// print statistics
	{
	fprintf(stderr, "\n");
	fprintf(stderr, "Prism       : %d\n", maxPrisms);
	fprintf(stderr, "Triangles   : %d\n", maxPrisms * 8);
	fprintf(stderr, "Rays        : %d\n", maxRays);
	fprintf(stderr, "GPU Blocks  : %d\n", blocks);
	fprintf(stderr, "GPU Threads : %d\n", threads);
	fprintf(stderr, "Runtime_GPU : %f s\n", runtimeGpu / 1000.0);
	fprintf(stderr, "Runtime_CPU : %f s\n", runtimeCpu / 1000.0);
	fprintf(stderr, "\n");
	}
	// Cleanup
	cudaFreeHost(hPrisms);


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
