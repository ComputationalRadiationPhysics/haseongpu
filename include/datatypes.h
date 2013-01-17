typedef float4 pointCu;
typedef float4 vectorCu;
typedef struct triangleCu{
	pointCu A;
	pointCu B;
	pointCu C;
} TRIANGLE_CU;

typedef struct prismCu{
	triangleCu t1;
//	float height; //OPTIMIZE: The height could be stored as 4th parameter of one of the Triangle-coordinates?
} PRISM_CU;

typedef struct planeCu {
	pointCu P;
	vectorCu normal;
} PLANE_CU;

// Describes one vertex of the input-Mesh
typedef struct vertexCu {
	pointCu P;		// the Position
//	float4 G;		// The ASE-Gain in this Point (values from the rays are added)

	// OPTIMIZE: distribute Writes of G over more than 1 position in this
	// variable (e.g. through modulo thread-ID)
	// -> could result in less concurrent write-operations
	// Alternatively, save G in 4th coordinate of P
} VERTEX_CU;

typedef struct rayCu {
	pointCu P;			// the random starting point. The ase_flux is saved as the 4th component
	vectorCu direction;  // the position of the vertexCu, where the ray is going to
	//float phiAse;		// the accumulated ASE-Flux for this ray. is now saved in P
} RAY_CU;
