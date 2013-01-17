typedef float4 PointCu;
typedef float4 VectorCu;

typedef struct TriangleCu {
  PointCu A;
  PointCu B;
  PointCu C;
} TRIANGLE_CU;

// height is saved in every point of t1 as 4th coordinate (w)
typedef struct PrismCu {
  TriangleCu t1;
} PRISM_CU;

typedef struct PlaneCu {
  PointCu P;
  VectorCu normal;
} PLANE_CU;

// Describes one vertex of the input-Mesh
typedef struct VertexCu {
  PointCu P;		// the Position
  // OPTIMIZE: distribute Writes of G over more than 1 position in this
  // variable (e.g. through modulo thread-ID)
  // -> could result in less concurrent write-operations
  // Alternatively, save G in 4th coordinate of P
} VERTEX_CU;

typedef struct RayCu {
  PointCu P;			// the random starting point. The ase_flux is saved as the 4th component
  VectorCu direction;  // the position of the vertexCu, where the ray is going to
} RAY_CU;
