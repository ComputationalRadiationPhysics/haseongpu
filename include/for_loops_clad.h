#include <mesh.h>

#ifndef for_loops_clad_H
#define for_loops_clad_H
float forLoopsClad(
	std::vector<double> *dndtAse,
	unsigned &raysPerSample,
	Mesh *mesh,
	std::vector<double> *betaCellsVector,
	float hostNTot,
	float hostSigmaA,
	float hostSigmaE,
	unsigned hostNumberOfPoints,
	unsigned hostNumberOfTriangles,
	unsigned hostNumberOfLevels,
	float hostThicknessOfPrism,
	float hostCrystalFluorescence	);

#endif /* for_loops_clad_H */
