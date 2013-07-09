#ifndef for_loops_clad_H
#define for_loops_clad_H
float forLoopsClad(
	std::vector<double> *dndtAse,
	unsigned &raysPerSample,
	std::vector<double> *betaValuesVector,
	std::vector<double> *xOfNormalsVector,
	std::vector<double> *yOfNormalsVector,
	std::vector<unsigned> *triangleIndicesVector,
	std::vector<int> *forbiddenVector,
	std::vector<int> *neighborsVector,
	std::vector<int> *positionsOfNormalsVector,
	std::vector<double> *pointsVector,
	std::vector<double> *betaCellsVector,
	std::vector<float> *surfacesVector,
	std::vector<double> *xOfTriangleCenterVector,
	std::vector<double> *yOfTriangleCenterVector,
	float hostNTot,
	float hostSigmaA,
	float hostSigmaE,
	unsigned hostNumberOfPoints,
	unsigned hostNumberOfTriangles,
	unsigned hostNumberOfLevels,
	float hostThicknessOfPrism,
	float hostCrystalFluorescence	);

#endif /* for_loops_clad_H */
