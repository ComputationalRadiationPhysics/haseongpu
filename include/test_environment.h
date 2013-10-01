#ifndef testEnvironment_T
#define testEnvironment_T

float testEnvironment (unsigned &hRaysPerSample,
		       const unsigned maxRaysPerSample,
		       const Mesh& dMesh,
		       const Mesh& hMesh,
		       const std::vector<double>& hSigmaA,
		       const std::vector<double>& hSigmaE,
		       const float expectationThreshold,
		       const bool useReflections,
		       std::vector<double> &dndtAse,
		       std::vector<float> &hPhiAse,
		       std::vector<double> &expectation
		       );

#endif /* testEnvironment_T */
