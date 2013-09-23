#ifndef THREAD_H
#define THREAD_H

#include <mesh.h>

#include <pthread.h>
#include <vector>


pthread_t calcDndtAseThreaded(unsigned &threads, 
			      unsigned &blocks, 
			      unsigned &hostRaysPerSample,
			      const unsigned maxRaysPerSample,
			      const Mesh& mesh,
			      const Mesh& hostMesh,
			      const std::vector<double>& sigmaA,
			      const std::vector<double>& sigmaE,
			      const float expectationThreshold,
			      const bool useReflections,
			      std::vector<double> &dndtAse,
			      std::vector<float> &phiAse,
			      std::vector<double> &expectation,
			      unsigned gpu_i,
			      unsigned minSample_i,
			      unsigned maxSample_i);

void joinAll(std::vector<pthread_t> threadIds);

#endif /* THREAD_H */
