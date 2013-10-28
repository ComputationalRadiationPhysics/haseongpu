#ifndef THREAD_H
#define THREAD_H

#include <mesh.h>

#include <pthread.h>
#include <vector>


pthread_t calcPhiAseThreaded( unsigned &hostRaysPerSample,
			      const unsigned maxRaysPerSample,
			      const Mesh& mesh,
			      const Mesh& hostMesh,
			      const std::vector<double>& sigmaA,
			      const std::vector<double>& sigmaE,
			      const std::vector<float>& mseThreshold,
			      const bool useReflections,
			      std::vector<float> &phiAse,
			      std::vector<double> &mse,
            std::vector<unsigned> &totalRays,
			      unsigned gpu_i,
			      unsigned minSample_i,
			      unsigned maxSample_i,
			      float &runtime);

void joinAll(std::vector<pthread_t> threadIds);

#endif /* THREAD_H */
