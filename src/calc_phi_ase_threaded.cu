#include <calc_phi_ase_threaded.h>
#include <vector>
#include <iostream>
#include <pthread.h>

#include <mesh.h>
#include <calc_phi_ase.h>

struct calcDndtAseArgs 
{
  calcDndtAseArgs(unsigned &phostRaysPerSample,
		  const unsigned pmaxRaysPerSample,
      const unsigned pmaxRepetitions,
		  const Mesh& pmesh,
		  const Mesh& phostMesh,
		  const std::vector<double>& psigmaA,
		  const std::vector<double>& psigmaE,
		  const std::vector<float>& pMseThreshold,
		  const bool puseReflections,
		  std::vector<float> &pphiAse,
		  std::vector<double> &pMse,
      std::vector<unsigned> &pTotalRays,
		  unsigned pgpu_i,
		  unsigned pminSample_i,
		  unsigned pmaxSample_i,
		  float &pruntime): hostRaysPerSample(phostRaysPerSample),
				    maxRaysPerSample(pmaxRaysPerSample),
            maxRepetitions(pmaxRepetitions),
				    mesh(pmesh),
				    hostMesh(phostMesh),
				    sigmaA(psigmaA),
				    sigmaE(psigmaE),
				    mseThreshold(pMseThreshold),
				    useReflections(puseReflections),
				    phiAse(pphiAse),
				    mse(pMse),
            totalRays(pTotalRays),
				    gpu_i(pgpu_i),
				    minSample_i(pminSample_i),
				    maxSample_i(pmaxSample_i),
				    runtime(pruntime){

  }
  unsigned &hostRaysPerSample;
  const unsigned maxRaysPerSample;
  const unsigned maxRepetitions;
  const Mesh& mesh;
  const Mesh& hostMesh;
  const std::vector<double>& sigmaA;
  const std::vector<double>& sigmaE;
  const std::vector<float>& mseThreshold;
  const bool useReflections;
  std::vector<float> &phiAse;
  std::vector<double> &mse;
  std::vector<unsigned> &totalRays;
  unsigned gpu_i;
  unsigned minSample_i;
  unsigned maxSample_i;
  float &runtime;
};

void *entryPoint(void* arg){
  calcDndtAseArgs *a = (calcDndtAseArgs*) arg;
  calcPhiAse( a->hostRaysPerSample,
   	      a->maxRaysPerSample,
          a->maxRepetitions,
   	      a->mesh,
   	      a->hostMesh,
   	      a->sigmaA,
   	      a->sigmaE,
   	      a->mseThreshold,
   	      a->useReflections,
   	      a->phiAse,
   	      a->mse,
          a->totalRays,
   	      a->gpu_i,
   	      a->minSample_i,
   	      a->maxSample_i,
	      a->runtime);

  return arg;
}

pthread_t calcPhiAseThreaded( unsigned &hostRaysPerSample,
			      const unsigned maxRaysPerSample,
            const unsigned maxRepetitions,
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
			      float &runtime){
  calcDndtAseArgs *args = new calcDndtAseArgs(hostRaysPerSample,
					      maxRaysPerSample,
                maxRepetitions,
					      mesh,
					      hostMesh,
					      sigmaA,
					      sigmaE,
					      mseThreshold,
					      useReflections,
					      phiAse,
					      mse,
                totalRays,
					      gpu_i,
					      minSample_i,
					      maxSample_i,
					      runtime);

  pthread_t threadId;
  pthread_create( &threadId, NULL, entryPoint, (void*) args);
  return threadId;

}

void joinAll(std::vector<pthread_t> threadIds){
  for(unsigned i = 0; i < threadIds.size(); ++i){
    pthread_join(threadIds[i], NULL);
  }

}

