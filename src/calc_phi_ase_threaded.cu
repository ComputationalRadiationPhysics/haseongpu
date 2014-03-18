#include <calc_phi_ase_threaded.h>
#include <vector>
#include <iostream>
#include <pthread.h>

#include <mesh.h>
#include <calc_phi_ase.h>

struct calcDndtAseArgs 
{
  calcDndtAseArgs(const unsigned pminRaysPerSample,
		  const unsigned pmaxRaysPerSample,
		  const unsigned pmaxRepetitions,
		  const Mesh& pmesh,
		  const std::vector<double>& psigmaA,
		  const std::vector<double>& psigmaE,
		  const double pMseThreshold,
		  const bool puseReflections,
		  std::vector<float> &pphiAse,
		  std::vector<double> &pMse,
		  std::vector<unsigned> &pTotalRays,
		  unsigned pgpu_i,
		  unsigned pminSample_i,
		  unsigned pmaxSample_i,
		  float &pruntime): minRaysPerSample(pminRaysPerSample),
				    maxRaysPerSample(pmaxRaysPerSample),
				    maxRepetitions(pmaxRepetitions),
				    mesh(pmesh),
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
  const unsigned minRaysPerSample;
  const unsigned maxRaysPerSample;
  const unsigned maxRepetitions;
  const Mesh& mesh;
  const std::vector<double>& sigmaA;
  const std::vector<double>& sigmaE;
  const double mseThreshold;
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
  calcPhiAse( a->minRaysPerSample,
   	      a->maxRaysPerSample,
	      a->maxRepetitions,
   	      a->mesh,
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

pthread_t calcPhiAseThreaded( const unsigned minRaysPerSample,
			      const unsigned maxRaysPerSample,
			      const unsigned maxRepetitions,
			      const Mesh& mesh,
			      const std::vector<double>& sigmaA,
			      const std::vector<double>& sigmaE,
			      const double mseThreshold,
			      const bool useReflections,
			      std::vector<float> &phiAse,
			      std::vector<double> &mse,
			      std::vector<unsigned> &totalRays,
			      const unsigned gpu_i,
			      const unsigned minSample_i,
			      const unsigned maxSample_i,
			      float &runtime){
  calcDndtAseArgs *args = new calcDndtAseArgs(minRaysPerSample,
					      maxRaysPerSample,
					      maxRepetitions,
					      mesh,
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

