#include <thread.h>
#include <vector>
#include <iostream>
#include <pthread.h>


#include <mesh.h>
#include <calc_dndt_ase.h>

struct testArgs {
  int a;
  int b;
  int c;
  int d;
};

struct calcDndtAseArgs 
{
  calcDndtAseArgs(unsigned &pthreads,
		  unsigned &pblocks,
		  unsigned &phostRaysPerSample,
		  const unsigned pmaxRaysPerSample,
		  const Mesh& pmesh,
		  const Mesh& phostMesh,
		  const std::vector<double>& psigmaA,
		  const std::vector<double>& psigmaE,
		  const float pexpectationThreshold,
		  const bool puseReflections,
		  std::vector<double> &pdndtAse,
		  std::vector<float> &pphiAse,
		  std::vector<double> &pexpectation,
		  unsigned pgpu_i,
		  unsigned pminSample_i,
		  unsigned pmaxSample_i): threads(pthreads),
					  blocks(pblocks),
					  hostRaysPerSample(phostRaysPerSample),
					  maxRaysPerSample(pmaxRaysPerSample),
					  mesh(pmesh),
					  hostMesh(phostMesh),
					  sigmaA(psigmaA),
					  sigmaE(psigmaE),
					  expectationThreshold(pexpectationThreshold),
					  useReflections(puseReflections),
					  dndtAse(pdndtAse),
					  phiAse(pphiAse),
					  expectation(pexpectation),
					  gpu_i(pgpu_i),
					  minSample_i(pminSample_i),
					  maxSample_i(pmaxSample_i){

  }
  unsigned &threads;
  unsigned &blocks; 
  unsigned &hostRaysPerSample;
  const unsigned maxRaysPerSample;
  const Mesh& mesh;
  const Mesh& hostMesh;
  const std::vector<double>& sigmaA;
  const std::vector<double>& sigmaE;
  const float expectationThreshold;
  const bool useReflections;
  std::vector<double> &dndtAse;
  std::vector<float> &phiAse;
  std::vector<double> &expectation;
  unsigned gpu_i;
  unsigned minSample_i;
  unsigned maxSample_i;
};

void *entryPoint(void* arg){
  calcDndtAseArgs *a = (calcDndtAseArgs*) arg;
  calcDndtAse(a->threads,
   	      a->blocks,
   	      a->hostRaysPerSample,
   	      a->maxRaysPerSample,
   	      a->mesh,
   	      a->hostMesh,
   	      a->sigmaA,
   	      a->sigmaE,
   	      a->expectationThreshold,
   	      a->useReflections,
   	      a->dndtAse,
   	      a->phiAse,
   	      a->expectation,
   	      a->gpu_i,
   	      a->minSample_i,
   	      a->maxSample_i);

  return arg;
}

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
			      unsigned maxSample_i)
{
  calcDndtAseArgs *a = new calcDndtAseArgs(threads,
		    blocks,
		    hostRaysPerSample,
		    maxRaysPerSample,
		    mesh,
		    hostMesh,
		    sigmaA,
		    sigmaE,
		    expectationThreshold,
		    useReflections,
		    dndtAse,
		    phiAse,
		    expectation,
		    gpu_i,
		    minSample_i,
		    maxSample_i);

  pthread_t threadId;
  pthread_create( &threadId, NULL, entryPoint, (void*) a);
  return threadId;

}

void joinAll(std::vector<pthread_t> threadIds){
  for(unsigned i = 0; i < threadIds.size(); ++i){
    pthread_join(threadIds[i], NULL);
  }

}

