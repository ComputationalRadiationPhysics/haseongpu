#include <calc_phi_ase_mpi.h>
#include <calc_phi_ase.h>
#include <mesh.h>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <logging.h>
#include <progressbar.h>
#include <algorithm>

// Nodes
#define HEAD_NODE 0

// Tags
#define RESULT_TAG 1
#define SAMPLE_REQUEST_TAG 2
#define SAMPLE_SEND_TAG 3
#define RUNTIME_TAG 4

// MSG SIZES
#define RESULT_MSG_LENGTH 5
#define SAMPLE_MSG_LENGTH 2

/**
 * @brief Manages communication with compute nodes.
 *        There are 4 possible messages :
 *        1. Compute node request sample
 *        2. Compute node sends results
 *        3. Head sends sample
 *        4. Head sends abort signal
 *
 * @param phiASE    return for phi ASE
 * @param mse       return for Mean squared error
 * @param totalRays return for raysPerSample
 */
void mpiHead(std::vector<float> &phiASE, 
	     std::vector<double> &mse,
	     std::vector<unsigned> &totalRays,
	     std::vector<float> &runtimes,
	     const Mesh& hMesh,
	     unsigned numberOfComputeNodes,
	     int sampleRange){
  MPI_Status status;
  float res[RESULT_MSG_LENGTH] = {0,0,0,0};
  int sample_i[SAMPLE_MSG_LENGTH] = {0,0};
  unsigned finishedComputeNodes = 0;
  unsigned sampleOffset = 0;
  sample_i[1] = sampleRange;

  while(finishedComputeNodes < numberOfComputeNodes){
    MPI_Recv(res, RESULT_MSG_LENGTH, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    switch(status.MPI_TAG){
      // Compute node finished and sends runtime
    case RUNTIME_TAG:
      runtimes.push_back(res[0]);
      finishedComputeNodes++;
      break;

      // Compute node sends its results for sample_i
    case RESULT_TAG:
      /**
       * res[0] : sample_i
       * res[1] : phiASE
       * res[2] : mse
       * res[3] : totalRays 
       **/
      sampleOffset = (unsigned)(res[0]);
      phiASE.at(sampleOffset)    = res[1];
      mse.at(sampleOffset)       = res[2];
      totalRays.at(sampleOffset) = (unsigned)res[3];
      fancyProgressBar(hMesh.numberOfSamples);
      break;

      // Compute node requests new sample point for computation
    case SAMPLE_REQUEST_TAG:
      if(sample_i[0] == (int)hMesh.numberOfSamples){
	// No sample points left, abort computation
	int abortMPI[2] = {-1,-1};
	MPI_Send(abortMPI, SAMPLE_MSG_LENGTH, MPI_INT, status.MPI_SOURCE, SAMPLE_SEND_TAG, MPI_COMM_WORLD);
      }
      else{
	// Send next sample range
	MPI_Send(sample_i, SAMPLE_MSG_LENGTH, MPI_INT, status.MPI_SOURCE, SAMPLE_SEND_TAG, MPI_COMM_WORLD);
	sample_i[0] = std::min(sample_i[0] + sampleRange, (int)hMesh.numberOfSamples); // min_sample_i
	sample_i[1] = std::min(sample_i[1] + sampleRange, (int)hMesh.numberOfSamples); // max_sample_i
	
      }
      break;

    default:
      break;

    }

  }

}

/**
 * @brief This MPI-node make phiASE computations.
 *        It will request a sample point range from
 *        MPI head node, and send their results back
 *        sequentially.
 *
 **/
void mpiCompute(unsigned &hostRaysPerSample,
		const unsigned maxRaysPerSample,
		const unsigned maxRepetitions,
		const Mesh& dMesh,
		const Mesh& hMesh,
		const std::vector<double>& hSigmaA,
		const std::vector<double>& hSigmaE,
		const double mseThreshold,
		const bool useReflections,
		std::vector<float> &hPhiAse,
		std::vector<double> &mse,
		std::vector<unsigned> &totalRays,
		unsigned gpu_i,
		float &runtime){
  while(true){
    MPI_Status status;
    int sample_i[SAMPLE_MSG_LENGTH] = {0,0};
    float totalRuntime = 0;
    float res[RESULT_MSG_LENGTH] = {0,0,0,0}; 
    MPI_Send(sample_i, SAMPLE_MSG_LENGTH, MPI_INT, HEAD_NODE, SAMPLE_REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(sample_i, SAMPLE_MSG_LENGTH, MPI_INT, HEAD_NODE, SAMPLE_SEND_TAG, MPI_COMM_WORLD, &status);


    if(sample_i[0] == -1){
      // Abort message received => send runtime
      res[0] = runtime;
      MPI_Send(res, RESULT_MSG_LENGTH, MPI_FLOAT, HEAD_NODE, RUNTIME_TAG, MPI_COMM_WORLD); 
      break;
    }
    else{
      // Sample range received => calculate
      calcPhiAse ( hostRaysPerSample,
		   maxRaysPerSample,
		   maxRepetitions,
		   dMesh,
		   hMesh,
		   hSigmaA,
		   hSigmaE,
		   mseThreshold,
		   useReflections,
		   hPhiAse,
		   mse,
		   totalRays,
		   gpu_i,
		   sample_i[0],
		   sample_i[1],
		   runtime);

      // Extract results and send it to head node
      for(int j=sample_i[0]; j < sample_i[1]; ++j){
	unsigned sampleOffset = (unsigned)(j);
	res[0] = (float)j; 
	res[1] = hPhiAse.at(sampleOffset);
	res[2] = mse.at(sampleOffset);
	res[3] = (float)totalRays.at(sampleOffset);
	totalRuntime += runtime;
	MPI_Send(res, RESULT_MSG_LENGTH, MPI_FLOAT, HEAD_NODE, RESULT_TAG, MPI_COMM_WORLD); 
      }

    }

  }

}

float calcPhiAseMPI ( unsigned &hRaysPerSample,
		      const unsigned maxRaysPerSample,
		      const unsigned maxRepetitions,
		      const Mesh& dMesh,
		      const Mesh& hMesh,
		      const std::vector<double>& hSigmaA,
		      const std::vector<double>& hSigmaE,
		      const double mseThreshold,
		      const bool useReflections,
		      std::vector<float> &hPhiAse,
		      std::vector<double> &mse,
		      std::vector<unsigned> &totalRays,
		      unsigned gpu_i){

  // Init MPI
  int mpiError = MPI_Init(NULL,NULL);
  if(mpiError != MPI_SUCCESS){
    dout(V_ERROR) << "Error starting MPI program." << std::endl;
    MPI_Abort(MPI_COMM_WORLD,mpiError);
    return 1;
  }

  // Get size and rank
  int rank;
  int size;
  float runtime;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<float> runtimes(size,0);

  // Rank 0 will be head node
  // all other ranks will be compute nodes
  switch(rank){
  case HEAD_NODE:
    mpiHead(hPhiAse, mse, totalRays, runtimes, hMesh, size-1, 1);
    cudaDeviceReset();   
    MPI_Finalize();
    break;


  default:
    // disable Information verbosity for other nodes than HEADNODE
    // (should have similar output anyway)
    verbosity &= ~V_PROGRESS;

    mpiCompute(hRaysPerSample,
	       maxRaysPerSample,
	       maxRepetitions,
	       dMesh,
	       hMesh,
	       hSigmaA,
	       hSigmaE,
	       mseThreshold,
	       useReflections,
	       hPhiAse,
	       mse,
	       totalRays,
	       gpu_i,
	       runtime);
    
    cudaDeviceReset();      
    MPI_Finalize();
    exit(0);
    break;
  };

  return size - 1;
}


