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
#define RESULT_MSG_LENGTH 4
#define SAMPLE_MSG_LENGTH 2

void mpiHead(std::vector<float> &results, 
	     std::vector<double> &mse,
	     std::vector<unsigned> &totalRays,
	     std::vector<float> &runtime,
	     unsigned numberOfComputeNodes,
	     int sampleRange){
  MPI_Status status;
  float res[RESULT_MSG_LENGTH] = {0,0,0,0};
  int sample_i[SAMPLE_MSG_LENGTH] = {0,0};
  unsigned finished = 0;
  sample_i[1] = sampleRange;

  while(finished < numberOfComputeNodes){
    MPI_Recv(res, RESULT_MSG_LENGTH, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    switch(status.MPI_TAG){
    case RUNTIME_TAG:
      runtime.push_back(res[0]);
      finished++;
      break;

    case RESULT_TAG:
      results.at((unsigned)res[0])   = res[1];
      mse.at((unsigned)res[0])       = res[2];
      totalRays.at((unsigned)res[0]) = (unsigned)res[3];
      //fileProgressBar(results.size(),"output/progress");
      break;

    case SAMPLE_REQUEST_TAG:
      if(sample_i[0] == (int)results.size()){
	int abortMPI[2] = {-1,-1};
	MPI_Send(abortMPI, SAMPLE_MSG_LENGTH, MPI_INT, status.MPI_SOURCE, SAMPLE_SEND_TAG, MPI_COMM_WORLD);
      }
      else{
	MPI_Send(sample_i, SAMPLE_MSG_LENGTH, MPI_INT, status.MPI_SOURCE, SAMPLE_SEND_TAG, MPI_COMM_WORLD);
	sample_i[0] = std::min(sampleRange + sample_i[0], (int) results.size());
	sample_i[1] = std::min(sampleRange + sample_i[1], (int) results.size());
	
      }
      break;

    default:
      break;

    }

  }

}

/**
 * This MPI-node make phiASE computations.
 * It will request a sample point range from
 * MPI head node, and send their results back
 * sequentially.
 *
 **/
void mpiCompute(unsigned &hostRaysPerSample,
		const unsigned maxRaysPerSample,
		const Mesh& dMesh,
		const Mesh& hMesh,
		const std::vector<double>& hSigmaA,
		const std::vector<double>& hSigmaE,
		const std::vector<float>& mseThreshold,
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
      res[0] = runtime;
      MPI_Send(res, RESULT_MSG_LENGTH, MPI_FLOAT, HEAD_NODE, RUNTIME_TAG, MPI_COMM_WORLD); 
      break;
    }
    else{

      calcPhiAse ( hostRaysPerSample,
		   maxRaysPerSample,
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

      for(int i=sample_i[0]; i < sample_i[1]; ++i){
	res[0] = i; 
	res[1] = hPhiAse.at(i);
	res[2] = mse.at(i);
	res[3] = totalRays.at(i);
	totalRuntime += runtime;
	MPI_Send(res, RESULT_MSG_LENGTH, MPI_FLOAT, HEAD_NODE, RESULT_TAG, MPI_COMM_WORLD); 
      }

    }

  }

}

float calcPhiAseMPI ( unsigned &hRaysPerSample,
		      const unsigned maxRaysPerSample,
		      const Mesh& dMesh,
		      const Mesh& hMesh,
		      const std::vector<double>& hSigmaA,
		      const std::vector<double>& hSigmaE,
		      const std::vector<float>& mseThreshold,
		      const bool useReflections,
		      std::vector<float> &hPhiAse,
		      std::vector<double> &mse,
		      std::vector<unsigned> &totalRays,
		      unsigned gpu_i,
		      unsigned maxSample_i){



  int mpiError = MPI_Init(NULL,NULL);
  if(mpiError != MPI_SUCCESS){
    std::cerr << "Error starting MPI program." << std::endl;
    MPI_Abort(MPI_COMM_WORLD,mpiError);
    return 1;
  }

  int rank;
  int size;
  float runtime;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<float> runtimes(size,0);

  switch(rank){
  case HEAD_NODE:
    mpiHead(hPhiAse, mse, totalRays, runtimes, size-1, ceil(maxSample_i / (float)(size-1)));
     // for(unsigned i = 0; i < hPhiAse.size(); ++i){
     //   dout(V_INFO) << i << " : " << hPhiAse.at(i) << std::endl;
     // }
    break;

  default:
    mpiCompute(hRaysPerSample,
	       maxRaysPerSample,
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


  MPI_Finalize();
  return runtime;
}


