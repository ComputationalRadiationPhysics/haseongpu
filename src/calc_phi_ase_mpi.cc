#include <calc_phi_ase_mpi.h>
#include <calc_phi_ase.h>
#include <mesh.h>
#include <vector>
#include <mpi.h>
#include <iostream>
#include <logging.h>
#include <progressbar.h>

#define HEAD_NODE 0
#define RESULT_TAG 1
#define SAMPLE_REQUEST_TAG 2
#define SAMPLE_SEND_TAG 3

void mpiHead(std::vector<float> &results, unsigned numberOfComputeNodes){
  MPI_Status status;
  float res[2] = {0,0};
  int sample_i[1] = {0};
  unsigned finished = 0;
  while(finished < numberOfComputeNodes){
    MPI_Recv(res, 2, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    switch(status.MPI_TAG){
    case RESULT_TAG:
      //std::cout << "Received result" << std::endl;
      results.at((unsigned)res[0]) = res[1];
      fileProgressBar(results.size(),"output/progress");
      break;

    case SAMPLE_REQUEST_TAG:
      if(sample_i[0] == (int)results.size()){
	int abortMPI[1] = {-1};
	//std::cout << "Send abort" << std::endl;
	MPI_Send(abortMPI, 1, MPI_INT, status.MPI_SOURCE, SAMPLE_SEND_TAG, MPI_COMM_WORLD);
	finished++;
      }
      else{
	//std::cout << "Send sample point " << sample_i[0]<< " to "<< status.MPI_SOURCE << std::endl;
	MPI_Send(sample_i, 1, MPI_INT, status.MPI_SOURCE, SAMPLE_SEND_TAG, MPI_COMM_WORLD);
	sample_i[0]++;
      }
      break;

    default:
      break;

    }

  }

}

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
    int rank;
    MPI_Status status;
    int sample_i[1] = {0};
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //std::cout << "Rank[" << rank << "] send sample request" << std::endl;
    MPI_Send(sample_i, 1, MPI_INT, HEAD_NODE, SAMPLE_REQUEST_TAG, MPI_COMM_WORLD);
    MPI_Recv(sample_i, 1, MPI_INT, HEAD_NODE, SAMPLE_SEND_TAG, MPI_COMM_WORLD, &status);
    //std::cout << "Rank[" << rank << "] received sample: " << sample_i[0] << std::endl;

    if(sample_i[0] == -1){
      break;
    }
    else{
      float res[2] = {0,0}; 
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
		   sample_i[0]+1,
		   runtime);

      res[0] = sample_i[0]; 
      res[1] = hPhiAse.at(sample_i[0]);
      //std::cout << "Rank[" << rank << "] send result " << sample_i[0] << " : " << res[1] << std::endl;
      MPI_Send(res, 2, MPI_FLOAT, HEAD_NODE, RESULT_TAG, MPI_COMM_WORLD); 

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


  int rank;
  int size;
  float runtime;

  int mpiError = MPI_Init(NULL,NULL);
  if(mpiError != MPI_SUCCESS){
    std::cerr << "Error starting MPI program." << std::endl;
    MPI_Abort(MPI_COMM_WORLD,mpiError);
    return 1;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  switch(rank){
  case HEAD_NODE:
    std::cerr << "Start Head node" << std::endl;
    mpiHead(hPhiAse, size-1);
    for(unsigned i = 0; i < hPhiAse.size(); ++i){
      dout(V_INFO) << i << " : " << hPhiAse.at(i) << std::endl;
    }
    break;

  default:
    std::cerr << "Start Compute node" << std::endl;
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
    break;

  };


  MPI_Finalize();
  return runtime;
}


