#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <thread>


#define RESULT_TAG 1
#define SAMPLE_REQUEST_TAG 2

void waitForResults(std::vector<float> *phiASE){
  MPI_Status status;
  float res[2] = {0,0};
  unsigned resultsRecv = 0;
  while(resultsRecv < phiASE->size()){
    MPI_Recv(res, 2, MPI_FLOAT, MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
    phiASE->at((unsigned)res[0]) = res[1];
    resultsRecv++;

  }

}

void schedulePoint(int numberOfSamples){
  MPI_Status status;
  int buf[1] = {0};
  int sampleSend = 0;
  while(sampleSend < numberOfSamples){
    MPI_Recv(buf, 1, MPI_INT, MPI_ANY_SOURCE, SAMPLE_REQUEST_TAG, MPI_COMM_WORLD, &status);
    buf[0] = sampleSend;
    MPI_Send(buf, 1, MPI_INT, status.MPI_SOURCE, SAMPLE_REQUEST_TAG, MPI_COMM_WORLD);
    sampleSend++;
  }

}

void testThread(){
  std::cout << "test" << std::endl;
}

int main(int argc, char** argv){

  using namespace std;


  int rank;
  int size;
  std::vector<float> *asePhi = new std::vector<float>(30,0.0);

  int mpiError = MPI_Init(NULL,NULL);
  if(mpiError != MPI_SUCCESS){
   cerr << "Error starting MPI program." << endl;
   MPI_Abort(MPI_COMM_WORLD,mpiError);
   return 1;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Manager Node
  if(rank == 0){
    //std::thread resultThread(waitForResults, asePhi);
    std::thread scheduleThread(schedulePoint, asePhi->size());

    resultThread.join();
    scheduleThread.join();

    // Finish compute nodes
    for(int i=0 ; i < size-1 ; ++i){
      int buf[1] = {0};
      MPI_Status status;
      MPI_Recv(buf, 1, MPI_INT, MPI_ANY_SOURCE, SAMPLE_REQUEST_TAG, MPI_COMM_WORLD, &status);
      int finished = -1; 
      MPI_Send(&finished, 1, MPI_INT, status.MPI_SOURCE, SAMPLE_REQUEST_TAG, MPI_COMM_WORLD);
    }

  }
  // Compute Node
  else{
    while(true){
      MPI_Status status;
      int sample_i = 0;
      int buf[1];
      MPI_Send(buf, 1, MPI_INT, 0, SAMPLE_REQUEST_TAG, MPI_COMM_WORLD);
      MPI_Recv(buf, 1, MPI_INT, 0, SAMPLE_REQUEST_TAG, MPI_COMM_WORLD, &status);
      sample_i = buf[0];
      if(sample_i == -1){
        break;
      }
      else{
	float res[2] = {0.0};
	res[0] = sample_i;
	res[1] = rank;
	MPI_Send(res, 2, MPI_FLOAT, 0, RESULT_TAG, MPI_COMM_WORLD);
        cout << "rank: " << rank << "   buf received: " << sample_i << endl;

      }
    }
  }

  MPI_Finalize();
  return 0;

}
