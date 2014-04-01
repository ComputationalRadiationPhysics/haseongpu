/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASENonGPU
 *
 * HASENonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASENonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASENonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


#include <mpi.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <cstdlib>     /* system, NULL, EXIT_FAILURE */
#include <iomanip>      /*setfill, setw*/
#include <parser.h>


#define HEAD_NODE 0
#define RESULT_TAG 1
#define SAMPLE_REQUEST_TAG 2
#define SAMPLE_SEND_TAG 3

float runSimulation(int sample_i){
  std::stringstream command;
  std::stringstream outputFile;
  float result = 0.0;

  outputFile << "./output/results/wavelength_000_sample" << std::setfill('0') << std::setw(6) << sample_i;
  command << "./utils/start_exp.sh " << sample_i;
  if(system(command.str().c_str())){
    std::cerr << "Fail to run Simulation for sample " << sample_i << std::endl;
    return -1;
  }
  
  fileToValue(outputFile.str(), result);

  return result;
  
}

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

void mpiCompute(){
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
      res[0] = sample_i[0]; 
      res[1] = runSimulation(sample_i[0]);
      //std::cout << "Rank[" << rank << "] send result " << sample_i[0] << " : " << res[1] << std::endl;
      MPI_Send(res, 2, MPI_FLOAT, HEAD_NODE, RESULT_TAG, MPI_COMM_WORLD); 

    }

  }

}

int main(int argc, char** argv){
  using namespace std;

  int rank;
  int size;
  std::vector<float> phiASE(20,0.0);

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
    mpiHead(phiASE, size-1);
    for(unsigned i=0; i < phiASE.size(); ++i){
      std::cout << i << " " << phiASE.at(i) << std::endl;
    }
    break;

  default:
    mpiCompute();
    break;
  };


  //std::cout << "Rank[" << rank << "] finalize" << std::endl;
  MPI_Finalize();
  return 0;

}
