#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>


int main(int argc, char** argv){

  using namespace std;

  int mpiError = MPI_Init(NULL,NULL);
  int rank;
  int size;
  vector<int> samplePoints(30);

  for(int i=0;i<samplePoints.size();++i){
    samplePoints[i] = i;
  }

  if(mpiError != MPI_SUCCESS){
   cerr << "Error starting MPI program." << endl;
   MPI_Abort(MPI_COMM_WORLD,mpiError);
   return 1;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 0){
    int buf[1];
    MPI_Status status;
    int sample_i[1];

    while(!samplePoints.empty()){
      MPI_Recv(buf,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      sample_i[0] = samplePoints.pop_back(); 
      MPI_Send(sample_i,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
    }

    for(unsigned i=0 ; i<size-1 ; ++i){
      MPI_Recv(buf,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      sample_i[0] = -1; 
      MPI_Send(sample_i,1,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
    }

  }else{
    int buf[1] = {0};
    while(true){
      MPI_Send(buf,1,MPI_INT,0,0,MPI_COMM_WORLD);
      MPI_Recv(buf,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      if(buf[0] == -1){
        MPI_Finalize();
        return 0;
      }else{
        ofstream myfile;
        myfile.open("~/mpi_out.txt",ios::ate);
        myfile << "rank: " << rank << "   buf received: " << buf[0] << endl;
        myfile.close();
      }
    }
  }


  return 0;

}
