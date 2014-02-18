#include "write_to_file.h"
#include <fstream> /* ofstream */
#include <sstream>
#include <iomanip>

int writeValueToFile(
    const float value, 
    const std::string path, 
    const std::string indexName1, 
    const int index1, 
    const std::string indexName2, 
    const int index2
    ){

  using namespace std;
  stringstream filenameStream;
  filenameStream << path << indexName1 << "_" << setfill('0') << setw(3) << index1 << "_" << indexName2 << setfill('0') << setw(6) << index2;

  ofstream oFile;
  oFile.open(filenameStream.str().c_str());
  oFile << value << endl;
  oFile.close();

  return 0;
}


void writeVectorToFile(std::vector<double> v, std::string pFilename){

  // Add time to filename
  time_t currentTime;
  time(&currentTime);
  std::stringstream filenameStream;
  filenameStream  << pFilename << "_" << (int) currentTime << ".dat";

  // Init filestream
  std::ofstream file;
  file.open(filenameStream.str().c_str());

  // Write vector data
  for(std::vector<double>::iterator it = v.begin(); it != v.end(); ++it){
    file << *it << std::endl;
  }
}
