#include "write_value_to_file.h"
#include <fstream> /* ofstream */
#include <sstream>
#include <iomanip>

//template <typename T>
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
