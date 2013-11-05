#include "write_value_to_file.h"
#include <fstream> /* ofstream */

template <typename T>
int writeValueToFile(const T value, const std::string filename){

  std::ofstream oFile;
  oFile.open(filename.c_str());
  oFile << value << std::endl;
  oFile.close();

  return 0;
}
