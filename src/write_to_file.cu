/**
 * Copyright 2013 Erik Zenker, Carlchristian Eckert, Marius Melzer
 *
 * This file is part of HASEonGPU
 *
 * HASEonGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HASEonGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HASEonGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */


#include <fstream> /* ofstream */
#include <sstream>
#include <iomanip>

#include <write_to_file.hpp>

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
