#include <string>  /* string */
#include <iostream>
#include <fstream> /* ifstream */
#include <cstdlib> /* atof */
#include <vector> 

template <class T>
int fileToVector(std::string filename, std::vector<T> *v){
  std::string line;
  std::ifstream fileStream; 
  T value = 0.0;

  fileStream.open(filename.c_str());
  if(fileStream.is_open()){
    while(fileStream.good()){
      std::getline(fileStream, line);
      value = (T) atof(line.c_str());
      v->push_back(value);
    }

  }
  else{
    fprintf(stderr, "Can't open file %s \n", filename.c_str());
    fileStream.close();
    return 1;
  }
  v->pop_back();
  fileStream.close();
  return 0;
  
}

template <class T>
int fileToValue(std::string filename, T &value){
  std::string line;
  std::ifstream fileStream; 

  fileStream.open(filename.c_str());
  if(fileStream.is_open()){
      std::getline(fileStream, line);
      value = (T) atof(line.c_str());
  }
  else{
    fprintf(stderr, "Can't open file %s \n", filename.c_str());
    fileStream.close();
    return 1;
  }
  fileStream.close();
  return 0;
  
}




