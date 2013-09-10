#ifndef PARSER_H
#define PARSER_H

#include <string>  /* string */
#include <iostream>
#include <fstream> /* ifstream */
#include <cstdlib> /* atof */
#include <vector> 

/**
 * @brief Parses a given file(filename) line by line.
 *        Each line should contain just one value
 *        and the value should be a number (short, unsigned,
 *        int, float, double).
 *
 * @param filename file to parse
 * @param vector contains the parsed values 
 *
 * @return 1 if parsing was succesful (file can be opened)
 * @return 0 otherwise
 **/
template <class T>
int fileToVector(const std::string filename, std::vector<T> *v){
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
    fprintf(stderr, "C Can't open file %s \n", filename.c_str());
    fileStream.close();
    return 1;
  }
  v->pop_back();
  fileStream.close();
  return 0;
  
}
/**
 * @brief Parses just one line(value) of a given file(filename).
 *        The value should be a number (short, unsigned,
 *        int, float, double).
 *
 * @param filename file to parse
 * @param value is the value which was parsed 
 *
 * @return 1 if parsing was succesful (file can be opened)
 * @return 0 otherwise
 **/
template <class T>
int fileToValue(const std::string filename, T &value){
  std::string line;
  std::ifstream fileStream; 

  fileStream.open(filename.c_str());
  if(fileStream.is_open()){
      std::getline(fileStream, line);
      value = (T) atof(line.c_str());
  }
  else{
    fprintf(stderr, "C Can't open file %s \n", filename.c_str());
    fileStream.close();
    return 1;
  }
  fileStream.close();
  return 0;
  
}


void parseCommandLine(
    const int argc,
    char** argv,
    unsigned *raysPerSample,
    unsigned *maxRaysPerSample,
    std::string *root,
    int *device,
    bool *silent,
    bool *writeVtk,
    std::string *compareLocation,
    int *mode,
    bool *useReflections,
    float *expectationThreshold
    );

int checkParameterValidity(
    const int argc,
    const unsigned raysPerSample,
    unsigned *maxRaysPerSample,
    const std::string root,
    int *device,
	  const unsigned deviceCount,
    const int mode,
    float *expectationThreshold
    );


#endif
