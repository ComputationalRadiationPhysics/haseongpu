#include "parser.h"
#include <string> /* string */
#include <vector> /* vector */
#include <stdio.h> /*fprintf*/
#include <logging.h> 

void parseCommandLine(
    const int argc,
    char** argv,
    unsigned *raysPerSample,
    unsigned *maxRaysPerSample,
    std::string *root,
    bool *silent,
    bool *writeVtk,
    std::string *compareLocation,
    RunMode *mode,
    bool *useReflections,
    unsigned *maxgpus,
    int *minSample_i,
    int *maxSample_i
    ) {

  std::vector<std::pair<std::string, std::string> > parameters;

  // Parse Commandline
  for (int i = 1; i < argc; ++i) {

    char* pos = strtok(argv[i], "=");
    std::pair < std::string, std::string > p(std::string(pos), std::string(""));
    pos = strtok(NULL, "=");
    if (pos != NULL) {
      p.second = std::string(pos);
    }
    parameters.push_back(p);
  }
  for (unsigned i = 0; i < parameters.size(); ++i) {
    std::pair < std::string, std::string > p = parameters.at(i);
    dout(V_INFO) << "arg[" << i << "]: (" << p.first << "," << p.second << ")" << std::endl;

    // Parse number of rays
    if (p.first == "--rays") {
      *raysPerSample = atoi(p.second.c_str());
    }

    if (p.first == "--maxrays"){
      *maxRaysPerSample = atoi(p.second.c_str());
    }

    if (p.first == "--experiment") {
      std::string temp_root(p.second);

      // Add slash at the end, if missing
      if ((temp_root)[temp_root.size() - 1] == 'w')
        temp_root.erase(temp_root.size() - 1, 1);
      else if (temp_root[temp_root.size() - 1] != '/')
        temp_root.append("/");

      *root = temp_root;
    }

    // Parse if we want less output
    if (p.first == "--silent") {
      *silent = true;
    }

    if (p.first == "--write-vtk") {
      *writeVtk = true;
    }

    // Parse what vtk file to compare with
    if (p.first == "--compare") {
      *compareLocation = p.second;
    }

    if (p.first == "--mode") {
      if (p.second == "ray_propagation_gpu")
        *mode = RAY_PROPAGATION_GPU;
      if (p.second == "for_loops")
        *mode = FOR_LOOPS;
      if (p.second == "test_environment")
        *mode = TEST;
      if (p.second == "mpi")
        *mode = RAY_PROPAGATION_MPI;

    }

    if (p.first == "--reflection"){
      *useReflections = true;
    }

    if (p.first == "--maxgpus"){
      *maxgpus = atoi(p.second.c_str());
    }

    if (p.first == "--min_sample_i"){
      *minSample_i = atoi(p.second.c_str());
    }
    if (p.first == "--max_sample_i"){
      *maxSample_i = atoi(p.second.c_str());
    }


  }
}

int checkParameterValidity(
    const int argc,
    const unsigned raysPerSample,
    unsigned *maxRaysPerSample,
    const std::string root,
    const unsigned deviceCount,
    const RunMode mode,
    unsigned *maxgpus,
    const int minSample_i,
    const int maxSample_i
    ) {

  if (argc <= 1) {
    dout(V_ERROR) << "No commandline arguments found" << std::endl;
    dout(V_ERROR) << "Usage    : ./octrace --mode=[runmode]" << std::endl;
    dout(V_ERROR) << "                     --rays=[number of rays]" << std::endl;
    dout(V_ERROR) << "                     --experiment=[location to experiment-data]" << std::endl;
    dout(V_ERROR) << "                     --compare=[location of vtk-file to compare with]" << std::endl;
    dout(V_ERROR) << "                     --maxrays=[max number of rays for adaptive sampling]" << std::endl;
    dout(V_ERROR) << "                     --maxgpus=[max number of gpus to use]" << std::endl;
    dout(V_ERROR) << "Runmodes : for_loops" << std::endl;
    dout(V_ERROR) << "           ray_propagation_gpu" << std::endl;
    dout(V_ERROR) << "           mpi" << std::endl;
    dout(V_ERROR) << "           test_environment" << std::endl;
    return 1;
  }
  if (mode == NONE) {
    dout(V_ERROR) << "Please specify the runmode with --mode=MODE" << std::endl;
    return 1;
  }
  if (raysPerSample == 0) {
    dout(V_ERROR) << "Please specify the number of rays per sample Point with --rays=RAYS" << std::endl;
    return 1;
  }
  if (root.size() == 0) {
    dout(V_ERROR) << "Please specify the experiment's location with --experiment=PATH_TO_EXPERIMENT" << std::endl;
    return 1;
  }

  *maxRaysPerSample = max(raysPerSample,*maxRaysPerSample);

  if(*maxgpus > deviceCount){
    dout(V_ERROR) << "You don't have so many devices, use --maxgpus=" << deviceCount << std::endl;
    return 1;
  }

  if(*maxgpus == 0){
    *maxgpus = deviceCount;
  }

  if(minSample_i < 0){
    dout(V_ERROR) << "--min_sample_i < 0!" << std::endl;
    return 1;
  }

  if(maxSample_i < minSample_i){
    dout(V_ERROR) << "maxSample_i < minSample_i!" << std::endl;
    return 1;
  }

  int samplesForNode = maxSample_i-minSample_i+1;
  if(samplesForNode < *maxgpus){
    dout(V_WARN) << "More GPUs requested than there are sample points. Number of used GPUs reduced to " << samplesForNode << std::endl;
     *maxgpus = samplesForNode;
  }
  return 0;
}
