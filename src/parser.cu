#include "parser.h"
#include <string> /* string */
#include <vector> /* vector */
#include <stdio.h> /*fprintf*/

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
    int *sample_i
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
    fprintf(stderr, "arg[%d]: (%s,%s)\n", i, p.first.c_str(), p.second.c_str());

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

    }

    if (p.first == "--reflection"){
      *useReflections = true;
    }

    if (p.first == "--maxgpus"){
      *maxgpus = atoi(p.second.c_str());
    }

    if (p.first == "--sample_i"){
      *sample_i = atoi(p.second.c_str());
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
    unsigned *maxgpus
    ) {

  if (argc <= 1) {
    fprintf(stderr, "C No commandline arguments found\n");
    fprintf(stderr, "C Usage    : ./octrace --mode=[runmode]\n"); 
    fprintf(stderr, "C                      --rays=[number of rays]\n"); 
    fprintf(stderr, "C                      --experiment=[location to experiment-data]\n");
    fprintf(stderr, "C                      --compare=[location of vtk-file to compare with]\n");
    fprintf(stderr, "C                      --maxrays=[max number of rays for adaptive sampling]\n");
    fprintf(stderr, "C                      --maxgpus=[max number of gpus to use]\n");
    fprintf(stderr, "C Runmodes : for_loops\n");
    fprintf(stderr, "             ray_propagation_gpu\n");
    fprintf(stderr, "             test_environment\n");
    return 1;
  }
  if (mode == NONE) {
    fprintf(stderr, "C Please specify the runmode with --mode=\n");
    return 1;
  }
  if (raysPerSample == 0) {
    fprintf(stderr, "C Please specify the number of rays per sample Point with --rays=\n");
    return 1;
  }
  if (root.size() == 0) {
    fprintf(stderr, "C Please specify the experiment's location with --experiment=\n");
    return 1;
  }

  *maxRaysPerSample = max(raysPerSample,*maxRaysPerSample);

  if(*maxgpus > deviceCount){
    fprintf(stderr, "C Warning: You don't have so many devices, use --maxgpus=%d", deviceCount);
    return 1;
  }

  if(*maxgpus == 0){
    *maxgpus = deviceCount;
  }
  return 0;
}
