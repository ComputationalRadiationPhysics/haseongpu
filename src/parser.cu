#include <string> /* string */
#include <vector> /* vector */
#include <stdio.h> /*fprintf*/

void parseCommandLine(
    const int argc,
    char** argv,
    unsigned *raysPerSample,
    std::string *root,
    int *device,
    bool *silent,
    std::string *compareLocation,
    int *mode
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

    if (p.first == "--experiment") {
      std::string temp_root(p.second);

      // Add slash at the end, if missing
      if ((temp_root)[temp_root.size() - 1] == 'w')
        temp_root.erase(temp_root.size() - 1, 1);
      else if (temp_root[temp_root.size() - 1] != '/')
        temp_root.append("/");

      *root = temp_root;
    }

    // Parse which cuda device to choose
    if (p.first == "--device") {
      *device = atoi(p.second.c_str());
    }

    // Parse if we want less output
    if (p.first == "--silent") {
      *silent = true;
    }

    // Parse what vtk file to compare with
    if (p.first == "--compare") {
      *compareLocation = p.second;
    }

    if (p.first == "--mode") {
      if (p.second == "ray_propagation_gpu")
        *mode = 0;
      if (p.second == "for_loops")
        *mode = 1;
    }
  }
}

int checkParameterValidity(
    int argc,
    unsigned raysPerSample,
    std::string root,
    int *device,
    int mode
    ) {

  if (argc <= 1) {
    fprintf(stderr, "C No commandline arguments found\n");
    fprintf(stderr, "C Usage    : ./octrace --mode=[runmode] --rays=[number of rays] --experiment=[location to experiment-data]\n");
    fprintf(stderr, "C Runmodes : for_loops\n");
    fprintf(stderr, "             ray_propagation_gpu\n");
    return 1;
  }
  if (mode == -1) {
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
  if (*device == -1) {
    *device = 0;
  }


  return 0;
}
