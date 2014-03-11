#include "parser.h"
#include <string> /* string */
#include <vector> /* vector */
#include <logging.h> 
#include <types.h>

void parseCommandLine(
    const int argc,
    char** argv,
    unsigned *raysPerSample,
    unsigned *maxRaysPerSample,
    std::string *inputPath,
    bool *writeVtk,
    std::string *compareLocation,
    RunMode *mode,
    bool *useReflections,
    unsigned *maxgpus,
    int *minSample_i,
    int *maxSample_i,
    unsigned *maxRepetitions,
    std::string *outputPath,
    double *mseThreshold
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

    if(p.first == "--input") {
      std::string temp_input(p.second);

      // Add slash at the end, if missing
      if ((temp_input)[temp_input.size() - 1] == 'w')
        temp_input.erase(temp_input.size() - 1, 1);
      else if (temp_input[temp_input.size() - 1] != '/')
        temp_input.append("/");

      *inputPath = temp_input;
    }

    if( p.first =="--output") {

      std::string temp_output(p.second);

      // Add slash at the end, if missing
      if ((temp_output)[temp_output.size() - 1] == 'w')
        temp_output.erase(temp_output.size() - 1, 1);
      else if (temp_output[temp_output.size() - 1] != '/')
        temp_output.append("/");

      *outputPath = temp_output;
    }

    if (p.first == "--write-vtk") {
      *writeVtk = true;
    }

    // Parse what vtk file to compare with
    if (p.first == "--compare") {
      *compareLocation = p.second;
    }

    if (p.first == "--runmode") {
      if (p.second == "threaded")
        *mode = GPU_THREADED;
      if (p.second == "cpu")
        *mode = CPU;
      if (p.second == "mpi")
        *mode = GPU_MPI;

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

    if (p.first == "--verbosity"){
      verbosity = unsigned(atoi(p.second.c_str()));
    }

    if(p.first == "--repetitions"){
      *maxRepetitions = unsigned(atoi(p.second.c_str()));
    }

    if(p.first == "--mse-threshold"){
      *mseThreshold = float(atof(p.second.c_str()));
    }

  }
}

int checkParameterValidity(
    const int argc,
    const unsigned raysPerSample,
    unsigned *maxRaysPerSample,
    const std::string inputPath,
    const unsigned deviceCount,
    const RunMode mode,
    unsigned *maxgpus,
    const int minSampleRange,
    const int maxSampleRange,
    const unsigned maxRepetitions,
    const std::string outputPath,
    double *mseThreshold
    ) {

  if (argc <= 1) {
    dout(V_ERROR) << "No commandline arguments found" << std::endl;
    dout(V_ERROR) << "Usage: ./calcPhiASE ARGS [OPTARGS]" << std::endl;
    dout(V_ERROR) << "" << std::endl;
    dout(V_ERROR) << "ARGS (required)" << std::endl;
    dout(V_ERROR) << "  --runmode=RUNMODE" << std::endl;
    dout(V_ERROR) << "  --rays=RAYS" << std::endl;
    dout(V_ERROR) << "  --maxrays=MAXRAYS" << std::endl;
    dout(V_ERROR) << "  --input=FOLDER" << std::endl;
    dout(V_ERROR) << "  --output=FOLDER" << std::endl;
    dout(V_ERROR) << "" << std::endl;
    dout(V_ERROR) << "OPTARGS (optional)" << std::endl;
    dout(V_ERROR) << "  --maxgpus=N" << std::endl;
    dout(V_ERROR) << "  --min_sample_i=<index of first sample>" << std::endl;
    dout(V_ERROR) << "  --max_sample_i=<index of last sample>" << std::endl;
    dout(V_ERROR) << "  --compare=<location of vtk-file to compare with>" << std::endl;
    dout(V_ERROR) << "  --verbosity=VERBOSITY_LEVEL" << std::endl;
    dout(V_ERROR) << "  --repetitions=MAX_REPETITIONS" << std::endl;
    dout(V_ERROR) << "  --mse_threshold=THRESHOLD" << std::endl;
    dout(V_ERROR) << "  --reflection" << std::endl;
    dout(V_ERROR) << "  --write-vtk" << std::endl;
    dout(V_ERROR) << "" << std::endl;
    dout(V_ERROR) << "Runmodes : cpu" << std::endl;
    dout(V_ERROR) << "           threaded" << std::endl;
    dout(V_ERROR) << "           mpi" << std::endl;
    dout(V_ERROR) << "Verbosity levels: 0 (quiet)" << std::endl; 
    dout(V_ERROR) << "                  1 (error)" << std::endl; 
    dout(V_ERROR) << "                  2 (warning)" << std::endl; 
    dout(V_ERROR) << "                  4 (info)" << std::endl; 
    dout(V_ERROR) << "                  8 (statistics)" << std::endl; 
    dout(V_ERROR) << "                 16 (debug)" << std::endl; 
    dout(V_ERROR) << "                 32 (progressbar)" << std::endl; 
    dout(V_ERROR) << "" << std::endl; 
    dout(V_ERROR) << "Please see README for more details!" << std::endl; 
    return 1;
  }
  if (mode == NONE) {
    dout(V_ERROR) << "Please specify the runmode with --runmode=MODE" << std::endl;
    return 1;
  }

  if (raysPerSample == 0) {
    dout(V_ERROR) << "Please specify the number of rays per sample Point with --rays=RAYS" << std::endl;
    return 1;
  }

  if (inputPath.size() == 0) {
    dout(V_ERROR) << "Please specify the experiment's location with --input=PATH_TO_EXPERIMENT" << std::endl;
    return 1;
  }

  if (outputPath.size() == 0) {
    dout(V_ERROR) << "Please specify the output location with --output=PATH_TO_EXPERIMENT" << std::endl;
    return 1;
  }


  if(*maxRaysPerSample < raysPerSample){
    dout(V_WARNING) << "maxRays < raysPerSample. Increasing maxRays to " << raysPerSample << " (will be non-adaptive!)" << std::endl;
    *maxRaysPerSample = raysPerSample;
  }

  if(*maxgpus > deviceCount){ dout(V_ERROR) << "You don't have so many devices, use --maxgpus=" << deviceCount << std::endl;
    return 1;
  }

  if(*maxgpus == 0){
    *maxgpus = deviceCount;
  }

  if(maxSampleRange < minSampleRange){
    dout(V_ERROR) << "maxSampleRange < minSampleRange!" << std::endl;
    return 1;
  }

  int samplesForNode = maxSampleRange-minSampleRange+1;
  if(samplesForNode < int(*maxgpus) && (minSampleRange != -1 || maxSampleRange != -1)){
    dout(V_WARNING) << "More GPUs requested than there are sample points. Number of used GPUs reduced to " << samplesForNode << std::endl;
     *maxgpus = unsigned(samplesForNode);
  }

  if(verbosity >= 64){
    verbosity = 63;
    dout(V_WARNING) << "Verbosity level should be between 0 (quiet) and 63 (all). Levels can be bitmasked together." << std::endl;
  }
  if(maxRepetitions < 1){
    dout(V_ERROR) << "At least 1 repetition is necessary!" << std::endl;
  }

  if(*mseThreshold == 0){
    *mseThreshold = 1000;
  }

  return 0;
}

void checkSampleRange(int* minSampleRange, int* maxSampleRange, const unsigned numberOfSamples){
  if(*maxSampleRange >= int(numberOfSamples)){
    dout(V_ERROR) << "maxSample_i is out of range! (There are only " << numberOfSamples << " samples)";
  }
  if(*minSampleRange == -1 && *maxSampleRange== -1){
    dout(V_WARNING) << "minSample_i/maxSample_i not set! Assuming a sample range of " << std::endl;
    dout(V_WARNING) << "0 to " << numberOfSamples-1 << std::endl;
    *minSampleRange = 0;
    *maxSampleRange = numberOfSamples-1;
  }
}
