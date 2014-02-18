#include "progressbar.h"
#include <unistd.h>
#include <logging.h>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sys/time.h>


/**
 * @brief prints a line of ascii-art in the style of a sine-wave
 *
 * @param stream the stream to which to write the output
 * @param tic continuously increasing value related to the real outside time
 * @param progress the progress of the whole process (i.e. 90, if the process ranges
 *        from 0-100 and the progress is at 90%
 *        progress must be <= length! (normalize progress in order to acheive this)
 * @param length the length of the finished wave.
 *        length must be at least as long as the progress!
 */
void printWave(std::ostream &stream, unsigned tic, int progress, int length){
  for(int i=0;i<progress ;++i){
    switch((tic+i) % 12){
      case 0: stream << "ø"; break;
      case 1: stream << "¤"; break;
      case 2: stream << "º"; break;
      case 3: stream << "°"; break;
      case 4: stream << "`"; break;
      case 5: stream << "°"; break;
      case 6: stream << "º"; break;
      case 7: stream << "¤"; break;
      case 8: stream << "ø"; break;
      case 9: stream << ","; break;
      case 10: stream << "¸"; break;
      case 11: stream << ","; break;
    }
  }
  for(int i=0; i < length-progress ; ++i){
    stream << " ";
  }
}

void simpleProgressBar(unsigned part, unsigned full){
  unsigned length = 80;

  float percentage = (float(part)+1) / float(full);

  dout(V_PROGRESS | V_NOLABEL) << "\r";
  dout(V_PROGRESS) << "Progress: [";
  for(int i=0 ; i < (percentage*length) ; i++){
    dout(V_PROGRESS | V_NOLABEL) << "#";
  }
  for(int i=0;i< length-(percentage*length) ;i++){
    dout(V_PROGRESS | V_NOLABEL) << " ";
  }
  dout(V_PROGRESS | V_NOLABEL) << "] " << int(percentage*100) << "% (" << part+1 << "/" << full << std::flush;
}


/**
 * @brief give the difference of two time values in milliseconds
 *
 * @param start the starting time of the process to time
 * @param end the time when the process did finish
 * @return the difference (end-start) in milliseconds
 *
 */
unsigned long long timevalDiffInMillis(timeval start, timeval end){
  unsigned long long t = 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000;
  return t;
}


void fancyProgressBar(const unsigned nTotal){

  const int length = 50;

  //use maxNTotal to find the global maximum between multiple calling threads
  static unsigned maxNTotal = 0;
  static timeval startTime;
  static unsigned part = 0;
  
  //find the starting time of the whole progress
  if(part==0){ gettimeofday(&startTime,NULL); }
  maxNTotal = max(maxNTotal, nTotal);
  static const unsigned fillwidthPart = unsigned(1+log10(maxNTotal));
  static unsigned tic  = 0;
  timeval now;
  gettimeofday(&now,NULL);
  ++part;

  //limit the update intervall (not faster than every 35ms, since that would be madness)
  unsigned long long millisSpent = timevalDiffInMillis(startTime,now); 
  if(millisSpent > 35*tic || part==maxNTotal){
    ++tic;

    const float percentage = float(part) / float(maxNTotal);
    const float timeSpent = float(millisSpent) / 1000;
    const float timeTotal = timeSpent/percentage;
    const int timeRemaining = timeTotal-timeSpent;

    dout(V_PROGRESS | V_NOLABEL) << "\r";
    dout(V_PROGRESS) << "Progress: [";
    printWave(dout(V_PROGRESS | V_NOLABEL), tic, int(percentage*length), length);
    dout(V_PROGRESS | V_NOLABEL) << "] ";

    dout(V_PROGRESS | V_NOLABEL) << std::setfill(' ') << std::setw(3) << int(percentage*100) << "%";
    dout(V_PROGRESS | V_NOLABEL) << " (" << std::setfill(' ') << std::setw(fillwidthPart) << part << "/" << maxNTotal << ")";
    dout(V_PROGRESS | V_NOLABEL) << " after " << int(timeSpent) << "s";
    dout(V_PROGRESS | V_NOLABEL) << " (" << int(timeTotal) << "s total, " << timeRemaining << "s remaining)";
    dout(V_PROGRESS | V_NOLABEL) << std::flush;
  }
}

void fancyProgressBar(const unsigned current, const unsigned nTotal){

  const int length = 50;
  static unsigned maxNTotal = 0;
  static timeval startTime;
  static unsigned part = 0;

  // get the starting time on the very first call
  if(part==0){ gettimeofday(&startTime,NULL); }
  maxNTotal = max(maxNTotal, nTotal);
  static const unsigned fillwidthPart = unsigned(1+log10(maxNTotal));
  static unsigned tic  = 0;
  timeval now;
  gettimeofday(&now,NULL);
  part=current;

  //limit the update intervall (not faster than every 35ms, since this would be madness)
  unsigned long long millisSpent = timevalDiffInMillis(startTime,now); 
  if(millisSpent > 35*tic || part==maxNTotal){
    ++tic;

    const float percentage = float(part) / float(maxNTotal);
    const float timeSpent = float(millisSpent) / 1000;
    const float timeTotal = timeSpent/percentage;
    const int timeRemaining = timeTotal-timeSpent;

    dout(V_PROGRESS | V_NOLABEL) << "\r";
    dout(V_PROGRESS) << "Progress: [";
    printWave(dout(V_PROGRESS | V_NOLABEL), tic, int(percentage*length), length);
    dout(V_PROGRESS | V_NOLABEL) << "] ";

    dout(V_PROGRESS | V_NOLABEL) << std::setfill(' ') << std::setw(3) << int(percentage*100) << "%";
    dout(V_PROGRESS | V_NOLABEL) << " (" << std::setfill(' ') << std::setw(fillwidthPart) << part << "/" << maxNTotal << ")";
    dout(V_PROGRESS | V_NOLABEL) << " after " << int(timeSpent) << "s";
    dout(V_PROGRESS | V_NOLABEL) << " (" << int(timeTotal) << "s total, " << timeRemaining << "s remaining)";
    dout(V_PROGRESS | V_NOLABEL) << std::flush;
  }
}


void fileProgressBar(unsigned nTotal, std::string path){

  const int length = 50;
  std::ofstream filestream;
  static timeval startTime;
  static unsigned part = 0;
  //if this is the first call, set the start time 
  if(part==0){ gettimeofday(&startTime,NULL); }
  static const unsigned fillwidthPart = unsigned(1+log10(nTotal));
  static unsigned tic  = 0;
  timeval now;
  gettimeofday(&now,NULL);
  ++part;

  //limit the update intervall (not faster than every 35ms, since that would be madness)
  unsigned long long millisSpent = timevalDiffInMillis(startTime,now); 
  if(millisSpent > 35*tic || part==nTotal){
    ++tic;

    if(!filestream.is_open()){
      filestream.open(path.c_str(),std::ofstream::trunc);
    }

    const float percentage = float(part) / float(nTotal);
    const float timeSpent = float(millisSpent) / 1000;
    const float timeTotal = timeSpent/percentage;
    const int timeRemaining = timeTotal-timeSpent;
    const time_t finish = now.tv_sec + timeRemaining;

    // write progressbar
    filestream << "Progress: [";
    printWave(filestream,tic,int(percentage*length),length);
    filestream << "] ";

    // write progress in numbers
    filestream << std::setfill(' ') << std::setw(3) << unsigned(percentage*100) << "%";
    filestream << " (" << std::setfill(' ') << std::setw(fillwidthPart) << part << "/" << nTotal << ")" << std::endl;
    filestream << "running " << int(timeSpent) << "s";

    // write remaining time estimates
    filestream << " (" << int(timeTotal) << "s total, " << timeRemaining << "s remaining)" << std::endl;
    filestream << "estimated completion: " << ctime(&finish);

    filestream.close();
  }
}

