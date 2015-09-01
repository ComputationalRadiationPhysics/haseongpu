/**
 * Copyright 2015 Erik Zenker, Carlchristian Eckert, Marius Melzer
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

#include <cmath> /* log10 */
#include <iomanip> /* std::setfill, std::setw */
#include <fstream> /* std::ofstream */
#include <algorithm> /* std::max */
#include <chrono> /* steady_clock, time_point */
#include <sstream> /* stringstream */

#include <boost/filesystem/fstream.hpp> /* fs::fstream */
#include <boost/filesystem/path.hpp> /* fs::path */

#include <progressbar.hpp>
#include <logging.hpp>

namespace fs = boost::filesystem;
namespace chr = std::chrono;

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


/**
* @brief prints a simplistic line of ascii-art
*
* @param stream the stream to which to write the output
* @param progress the progress of the whole process (i.e. 90, if the process ranges
*        from 0-100 and the progress is at 90%
*        progress must be <= length! (normalize progress in order to acheive this)
* @param length the length of the finished bar.
*        length must be at least as long as the progress!
*/
void printBar(std::ostream &stream, int progress, int length){
    for (int i = 0; i<progress; ++i){
        stream << "#";
    }
    for (int i = 0; i < length - progress; ++i){
        stream << " ";
    }
}


/** 
 * Writes the time in a human-friendly way
 *
 * The returned string will be formatted to display the time in terms of hours,
 * minutes and seconds. If the duration is small enough, some of those will be
 * omitted. A template argument allows to include also milliseconds
 * Example: 5000s -> 1h 23m 20s
 *          1000s -> 16m 40s
 *
 * @param time the duration to write
 */
std::string humanRepresentation(chr::duration<float> const time){

  typedef chr::duration<int, std::ratio<1,1>> seconds;
  typedef chr::duration<int, std::ratio<60,1>> minutes;
  typedef chr::duration<int, std::ratio<3600,1>> hours;

  auto const tSec  = chr::duration_cast<seconds>(time).count() % 60; 
  auto const tMin  = chr::duration_cast<minutes>(time).count() % 60; 
  auto const tHour = chr::duration_cast<hours>(time).count();

  std::stringstream ss; 

  if(tHour) ss << tHour << "h ";
  if(tMin)  ss << tMin  << "m "; 
  if(tSec)  ss << tSec  << "s";

  return ss.str();

}


void fancyProgressBar(const unsigned nTotal){

#ifdef _WIN32
  const int length = 16;
#else
  const int length = 50;
#endif

  //use maxNTotal to find the global maximum between multiple calling threads
  static unsigned maxNTotal = 0;
  static chr::time_point < chr::steady_clock::time_point::clock > startTime;
  static unsigned part = 0;
  
  //find the starting time of the whole progress
  if(part==0){ startTime=chr::steady_clock::now(); }
  maxNTotal = std::max(maxNTotal, nTotal);
  static const unsigned fillwidthPart = unsigned(1+log10(maxNTotal));
  static unsigned tic  = 0;
  auto now = chr::steady_clock::now();
  ++part;

  //limit the update intervall (not faster than every 35ms, since that would be madness)
  chr::duration<float> const timeSpent = now - startTime; 
  if(timeSpent.count() > 0.035*tic || part==maxNTotal){
    ++tic;

    const float percentage = float(part) / float(maxNTotal);
    const auto timeTotal = timeSpent/percentage;
    const auto timeRemaining = timeTotal-timeSpent;

    dout(V_PROGRESS | V_NOLABEL) << "\r";
    dout(V_PROGRESS) << "[";
#ifdef _WIN32
    printBar(dout(V_PROGRESS | V_NOLABEL), int(percentage*length), length);
#else
    printWave(dout(V_PROGRESS | V_NOLABEL), tic, int(percentage*length), length);
#endif
    dout(V_PROGRESS | V_NOLABEL) << "] ";

    dout(V_PROGRESS | V_NOLABEL) << std::setfill(' ') << std::setw(3) << int(percentage*100) << "%";
    dout(V_PROGRESS | V_NOLABEL) << " (" << std::setfill(' ') << std::setw(fillwidthPart) << part << "/" << maxNTotal << ")";
    dout(V_PROGRESS | V_NOLABEL) << " after " << humanRepresentation(timeSpent) ;
    dout(V_PROGRESS | V_NOLABEL) << " (" << humanRepresentation(timeTotal) << " total, " << humanRepresentation(timeRemaining) << " remaining)";
    dout(V_PROGRESS | V_NOLABEL) << std::flush;
  }
}
