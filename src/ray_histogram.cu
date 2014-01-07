#include <iostream>
#include <iomanip>
#include <logging.h>
#include <vector>
#include <types.h>
#include <map>

void ray_histogram(const std::vector<unsigned> totalRays, const unsigned max, const double mseThreshold, const std::vector<double> mseValues){
  // length of the maximum number of samples (e.g. max==4210)
  int fillwidth = log10(max)+4;

  // maximum length of the filling bar
  unsigned maxLength=50;

  // necessary size of the histogram
  std::map<unsigned,unsigned> histGreen;
  std::map<unsigned,unsigned> histRed;
  // if the entry doesn't exist, create it
  for(unsigned j=0; j<totalRays.size() ; ++j){
    //std::map<unsigned,unsigned>::iterator it = hist.find(totalRays.at(j));
    if(histGreen.find(totalRays.at(j)) == histGreen.end()){
      histGreen.insert(std::pair<unsigned,unsigned>(totalRays.at(j),0));
      histRed.insert(std::pair<unsigned,unsigned>(totalRays.at(j),0));
    }
    if(mseValues.at(j) <= mseThreshold){
      histGreen.find(totalRays.at(j))->second++;
      //itG->second++;
    } else{
      histRed.find(totalRays.at(j))->second++;
      //itR->second++; 
    }
  }


  std::map<unsigned,unsigned>::iterator itG;
  std::map<unsigned,unsigned>::iterator itR;
  for(itG=histGreen.begin(), itR=histRed.begin(); itG!=histGreen.end(); ++itG, ++itR){
    dout(V_STAT) << std::setw(fillwidth) << std::setfill(' ') <<  itG->first << " (";
    dout(V_STAT | V_NOLABEL) << "\033[0;32m" << std::setw(log10(totalRays.size())+3) << itG->second << "x";
    dout(V_STAT | V_NOLABEL) << "\033[0m" << " / ";
    dout(V_STAT | V_NOLABEL) << "\033[0;31m" << std::setw(log10(totalRays.size())+3) << itR->second << "x";
    dout(V_STAT | V_NOLABEL) << "\033[0m" << "):";

    // set color = green
    dout(V_STAT | V_NOLABEL) << "\033[0;32m";
    for(unsigned j=0;j< ceil(maxLength*(float(itG->second)/totalRays.size())) ; ++j){
      dout(V_STAT | V_NOLABEL) << "#"; 
    }

    // set color = red
    dout(V_STAT | V_NOLABEL) << "\033[0;31m";
    for(unsigned j=0;j< ceil(maxLength*(float(itR->second)/totalRays.size())) ; ++j){
      dout(V_STAT | V_NOLABEL) << "#"; 
    }
    dout(V_STAT | V_NOLABEL) << std::endl;
  }
}
