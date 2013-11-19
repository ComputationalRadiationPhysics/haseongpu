#include <iostream>
#include <iomanip>
#include <logging.h>
#include <vector>

void ray_histogram(const std::vector<unsigned> totalRays, const unsigned min, const unsigned max, const unsigned tooHighMse){
  int fillwidth = log10(max)+4;

  std::vector<unsigned> hist(log10(max/min)+1,0);
  for(unsigned j=0; j<totalRays.size() ; ++j){
    hist.at(log10(ceil(totalRays[j]/float(min))))++; 
  }
  int k=0;
  for(unsigned i=min ; i<=max ; i*=10){
    dout(V_STAT) << std::setw(fillwidth) << std::setfill(' ') <<  i << " (" << std::setw(log10(totalRays.size())+3) << hist.at(k) << "x):";

    dout(V_STAT | V_NOLABEL) << "\033[0;32m";
    if( i<max){
      for(unsigned j=0;j< ceil(50*(float(hist.at(k))/totalRays.size())); ++j){
        dout(V_STAT | V_NOLABEL) << "#"; 
      }
    }else{
      for(unsigned j=0;j< ceil(50*(float(hist.at(k)-tooHighMse)/totalRays.size())); ++j){
        dout(V_STAT | V_NOLABEL) << "#"; 
      }
      
      dout(V_STAT | V_NOLABEL) << "\033[0;31m";
      for(unsigned j=0;j< ceil(50*(float(tooHighMse)/totalRays.size())); ++j){
        dout(V_STAT | V_NOLABEL) << "#"; 
      }
    }
    dout(V_STAT | V_NOLABEL) << std::endl;
    k++;
  }


}
