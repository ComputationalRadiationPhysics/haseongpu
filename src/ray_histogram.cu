#include <iostream>
#include <iomanip>
#include <logging.h>
#include <vector>
#include <types.h>
#include <map>

void ray_histogram_old(const std::vector<unsigned> totalRays, const unsigned min, const unsigned max, const unsigned tooHighMse){
  // length of the maximum number of samples (e.g. max==4210)
  int fillwidth = log10(max)+4;

  // maximum length of the filling bar
  unsigned maxLength=50;

  // necessary size of the histogram
  std::vector<unsigned> hist(log10(max/min)+1,0);

  // fill the histogram
  for(unsigned j=0; j<totalRays.size() ; ++j){
    // get logarithmic position in hist, based on every entry in totalRays 
    hist.at(log10(ceil(totalRays[j]/float(min))))++; 
  }
  int k=0;
  // print every line (e.g. 10k, 100k, 1M, ...)
  for(unsigned i=min ; i<=max ; i*=RAY_MULTIPLICATOR){
    dout(V_STAT) << std::setw(fillwidth) << std::setfill(' ') <<  i << " (" << std::setw(log10(totalRays.size())+3) << hist.at(k) << "x):";

    // set color = green
    dout(V_STAT | V_NOLABEL) << "\033[0;32m";
    // do this for every entry which can't have an exceeding MSE value
    if( i<max){
      for(unsigned j=0;j< ceil(maxLength*(float(hist.at(k))/totalRays.size())); ++j){
        dout(V_STAT | V_NOLABEL) << "#"; 
        if(j>maxLength){
          dout(V_ERROR | V_NOLABEL) << std::endl;
          dout(V_ERROR) << "Exceeded maxLength in Loop 1!" << std:: endl;
          break;
        }
      }

      // the last line can have exceeding MSE values -> print those in red
    }else{
      for(unsigned j=0;j< ceil(maxLength*(float(hist.at(k)-tooHighMse)/totalRays.size())) ; ++j){
        dout(V_STAT | V_NOLABEL) << "#"; 
        if(j>maxLength){
          dout(V_ERROR | V_NOLABEL) << std::endl;
          dout(V_ERROR) << "Exceeded maxLength in Loop 2!" << std:: endl;
          dout(V_ERROR) << "hist.at(" << k << ")= " << hist.at(k) << " tooHighMse=" << tooHighMse << std::endl;
          break;
        }
      }

      // set color = red for the rays which didn't hit the MSE threshold
      dout(V_STAT | V_NOLABEL) << "\033[0;31m";
      for(unsigned j=0;j< ceil(maxLength*(float(tooHighMse)/totalRays.size())) ; ++j){
        dout(V_STAT | V_NOLABEL) << "#"; 
        if(j>maxLength){
          dout(V_ERROR | V_NOLABEL) << std::endl;
          dout(V_ERROR) << "Exceeded maxLength in Loop 3!" << std:: endl;
          break;
        }
      }
    }
    dout(V_STAT | V_NOLABEL) << std::endl;
    k++;
  }

}

void ray_histogram(const std::vector<unsigned> totalRays, const unsigned max, const unsigned tooHighMse){
  // length of the maximum number of samples (e.g. max==4210)
  int fillwidth = log10(max)+4;

  // maximum length of the filling bar
  unsigned maxLength=50;

  // necessary size of the histogram
  std::map<unsigned,unsigned> hist;
  // if the entry doesn't exist, create it
  for(unsigned j=0; j<totalRays.size() ; ++j){
    std::map<unsigned,unsigned>::iterator it = hist.find(totalRays.at(j));
        if(it == hist.end()){
      hist.insert(std::pair<unsigned,unsigned>(totalRays.at(j),1));
    }else{
      it->second++; 
    }
  }

  int k=0;
  // print every line (e.g. 10k, 100k, 1M, ...)

  
  for(std::map<unsigned,unsigned>::iterator it = hist.begin(); it!=hist.end(); ++it){
    dout(V_STAT) << std::setw(fillwidth) << std::setfill(' ') <<  it->first << " (" << std::setw(log10(totalRays.size())+3) << it->second << "x):";

    // set color = green
    dout(V_STAT | V_NOLABEL) << "\033[0;32m";
    // do this for every entry which can't have an exceeding MSE value
    if( it->first < max){
      for(unsigned j=0;j< ceil(maxLength*(float(it->second)/totalRays.size())); ++j){
        dout(V_STAT | V_NOLABEL) << "#"; 
        if(j>maxLength){
          dout(V_ERROR | V_NOLABEL) << std::endl;
          dout(V_ERROR) << "Exceeded maxLength in Loop 1!" << std:: endl;
          break;
        }
      }

      // the last line can have exceeding MSE values -> print those in red
    }else{
      for(unsigned j=0;j< ceil(maxLength*(float(it->second - tooHighMse)/totalRays.size())) ; ++j){
        dout(V_STAT | V_NOLABEL) << "#"; 
        if(j>maxLength){
          dout(V_ERROR | V_NOLABEL) << std::endl;
          dout(V_ERROR) << "Exceeded maxLength in Loop 2!" << std:: endl;
          dout(V_ERROR) << "hist(" << it->first << "= " << it->second << " tooHighMse=" << tooHighMse << std::endl;
          break;
        }
      }

      // set color = red for the rays which didn't hit the MSE threshold
      dout(V_STAT | V_NOLABEL) << "\033[0;31m";
      for(unsigned j=0;j< ceil(maxLength*(float(tooHighMse)/totalRays.size())) ; ++j){
        dout(V_STAT | V_NOLABEL) << "#"; 
        if(j>maxLength){
          dout(V_ERROR | V_NOLABEL) << std::endl;
          dout(V_ERROR) << "Exceeded maxLength in Loop 3!" << std:: endl;
          break;
        }
      }
    }
    dout(V_STAT | V_NOLABEL) << std::endl;
    k++;
  }


}
