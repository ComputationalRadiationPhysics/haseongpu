#include "progressbar.h"
#include <unistd.h>
#include <ctime>
#include <logging.h>
#include <cmath>
#include <iomanip>
#include <fstream>

void simpleProgressBar(unsigned part, unsigned full){
	unsigned length = 80;

	float percentage = (float(part)+1) / float(full);

	dout(V_INFO | V_NOLABEL) << "\r";
	dout(V_INFO) << "Progress: [";
	for(int i=0 ; i < (percentage*length) ; i++){
		dout(V_INFO | V_NOLABEL) << "#";
	}
	for(int i=0;i< length-(percentage*length) ;i++){
		dout(V_INFO | V_NOLABEL) << " ";
	}
	dout(V_INFO | V_NOLABEL) << "] " << int(percentage*100) << "%% (" << part+1 << "/" << full << std::flush;
}

void fancyProgressBar(unsigned part, unsigned full, unsigned length, time_t starttime){

	float percentage = (float(part)+1) / float(full);

	dout(V_INFO | V_NOLABEL) << "\r";
	dout(V_INFO) << "Progress: [";
	for(int i=0 ; i < (percentage*length) ; i++){
		dout(V_INFO | V_NOLABEL) << "#";
	}
	for(int i=0;i< (length-(percentage*length)-1) ;i++){
		dout(V_INFO | V_NOLABEL) << " ";
	}
	time_t now = time(0);
	double timeSpent = difftime(now,starttime);
	int timeTotal = timeSpent/percentage;
	int timeRemaining = timeTotal-timeSpent;
	dout(V_INFO | V_NOLABEL) << "] " << int(percentage*100) << "% (" << part+1 << "/" << full << ") after " << int(timeSpent) << "s (" << timeTotal << "s total, " << timeRemaining << "remaining)" << std::flush; 
}


void fileProgressBar(unsigned nTotal, std::string path){
	int length = 50;
	static unsigned part = 0;
	static const time_t starttime = time(0);
	static const unsigned fillwidthPart = unsigned(1+log10(nTotal));
	static std::ofstream filestream;

	++part;
	const float percentage = float(part) / float(nTotal);

	// write progressbar
	//if(unsigned(percentage*length) > progress || part%10==0){
		if(!filestream.is_open()){
			filestream.open(path.c_str(),std::ofstream::trunc);
		}
		int progress = int(percentage*length);
		filestream << "Progress: [";
		for(int i=0;i<progress-1 ;++i){
			filestream << "#";
		}
    switch(part % 4){
      case 0: filestream << "-"; break;
      case 1: filestream << "\\"; break;
      case 2: filestream << "|"; break;
      case 3: filestream << "/"; break;
    }
		for(int i=0; i < length-progress ; ++i){
			filestream << " ";
		}
		filestream << "] ";

		// write progress in percent
		filestream << std::setfill(' ') << std::setw(3) << unsigned(percentage*100) << "% (" << std::setfill(' ') << std::setw(fillwidthPart) << part << "/" << nTotal << ")" << std::endl;

		// go to next line and write time
		const double timeSpent = difftime(time(0),starttime);
		const int timeTotal = timeSpent/percentage;
		const int timeRemaining = timeTotal-timeSpent;
		filestream << "Runtime " << int(timeSpent) << "s (" << timeTotal << "s total, " << timeRemaining << "s remaining)" << std::flush;

		filestream.close();
	//}
}

