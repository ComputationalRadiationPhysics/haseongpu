#include "progressbar.h"
#include <stdio.h>
#include <unistd.h>
#include <ctime>
#include <logging.h>

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
