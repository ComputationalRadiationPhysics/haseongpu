#include "progressbar.h"
#include <stdio.h>
#include <unistd.h>
#include <ctime>

void simpleProgressBar(unsigned part, unsigned full){
	unsigned length = 80;

	float percentage = (float(part)+1) / float(full);

	printf("\r Progress: [");
	for(int i=0 ; i < (percentage*length) ; i++){
		printf("#");
	}
	for(int i=0;i< length-(percentage*length) ;i++){
		printf(" ");
	}
	printf("] %d%% (%d/%d)",int(percentage*100),part+1,full);
	fflush(stdout);
}

void fancyProgressBar(unsigned part, unsigned full, unsigned length, time_t starttime){

	float percentage = (float(part)+1) / float(full);

	printf("\r Progress: [");
	for(int i=0 ; i < (percentage*length) ; i++){
		printf("#");
	}
	for(int i=0;i< (length-(percentage*length)-1) ;i++){
		printf(" ");
	}
	time_t now = time(0);
	double timeSpent = difftime(now,starttime);
	int timeTotal = timeSpent/percentage;
	int timeRemaining = timeTotal-timeSpent;
	printf("] %.1f%% (%d/%d) after %ds (%ds total, %ds remaining)",percentage*100,part+1,full,int(timeSpent),timeTotal,timeRemaining);
	fflush(stdout);
}
