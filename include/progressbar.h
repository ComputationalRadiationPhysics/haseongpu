#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

//#include <ctime>
#include <string>
void simpleProgressBar(unsigned part, unsigned full);
//void fancyProgressBar(unsigned part, unsigned full, unsigned length, time_t starttime);
void fancyProgressBar(const unsigned nTotal);
void fileProgressBar(const unsigned nTotal, const std::string path);


#endif /* PROGRESS_BAR_H */
