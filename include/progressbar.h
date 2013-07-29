#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

#include <ctime>
void simpleProgressBar(unsigned part, unsigned full);
void fancyProgressBar(unsigned part, unsigned full, unsigned length, time_t starttime);


#endif /* PROGRESS_BAR_H */
