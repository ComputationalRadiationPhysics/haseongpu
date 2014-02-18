/**
 * @brief A collection of functions to print the progress
 *        of a process as a nice progressbar
 *
 * @author Erik Zenker
 * @author Carlchristian Eckert
 * @licence GPLv3
 *
 **/

#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

#include <string>

/**
 * @brief writes the progress of an operation to dout(V_PROGRESS) (see logging.h),
 *        updating the progress every time the function is called.
 *
 * @param part the current progress (must be <= 'full')
 * @param fullthe maximum of the progress (i.e. 100, if you have 100 steps)
 *
 */
void simpleProgressBar(unsigned part, unsigned full);

/**
 * @brief writes the progress of an operation to dout(V_PROGRESS) (see logging.h),
 *        updating the progress every time the function is called.
 *        Works withh a threaded approach or if called directly from a single managing thread
 *
 * @param nTotal the maximum of the progress (i.e. 100, if you have 100 steps)
 *
 */
void fancyProgressBar(const unsigned nTotal);


/**
 * @brief DEPRECATED!
 *        writes the progress of an operation into dout(V_PROGRESS),
 *        updating  the progress every time the function
 *        is called. works with multiple non-threaded callers.
 *
 * @param nTotal the maximum of the progress (i.e. 100, if you have 100 steps)
 * @param path the name of the file to write
 *
 */
void fancyProgressBar(const unsigned current,const unsigned nTotal);


/**
 * @brief writes the progress of an operation into a given file,
 *        updating the file and the progress every time the function
 *        is called.
 *
 * @param nTotal the maximum of the progress (i.e. 100, if you have 100 steps)
 * @param path the name of the file to write
 *
 */
void fileProgressBar(const unsigned nTotal, const std::string path);


#endif /* PROGRESS_BAR_H */
