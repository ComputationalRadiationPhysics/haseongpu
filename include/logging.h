#ifndef OCTRACE_LOGGING_H
#define OCTRACE_LOGGING_H

#include <iostream>

#define V_QUIET 0
#define V_ERROR 1
#define V_WARNING 2
#define V_INFO 4
#define V_STAT 8
#define V_DEBUG 16

#define COLOR_ERROR 31
#define COLOR_WARN 31
#define COLOR_INFO 32
#define COLOR_STATISTIC 32
#define COLOR_DEBUG 36

extern unsigned verbosity;

struct nullstream : std::ostream {
  nullstream() : std::ostream(0) { }
};

std::ostream& dout(unsigned activation_level);

#endif
