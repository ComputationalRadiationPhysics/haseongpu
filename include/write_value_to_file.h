#ifndef WRITE_ARRAY_FILE_H
#define WRITE_ARRAY_FILE_H

#include <string>

template <typename T>
int writeValueToFile(
	       const T value,
	       const std::string filename);

#endif /* WRITE_ARRAY_FILE_H*/
