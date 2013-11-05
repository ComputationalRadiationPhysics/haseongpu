#ifndef WRITE_ARRAY_FILE_H
#define WRITE_ARRAY_FILE_H

#include <string>

//template < typename T >
int writeValueToFile(
    const float value, 
    const std::string path, 
    const std::string indexName1, 
    const int index1, 
    const std::string indexName2, 
    const int index2
    );

#endif /* WRITE_ARRAY_FILE_H*/
