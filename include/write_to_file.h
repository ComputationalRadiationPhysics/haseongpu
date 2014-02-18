/**
 * @author Carlchristan Eckert
 * @author Erik Zenker 
 * @license GPLv3
 */

#pragma once
#include <string>
#include <vector>

/** 
 * @brief writes a value to a file, where the filename is appended with 2 longish indices
 *
 * @param value the value to write to a file
 * @param path the path where to create the file
 * @param indexName1 identifier of the first index
 * @param index1 the value of the first index
 * @param indexName2 identifier of the second index
 * @param index2 the value of the second index
 */
int writeValueToFile(
    const float value, 
    const std::string path, 
    const std::string indexName1, 
    const int index1, 
    const std::string indexName2, 
    const int index2
    );


/**
 * @brief writes a vector to a file
 *
 * @param v the vector to write
 * @param pFilename the name of the output file
 *
 */
void writeVectorToFile(std::vector<double> v, std::string pFilename);
