#ifndef WRITE_DNDT_ASE_H
#define WRITE_DNDT_ASE_H

#include <vector>

/**
 * @brief creates a textfile containing all the dnd_ase values
 *
 * @param *ase the values to print
 *
 * @author Carlchristian Eckert
 * @author Erik Zenker
 * @author Marius Melzer
 *
 * @license GPLv3
 */
void writeDndtAse(std::vector<float>* ase);

#endif 
