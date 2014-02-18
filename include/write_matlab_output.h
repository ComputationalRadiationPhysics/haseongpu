#ifndef WRITE_DNDT_ASE_H
#define WRITE_DNDT_ASE_H

#include <vector>
#include <string>

/**
 * @brief creates textfiles containing all the results as a 2D matrix
 *        (needs to be reshaped in matlab, see calcPhiASE.m)
 *
 * @param experimentPath the path where to create the output files
 * @param ase the phi_ASE values to print
 * @param N_rays the number of rays used per sample point
 * @param mse_values the error for each sample point
 * @param numberOfSamples the amount of vertices that were sampled in the mesh (over all levels)
 * @param numberOfLevels the amount of levels in which the gain medium was split 
 *
 * @author Carlchristian Eckert
 * @author Erik Zenker
 * @author Marius Melzer
 *
 * @license GPLv3
 */
void writeMatlabOutput(
    const std::string experimentPath,
    const std::vector<float> ase,
    const std::vector<unsigned> N_rays, 
    const std::vector<double> mse_values,
    const unsigned numberOfSamples,
    const unsigned numberOfLevels
    );

#endif 
