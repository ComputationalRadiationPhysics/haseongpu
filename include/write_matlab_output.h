#ifndef WRITE_DNDT_ASE_H
#define WRITE_DNDT_ASE_H

#include <vector>
#include <string>

/**
 * @brief creates textfiles containing all the results (one column per wavelength!)
 *
 * @param *ase the phi_ASE values to print
 * @param *N_rays the number of rays used per wavelength
 * @param *expectedValues the error for each sample point
 * @param numberOfWavelengths how many wavelengths were used
 * @param numberOfSamples the amount of vertices that were sampled in the mesh
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
    const std::vector<double> expectedValues,
    const unsigned numberOfSamples,
    const unsigned numberOfLevels
    );

#endif 
