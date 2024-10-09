####################################################################################
# Copyright 2014 Erik Zenker, Carlchristian Eckert
#
# This file is part of HASEonGPU
#
# HASEonGPU is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HASEonGPU is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HASEonGPU.
# If not, see <http://www.gnu.org/licenses/>.
###################################################################################

import numpy as np
import os
import warnings
import shutil 

############################## clean_IO_files #################################
#
# deletes the temporary folder and the dndt_ASE.txt
# 
# @param TMP_FOLDER the folder to remove
#
#
###############################################################################
def clean_IO_files (TMP_FOLDER):

  if os.path.exists(TMP_FOLDER) and os.path.isdir(TMP_FOLDER):
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      shutil.rmtree(TMP_FOLDER)


# ################## parse_calcPhiASE_output #############################
#
# takes all the variables and puts them into textfiles
# so the CUDA code can parse them. Take care that the
# names of the textfiles match those in the parsing function
# 
# for most parameters, see calcPhiASE (above)
# @param FOLDER the folder in which to create the input files (usually a
#               temporary folder visible by all the participating nodes)
#
#
##########################################################################

def create_calcPhiASE_input (points,
                            triangleNormalsX,  
                            triangleNormalsY,   
                            forbiddenEdge,  
                            triangleNormalPoint,    
                            triangleNeighbors,  
                            trianglePointIndices,   
                            thickness,  
                            numberOfLevels, 
                            nTot,
                            betaVolume, 
                            laserParameter, 
                            crystal,    
                            betaCells,  
                            triangleSurfaces,   
                            triangleCenterX,    
                            triangleCenterY,    
                            claddingCellTypes,  
                            claddingNumber, 
                            claddingAbsorption, 
                            refractiveIndices,  
                            reflectivities, 
                            FOLDER):
    
    CURRENT_DIR = os.getcwd()

    os.mkdir(FOLDER)
    os.chdir(FOLDER) 

    points=                   np.transpose(points) 
    triangleNormalsX=         np.transpose(triangleNormalsX)
    triangleNormalsY=         np.transpose(triangleNormalsY)
    forbiddenEdge=            np.transpose(forbiddenEdge)
    triangleNormalPoint=      np.transpose(triangleNormalPoint)
    triangleNeighbors=        np.transpose(triangleNeighbors)
    trianglePointIndices=     np.transpose(trianglePointIndices)
    betaVolume=               np.transpose(betaVolume)
    laserParameter['s_abs']=  np.transpose(laserParameter['s_abs'])
    laserParameter['s_ems']=  np.transpose(laserParameter['s_ems'])
    laserParameter['l_abs']=  np.transpose(laserParameter['l_abs'])
    laserParameter['l_ems']=  np.transpose(laserParameter['l_ems'])
    betaCells=                np.transpose(betaCells)
    triangleSurfaces=         np.transpose(triangleSurfaces)
    triangleCenterX=          np.transpose(triangleCenterX)
    triangleCenterY=          np.transpose(triangleCenterY)
    claddingCellTypes=        np.transpose(claddingCellTypes)
    refractiveIndices=        np.transpose(refractiveIndices)
    reflectivities=           np.transpose(reflectivities)

    #with open('points.txt','w') as f:
    #  f.write('\n'.join(map(str, points[:,0]))+'\n'+'\n'.join(map(str, points[:,1])))
    np.savetxt('points.txt', points, delimiter='\n', fmt='%.50f')

    #with open('triangleNormalsX.txt','w') as f:
    np.savetxt('triangleNormalsX.txt', triangleNormalsX, delimiter='\n', fmt='%.50f')

    #with open('triangleNormalsY.txt','w') as f:
    #  f.write('\n'.join(map(str, triangleNormalsY)))
    np.savetxt('triangleNormalsY.txt', triangleNormalsY, delimiter='\n', fmt='%.50f')

    #with open('forbiddenEdge.txt','w') as f:
    #  f.write('\n'.join(map(str, forbiddenEdge)))
    np.savetxt('forbiddenEdge.txt', forbiddenEdge, delimiter='\n', fmt='%d')

    #with open('triangleNormalPoint.txt','w') as f:
    #  f.write('\n'.join(map(str, triangleNormalPoint)))
    np.savetxt('triangleNormalPoint.txt', triangleNormalPoint, delimiter='\n', fmt='%d')

    #with open('triangleNeighbors.txt','w') as f:
    #  f.write('\n'.join(map(str, triangleNeighbors)))
    np.savetxt('triangleNeighbors.txt', triangleNeighbors, delimiter='\n', fmt='%d')

    #with open('trianglePointIndices.txt','w') as f:
    #  f.write('\n'.join(map(str, trianglePointIndices)))
    np.savetxt('trianglePointIndices.txt', trianglePointIndices, delimiter='\n', fmt='%d')

    # thickness of one slice!
    with open('thickness.txt','w') as f:
      f.write(str(thickness) + '\n')

    # number of slices
    with open('numberOfLevels.txt','w') as f:
      f.write(str(numberOfLevels) + '\n')
 
    with open('numberOfTriangles.txt','w') as f:
      a = len(trianglePointIndices[0])
      f.write(str(a) + '\n')

    with open('numberOfPoints.txt','w') as f:
      a = len(points[0])
      f.write(str(a) + '\n')

    with open('nTot.txt','w') as f:
      nTot= float(nTot)
      f.write(str(nTot) + '\n')

    #with open('betaVolume.txt','w') as f:
    #  f.write('\n'.join(map(str, betaVolume)))
    np.savetxt('betaVolume.txt', betaVolume, delimiter='\n', fmt='%.50f')

    #with open('sigmaA.txt','w') as f:
    #  f.write('\n'.join(map(str, laserParameter['s_abs'])))
    np.savetxt('sigmaA.txt', laserParameter['s_abs'], delimiter='\n', fmt='%.50f')

    #with open('sigmaE.txt','w') as f:
    #  f.write('\n'.join(map(str, laserParameter['s_ems'])))
    np.savetxt('sigmaE.txt', laserParameter['s_ems'], delimiter='\n', fmt='%.50f')
    
    #with open('lambdaA.txt','w') as f:
    #  f.write('\n'.join(map(str, laserParameter['l_abs'])))
    np.savetxt('lambdaA.txt', laserParameter['l_abs'], delimiter='\n', fmt='%.50f')

    #with open('lambdaE.txt','w') as f:
    #  f.write('\n'.join(map(str, laserParameter['l_ems'])))
    np.savetxt('lambdaE.txt', laserParameter['l_ems'], delimiter='\n', fmt='%.50f')

    with open('crystalTFluo.txt','w') as f:
      f.write(str(crystal['tfluo']) + '\n')
    
    #with open('betaCells.txt','w') as f:
    #  f.write('\n'.join(map(str, betaCells)))
    np.savetxt('betaCells.txt', betaCells, delimiter='\n', fmt='%.50f')

    #with open('triangleSurfaces.txt','w') as f:
    #  f.write('\n'.join(map(str, triangleSurfaces)))
    np.savetxt('triangleSurfaces.txt', triangleSurfaces, delimiter='\n', fmt='%.50f')  

    #with open('triangleCenterX.txt','w') as f:
    #  f.write('\n'.join(map(str, triangleCenterX)))
    np.savetxt('triangleCenterX.txt', triangleCenterX, delimiter='\n', fmt='%.50f')

    #with open('triangleCenterY.txt','w') as f:
    #  f.write('\n'.join(map(str, triangleCenterY)))
    np.savetxt('triangleCenterY.txt', triangleCenterY, delimiter='\n', fmt='%.50f')

    #with open('claddingCellTypes.txt','w') as f:
    #  f.write('\n'.join(map(str, claddingCellTypes)))
    np.savetxt('claddingCellTypes.txt', claddingCellTypes, delimiter='\n', fmt='%d')

    with open('claddingNumber.txt','w') as f:
      f.write(str(claddingNumber) + '\n')

    with open('claddingAbsorption.txt','w') as f:
      f.write(str(claddingAbsorption) + '\n')

    #with open('refractiveIndices.txt','w') as f:
    #  f.write('\n'.join(map(str, refractiveIndices)))
    np.savetxt('refractiveIndices.txt', refractiveIndices, delimiter='\n', fmt='%3.5f')

    #with open('reflectivities.txt','w') as f:
    #  f.write('\n'.join(map(str, reflectivities)))
    np.savetxt('reflectivities.txt', reflectivities, delimiter='\n', fmt='%.50f')
    
    os.chdir(CURRENT_DIR)

    triangleNormalsX=         np.transpose(triangleNormalsX)
    triangleNormalsY=         np.transpose(triangleNormalsY)
    forbiddenEdge=            np.transpose(forbiddenEdge)
    triangleNormalPoint=      np.transpose(triangleNormalPoint)
    triangleNeighbors=        np.transpose(triangleNeighbors)
    trianglePointIndices=     np.transpose(trianglePointIndices)
    betaVolume=               np.transpose(betaVolume)
    laserParameter['s_abs']=  np.transpose(laserParameter['s_abs'])
    laserParameter['s_ems']=  np.transpose(laserParameter['s_ems'])
    laserParameter['l_abs']=  np.transpose(laserParameter['l_abs'])
    laserParameter['l_ems']=  np.transpose(laserParameter['l_ems'])
    betaCells=                np.transpose(betaCells)
    triangleSurfaces=         np.transpose(triangleSurfaces)
    triangleCenterX=          np.transpose(triangleCenterX)
    triangleCenterY=          np.transpose(triangleCenterY)
    claddingCellTypes=        np.transpose(claddingCellTypes)
    refractiveIndices=        np.transpose(refractiveIndices)
    reflectivities=           np.transpose(reflectivities)

######################### parse_calcPhiASE_output #############################
#
# takes the output from the CUDA code and fills it into a variable
# assumes that the matrix is saved as a 3D-matrix where the first line 
# denotes the dimensions and the second line is the whole data
#
# @param FOLDER the folder which contains the output files
#
###############################################################################

def parse_calcPhiASE_output (FOLDER):
  CURRENT_DIR = os.getcwd()
  os.chdir(FOLDER)
  with open('phi_ASE.txt','r') as fid:
    arraySize = [int(x) for x in next(fid).split()]
    phiASE = [float(x) for x in next(fid).split()] 
    phiASE = np.reshape(phiASE,arraySize, order= 'F') #order= 'F' beacause file was adapted from Matlab

  with open('mse_values.txt', 'r') as fid: 
    arraySize= [int(x) for x in next(fid).split()]
    mseValues = [float(x) for x in next(fid).split()]
    mseValues = np.reshape(mseValues,arraySize, order= 'F')


  with open('N_rays.txt', 'r') as fid: 
    arraySize = [int(x) for x in next(fid).split()]
    raysUsedPerSample = [float(x) for x in next(fid).split()]
    raysUsedPerSample = np.reshape(raysUsedPerSample,arraySize, order= 'F')
  
  os.chdir(CURRENT_DIR)
  return([mseValues, raysUsedPerSample, phiASE])

################################ calcPhiASE ########################################
# % calculates the phiASE values for a given input
# most meshing paramers are given through the function parameters.
# However, many parameters for optimization of the computation are
# set in the beginning of the function (adjust as needed)
# 
# for most mesh parameters see README file
#
# @return phiASE the ASE-Flux for all the given sample points
# @return mseValues the MeanSquaredError values corresponding to phiASE
# @return raysUsedPerSample the number of rays used to calculate each phiASE value
#
####################################################################################

def calcPhiASE(
  points,
  trianglePointIndices,
  betaCells,
  betaVolume,
  claddingCellTypes,
  claddingNumber,
  claddingAbsorption,
  useReflections,
  refractiveIndices,
  reflectivities,
  triangleNormalsX,
  triangleNormalsY,
  triangleNeighbors,
  triangleSurfaces,
  triangleCenterX,
  triangleCenterY,
  triangleNormalPoint,
  forbiddenEdge,
  minRaysPerSample,
  maxRaysPerSample,
  mseThreshold,
  repetitions,
  nTot,
  thickness,
  laserParameter,
  crystal,
  numberOfLevels,
  deviceMode,
  parallelMode,
  maxGPUs,
  nPerNode
  ):

  ############# auto-generating some more input #############
  minSample=0
  nP = len(points)
  maxSample=(numberOfLevels*nP)-1

  if useReflections == True:
      REFLECT=' --reflection=1'
  else:
      REFLECT=' --reflection=0'

  if parallelMode=='mpi':
    Prefix=[ 'mpiexec -npernode ', str(nPerNode), ' ' ]
    maxGPUs=1
  else:
    Prefix=''

  # create a tmp-folder in the same directory as this script
  # before filedirectory, name and extension, but only filddirectory is used
  CALCPHIASE_DIR= os.getcwd() 
  TMP_FOLDER= os.path.join(CALCPHIASE_DIR,'input_tmp')

  ################## doing the computation ##################
  # make sure that the temporary folder is clean 
  clean_IO_files(TMP_FOLDER)

  # create new input in the temporary folder
  create_calcPhiASE_input(points,
                          triangleNormalsX,
                          triangleNormalsY,
                          forbiddenEdge,
                          triangleNormalPoint,
                          triangleNeighbors,
                          trianglePointIndices,
                          thickness,
                          numberOfLevels,
                          nTot,
                          betaVolume,
                          laserParameter,
                          crystal,
                          betaCells,
                          triangleSurfaces,
                          triangleCenterX,
                          triangleCenterY,
                          claddingCellTypes,
                          claddingNumber,
                          claddingAbsorption,
                          refractiveIndices,
                          reflectivities,
                          TMP_FOLDER)

  # do the propagation
  status = os.system(Prefix +CALCPHIASE_DIR+'/../../build/calcPhiASE '
                           + '--parallel-mode='+ parallelMode
                           + ' --device-mode='+ deviceMode
                           + ' --min-rays='+ str(int(minRaysPerSample))
                           + ' --max-rays='+ str(int(maxRaysPerSample))
                           + REFLECT
                           + ' --input-path='+ TMP_FOLDER
                           + ' --output-path='+ TMP_FOLDER
                           + ' --min-sample-i='+ str(minSample)
                           + ' --max-sample-i='+ str(maxSample)
                           + ' --ngpus='+ str(maxGPUs)
                           + ' --repetitions='+ str(repetitions)
                           + ' --mse-threshold='+ str(mseThreshold)
                           + ' --spectral-resolution='+ str(laserParameter['l_res']) )

  if status != 0:
      print('this step of the raytracing computation did NOT finish successfully. Aborting.')
      exit()


  # get the result
  [ mseValues, raysUsedPerSample, phiASE ] = parse_calcPhiASE_output(TMP_FOLDER)

  # cleanup
  clean_IO_files(TMP_FOLDER)


