 ##################################################################################
 # Copyright 2013 Daniel Albach, Erik Zenker, Carlchristian Eckert
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
 #################################################################################

#  function corresponding to the gain calculation
# this will be more versatile, because with the proper values given, you
# can call pump and amplification with the same routine
#                   Daniel Albach       2009/07/02

# calling has to provide the beta distribution before interaction (at least
# zero-array), the temporal pulse-shape in terms of intensity
# (mode.extr==1) or a given (constant) pump-intensity (mode.extr==0), all
# the necessary constants also for the crystal and steps (temporal and
# spatial), informations about the intensity filed (pump or pulse) and the
# doping-gradient (given externally)

import numpy as np

def beta_int3(beta_crystal, pulse, const, crystal, steps, int_field, mode, Ntot_gradient):
    
    # declarations
    int_field['max_ems'] = np.max(int_field['s_ems']) 
    int_field['max_abs'] = np.max(int_field['s_abs'])
    sigma_abs = int_field['max_abs'] # cm^2
    sigma_ems = int_field['max_ems'] # cm^2
    
    # discretization
    steps_time = steps['time']
    steps_crystal = steps['crys']
    
    # extracting the constants
    c = const['c']
    h = const['h']
    N_1percent = const['N1per']
    
    # extracting the "pump" constants
    I_pump = int_field['I'] # W/cm²
    tau_pump = int_field['T']
    wavelength = int_field['wavelength'] # m
    
    # extracting the crystal constants
    doping = crystal['doping']
    tau_fluo = crystal['tfluo']
    crystal_length = crystal['length'] # cm
    exp_factor = crystal['nlexp']
    
    # total doping concentration
    Ntot = N_1percent * doping
    
    time_step = tau_pump / (steps_time - 1)
    crystal_step = crystal_length / (steps_crystal - 1)
    
    # prepare the vectors with zeros
    beta_store = np.zeros((steps_crystal, steps_time))
    pump_l = np.zeros((steps_crystal, 1))
    pump_BRM = np.zeros((steps_crystal, 1))
    pump = np.zeros((steps_crystal, 1))
    
    # Ntot_gradient = np.zeros((steps_crystal, 1))
    
    # exponential gradient
    # if mode['dop'] == 0:
    #     for igradient in range(steps_crystal):
    #         Ntot_gradient[igradient] = Ntot * np.exp(np.log(exp_factor) / crystal_length * (igradient - 1) * crystal_step)
    # else:
    #     doping gradient given from the outside
    
    for itime in range(steps_time):
        # now the first slice moves
        # go with it into the slices of the crystal
        # for the first slice it is always I_pump, the second gets the
        # estimation with an average beta from m -> m + 1 with an exponential
        # function and so on - boundaries!
        if mode['extr'] == 0:
            pump[0] = I_pump
        else:
            pump[0] = pulse[itime]
        
        # this is the positive direction
        for icrys in range(steps_crystal - 1):
            # step one is from point one to two for I_pump
            beta_average = (beta_crystal[icrys] + beta_crystal[icrys + 1]) / 2
            pump[icrys + 1] = pump[icrys] * np.exp(-(sigma_abs - beta_average * (sigma_abs + sigma_ems)) * Ntot_gradient[icrys] * crystal_step)
        
        # now make the case of Backreflection - rough approximation, that the
        # beta hasn't changed during the roundtrip - valid for the pump
        # (integration step is ~5 orders of magnitude longer than the roundtrip),
        # but for the pulse it might get a pity - solution:

        if mode['BRM']== 1:
            beta_crystal = np.flipud(beta_crystal)
            
            pump_BRM[0] = pump[-1]*mode['R']
            Ntot_gradient = np.flipud(Ntot_gradient)
        
    #   this is the negative direction
            for jcrys in range(steps_crystal-1): 
    #           step one is from point one to two for I_pump
                beta_average = (beta_crystal[jcrys]+beta_crystal[jcrys+1])/2
                pump_BRM[jcrys+1] = pump_BRM[jcrys] * np.exp(-(sigma_abs - beta_average*(sigma_abs+sigma_ems))*Ntot_gradient[jcrys]*crystal_step)
    #         now turn the second pumppart and the beta again
            pump_BRM= np.rot90(pump_BRM,2)
            beta_crystal = np.flipud(beta_crystal)
            
    #         full pump intensity is I+ + I-
            Ntot_gradient = np.flipud(Ntot_gradient)
    #         in the case of BRM the return value has to be I-!
            pulse[itime]=pump_BRM[0]
        else:
            pulse[itime] = pump_l[icrys+1]
            
        pump_l = pump + pump_BRM

        #   now calculate the local beta
        for ibeta in range(steps_crystal):
            A1 = sigma_abs*pump_l[ibeta]/(h*c/wavelength)
            C1 = (sigma_abs+sigma_ems)*pump_l[ibeta]/(h*c/wavelength)+1/tau_fluo
            
            beta_crystal[ibeta] = A1/C1*(1-np.exp(-C1*time_step))+ beta_crystal[ibeta]*np.exp(-C1*time_step)
        
        #     if icrys or jcrys makes no difference
        beta_store[:,itime]=beta_crystal

    return([beta_crystal,beta_store,pulse,Ntot_gradient])