
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
 #  GNU General Public License for more details.
 # 
 #  You should have received a copy of the GNU General Public License
 #  along with HASEonGPU.
 #  If not, see <http://www.gnu.org/licenses/>.
 ##################################################################################


#   function corresponding to the gain calculation
#  this will be more versatile, because with the proper values given, you
#  can call pump and amplification with the same routine
#                    Daniel Albach       2008/03/24

    
import numpy as np

def beta_int(beta_crystal, pulse, const, crystal, steps, int_field, mode, Ntot_gradient):
    # Ensure int_field['s_ems'] is an array
    if not isinstance(int_field['s_ems'], (np.ndarray, list)):
        raise ValueError("int_field['s_ems'] must be an array (list or numpy array)")
    
    # Ensure int_field['s_abs'] is an array
    if not isinstance(int_field['s_abs'], (np.ndarray, list)):
        raise ValueError("int_field['s_abs'] must be an array (list or numpy array)")
    
    # Extracting maximum emission cross section
    int_field = {key: np.array(val) for key, val in int_field.items()}  # Convert lists to numpy arrays if needed
    max_ems_index = np.argmax(int_field['s_ems'])
    int_field['max_ems'] = int_field['s_ems'][max_ems_index]
    int_field['max_abs'] = int_field['s_abs'][max_ems_index]
    
    sigma_abs = int_field['max_abs']  # cm^2
    sigma_ems = int_field['max_ems']  # cm^2
    
    # Discretization
    steps_time = steps['time']
    steps_crystal = steps['crys']
    
    # Extracting constants
    c = const['c']
    h = const['h']
    N_1percent = const['N1per']
    
    # Extracting pump constants
    I_pump = int_field['I']  # W/cm^2
    tau_pump = int_field['T']
    wavelength = int_field['wavelength']  # m
    
    # Extracting crystal constants
    doping = crystal['doping']
    tau_fluo = crystal['tfluo']
    crystal_length = crystal['length']  # cm
    exp_factor = crystal['nlexp']
    
    # Total doping concentration
    Ntot = N_1percent * doping
    
    time_step = tau_pump / (steps_time - 1)
    crystal_step = crystal_length / (steps_crystal - 1)
    
    # Prepare the vectors with zeros
    beta_store = np.zeros((steps_crystal, steps_time))
    pump = np.zeros(steps_crystal)
    pump_l = np.zeros(steps_crystal)
    Ntot_gradient = np.zeros(steps_crystal)
    
    # Exponential gradient
    for igradient in range(steps_crystal):
        Ntot_gradient[igradient] = Ntot * np.exp(np.log(exp_factor) / crystal_length * (igradient - 1) * crystal_step)
    
    for itime in range(steps_time):
        # Initialize pump for the first slice
        pump[0] = I_pump
        
        # Forward direction
        for icrys in range(steps_crystal - 1):
            beta_average = (beta_crystal[icrys] + beta_crystal[icrys + 1]) / 2
            pump[icrys + 1] = pump[icrys] * np.exp(- (sigma_abs - beta_average * (sigma_abs + sigma_ems)) * Ntot_gradient[icrys] * crystal_step)
        
        # Handle Backreflection mode
        if mode['BRM'] == 1:
            beta_crystal = np.flipud(beta_crystal)
            
            pump_BRM = np.zeros_like(pump)
            pump_BRM[0] = pump[-1] * mode['R']
            Ntot_gradient = np.flipud(Ntot_gradient)
            
            # Backward direction
            for jcrys in range(steps_crystal - 1):
                beta_average = (beta_crystal[jcrys] + beta_crystal[jcrys + 1]) / 2
                pump_BRM[jcrys + 1] = pump_BRM[jcrys] * np.exp(- (sigma_abs - beta_average * (sigma_abs + sigma_ems)) * Ntot_gradient[jcrys] * crystal_step)
            
            # Rotate pump_BRM 180 degrees
            pump_BRM = np.rot90(pump_BRM.reshape((steps_crystal, 1)), 2).flatten()
            beta_crystal = np.flipud(beta_crystal)
            
            # Sum pump intensities
            pump_l = pump + pump_BRM
            Ntot_gradient = np.flipud(Ntot_gradient)
        else:
            pump_l = pump.copy()
        
        # Calculate local beta
        for ibeta in range(steps_crystal):
            A1 = sigma_abs * pump_l[ibeta] / (h * c / wavelength)
            C1 = (sigma_abs + sigma_ems) * pump_l[ibeta] / (h * c / wavelength) + 1 / tau_fluo
            beta_crystal[ibeta] = A1 / C1 * (1 - np.exp(-C1 * time_step)) + beta_crystal[ibeta] * np.exp(-C1 * time_step)
        
        pulse[itime] = pump[-1]
        beta_store[:, itime] = beta_crystal.copy()
    
    return beta_crystal, beta_store, pulse, Ntot_gradient

