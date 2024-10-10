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

############################### Laser pump routine ################################
# ASE calulations main file 16kW
# longer calculation without pump action possibel for heat generation
# estimation in the cladding
#
# @author Daniel Albach
# @date   2011-05-03
# @licence: GPLv3
#
###################################################################################


# Import necessary libraries
import numpy as np
import scipy.io as sio
from scipy.interpolate import griddata
import time
import os
from  vtk_wedge import vtk_wedge
from set_variables import set_variables
from calcPhiASE import calcPhiASE


# Constants
c = 3e8  # Speed of light in m/s
h = 6.626e-34  # Planck's constant in Js

# Check if the code is running in Octave
def is_octave():
    try:
        from octave import __version__  # This is a placeholder, replace with actual Octave check if necessary
        return True
    except ImportError:
        return False

# Function to load .vtk files
def vtk_wedge(filename, data, p, t, mesh_z, z_mesh):
    # Implement vtk_wedge functionality based on the original MATLAB code
    pass

# Function to calculate ASE (Amplified Spontaneous Emission)
def calc_phi_ASE(p, t_int, beta_cell, beta_vol, clad_int, clad_number, clad_abs, use_reflections, refractive_indices, reflectivities, normals_x, normals_y, sorted_int, surface, x_center, y_center, normals_p, forbidden, min_rays_per_sample, max_rays_per_sample, mse_threshold, repetitions, n_tot, z_mesh, laser, crystal, mesh_z, device_mode, parallel_mode, max_gpus, n_per_node):
    # Implement ASE calculation functionality based on the original MATLAB code
    pass

# Function to perform beta integration
def beta_int3(beta_crystal, pulse, const, crystal, steps, pump, mode, n_tot_gradient):
    # Implement beta_int3 functionality based on the original MATLAB code
    return beta_crystal, pulse, n_tot_gradient

# Initialization
start_time = time.time()
is_octave_mode = is_octave()

# Crystal parameters
crystal = {
    'doping': 2,
    'length': 0.7,  # cm
    'tfluo': 9.5e-4,  # 1/cm
    'nlexp': 1,
    'levels': 10
}

# Timesteps
steps = {
    'time': 100,
    'crys': crystal['levels']
}

# Pump parameters
pump = {
    'stretch': 1,
    'aspect': 1,
    'diameter': 3,
    's_abs': 0.778e-20,  # Absorption cross section in cm²
    's_ems': 0.195e-20,  # Emission cross section in cm²
    'I': 16e3,  # Pump intensity in W/cm²
    'T': 1e-3,  # Pump duration in s
    'wavelength': 940e-9,  # Pump wavelength in m
    'ry': 3 / 2 * 1,
    'rx': 3 / 2,  # Pump radius
    'exp': 40  # Pump exponent
}

# Laser parameters
laser = {
    's_abs': np.loadtxt('sigma_a.txt'),  # Absorption spectrum cm²
    's_ems': np.loadtxt('sigma_e.txt'),  # Emission spectrum in cm²
    'l_abs': np.loadtxt('lambda_a.txt'),  # Wavelengths absorption spectrum in nm
    'l_ems': np.loadtxt('lambda_e.txt'),  # Wavelengths emission spectrum in nm
    'l_res': 1000,  # Resolution of linear interpolated spectrum
    'I': 1e6,  # Laser intensity
    'T': 1e-8,  # Laser duration
    'wavelength': 1030e-9  # Lasing wavelength in m
}
laser['max_ems'] = np.max(laser['s_ems'])
laser['max_abs'] = laser['s_abs'][np.argmax(laser['s_ems'])]

# Mode parameters
mode = {
    'BRM': 1,  # 1 Backreflection, 0 Transmission mode
    'R': 1,    # Reflectivity of the coating
    'extr': 0  # No extraction
}

# Constants
const = {
    'N1per': 1.38e20,
    'c': c,
    'h': h
}
N_tot = const['N1per'] * crystal['doping']
Ntot_gradient = np.full((crystal['levels'], 1), crystal['doping'] * const['N1per'])
mesh_z = crystal['levels']
z_mesh = crystal['length'] / (mesh_z - 1)
timeslice = 50
timeslice_tot = 150
timetotal = 1e-3  # s
time_t = timetotal / timeslice

# ASE application parameters
ase_params = {
    'maxGPUs': 2,
    'nPerNode': 4,
    'deviceMode': 'gpu',
    'parallelMode': 'threaded',
    'useReflections': True,
    'refractiveIndices': [1.83, 1, 1.83, 1],
    'repetitions': 4,
    'minRaysPerSample': 1e5,
    'maxRaysPerSample': 1e7,
    'mseThreshold': 0.005
}

# Load grid from file
grid_data = sio.loadmat('pt.mat')
p = grid_data['p']
t = grid_data['t']
set_variables(p, t)  # You need to define this function or replace it with equivalent logic
variables_data = sio.loadmat('variable.mat')
set_variables(p, t)  # You need to define this function or replace it with equivalent logic

# Mesh dependent definitions
N_cells = t.shape[0]
beta_cell = np.zeros((p.shape[0], mesh_z))
beta_vol = np.zeros((N_cells, mesh_z - 1))
sorted_int = np.zeros((1, 1))  # Placeholder
reflectivities = np.zeros((1, 2 * sorted_int.shape[1]))

# Initialize variables
phi_ASE = np.zeros((p.shape[0], mesh_z))
dndt_ASE = np.zeros_like(phi_ASE)
flux_clad = np.zeros_like(phi_ASE)
dndt_pump = np.zeros_like(beta_cell)

# Save initial values
vtk_wedge('beta_cell_0.vtk', beta_cell, p, t, mesh_z, z_mesh)
vtk_wedge('dndt_pump_0.vtk', dndt_pump, p, t, mesh_z, z_mesh)
vtk_wedge('dndt_ASE_0.vtk', dndt_ASE, p, t, mesh_z, z_mesh)
sio.savemat('save_0.mat', {'beta_cell': beta_cell, 'dndt_pump': dndt_pump, 'dndt_ASE': dndt_ASE})

temp = pump['T']
time_beta = 1e-6
pump['T'] = time_beta
temp_f = crystal['tfluo']
crystal['tfluo'] = 1
beta_c_2 = np.zeros_like(beta_cell)
intensity = pump['I']

# Perform beta integration for each point
for i_p in range(p.shape[0]):
    beta_crystal = beta_cell[i_p, :]
    pulse = np.zeros(steps['time'])
    pump['I'] = intensity * np.exp(-np.sqrt(p[i_p, 0]**2 / pump['ry']**2 + p[i_p, 1]**2 / pump['rx']**2)**pump['exp'])
    beta_crystal, _, _, Ntot_gradient = beta_int3(beta_crystal, pulse, const, crystal, steps, pump, mode, Ntot_gradient)
    beta_c_2[i_p, :] = beta_crystal

# Update pump rate and crystal fluorescence lifetime
dndt_pump = (beta_c_2 - beta_cell) / time_beta
pump['I'] = intensity
pump['T'] = temp
crystal['tfluo'] = temp_f

# Integrate for the first time
for i_p in range(p.shape[0]):
    for i_z in range(mesh_z):
        beta_cell[i_p, i_z] = crystal['tfluo'] * (dndt_pump[i_p, i_z] - dndt_ASE[i_p, i_z]) * (1 - np.exp(-time_t / crystal['tfluo'])) + beta_cell[i_p, i_z] * np.exp(-time_t / crystal['tfluo'])

vtk_wedge('beta_cell_1.vtk', beta_cell, p, t, mesh_z, z_mesh)

# Cladding definition
clad_number = 1
number_of_triangles, _ = t.shape
clad = np.zeros(number_of_triangles, dtype=int)
clad_abs = 5.5

# Prepare cladding and Yb-doped parts
clad_index = np.where(clad == clad_number)[0]
yb_index = np.where(clad != clad_number)[0]

yb_points_t = np.delete(t, clad_index, axis=0)
yb_points_t2 = np.reshape(yb_points_t, -1)
yb_points_t = np.unique(yb_points_t2)

clad_points_t = np.delete(t, yb_index, axis=0)
clad_points_t2 = np.reshape(clad_points_t, -1)
clad_points_t = np.unique(clad_points_t2)

length_clad = clad_points_t.shape[0]
length_yb = yb_points_t.shape[0]
clad_int = clad.astype(int)

# Main pump loop
for i_slice in range(timeslice_tot - 1):
    print(f'TimeSlice {i_slice} of {timeslice_tot - 1} started')

    # File names
    file_b = f'beta_cell_{i_slice}.vtk'
    file_p = f'dndt_pump_{i_slice}.vtk'
    file_A = f'dndt_ASE_{i_slice}.vtk'
    file_C = f'flux_clad_{i_slice}.vtk'
    
    vtk_wedge(file_b, beta_cell, p, t, mesh_z, z_mesh)
    vtk_wedge(file_p, dndt_pump, p, t, mesh_z, z_mesh)

    # Interpolate beta_vol from beta_cell
    x_1 = p[:, 0]
    y_1 = p[:, 1]
    x_2 = p[:, 0]
    y_2 = p[:, 1]

    x = np.concatenate([x_1, x_2])
    y = np.concatenate([y_1, y_2])

    z_1 = np.zeros_like(x_1)
    z_2 = np.zeros_like(x_2) + z_mesh

    z = np.concatenate([z_1, z_2])

    xi = np.zeros_like(x)  # Replace with actual x_center
    yi = np.zeros_like(y)  # Replace with actual y_center
    zi = np.zeros_like(xi) + z_mesh / 2

    for i_z in range(mesh_z - 1):
        v_1 = beta_cell[:, i_z]
        v_2 = beta_cell[:, i_z + 1]

        v = np.concatenate([v_1, v_2])

        beta_vol[:, i_z] = griddata((x, y, z), v, (xi, yi, zi), method='linear')

        z += z_mesh
        zi += z_mesh

    # Call external ASE application
    phi_ASE, mse_values, N_rays = calc_phi_ase(
        p, t, beta_cell, beta_vol, clad_int, clad_number, clad_abs, 
        ase_params['useReflections'], ase_params['refractiveIndices'], 
        ase_params['reflectivities'], None, None, None, None, None, 
        ase_params['minRaysPerSample'], ase_params['maxRaysPerSample'], 
        ase_params['mseThreshold'], ase_params['repetitions'], N_tot, 
        z_mesh, laser, crystal, mesh_z, ase_params['deviceMode'], 
        ase_params['parallelMode'], ase_params['maxGPUs'], ase_params['nPerNode']
    )

    # Calculate dn/dt ASE for Yb:YAG points
    for i_p in range(length_yb):
        for i_z in range(mesh_z):
            pos_yb = yb_points_t[i_p]
            g_l = -(N_tot * laser['max_abs'] - N_tot * beta_cell[pos_yb, i_z] * (laser['max_ems'] + laser['max_abs']))
            dndt_ASE[pos_yb, i_z] = g_l * phi_ASE[pos_yb, i_z] / crystal['tfluo']

    # Calculate flux for cladding points
    for i_p in range(length_clad):
        for i_z in range(mesh_z):
            pos_clad = clad_points_t[i_p]
            flux_clad[pos_clad, i_z] = clad_abs * phi_ASE[pos_clad, i_z] / crystal['tfluo']

    vtk_wedge(file_A, dndt_ASE, p, t, mesh_z, z_mesh)
    vtk_wedge(file_C, flux_clad, p, t, mesh_z, z_mesh)
    sio.savemat(f'save_{i_slice}.mat', {'beta_cell': beta_cell, 'dndt_pump': dndt_pump, 'dndt_ASE': dndt_ASE, 'flux_clad': flux_clad})

    # Prepare next timeslice
    temp = pump['T']
    pump['T'] = time_beta
    temp_f = crystal['tfluo']
    crystal['tfluo'] = 1
    beta_c_2 = np.zeros_like(beta_cell)
    intensity = pump['I']

    for i_p in range(p.shape[0]):
        beta_crystal = beta_cell[i_p, :]
        pulse = np.zeros(steps['time'])
        pump['I'] = intensity * np.exp(-np.sqrt(p[i_p, 0]**2 / pump['ry']**2 + p[i_p, 1]**2 / pump['rx']**2)**pump['exp'])
        beta_crystal, _, _, Ntot_gradient = beta_int3(beta_crystal, pulse, const, crystal, steps, pump, mode, Ntot_gradient)
        beta_c_2[i_p, :] = beta_crystal

    if i_slice < timeslice_tot - 1:
        dndt_pump = (beta_c_2 - beta_cell) / time_beta
    else:
        dndt_pump = np.zeros_like(beta_cell)

    pump['I'] = intensity
    pump['T'] = temp
    crystal['tfluo'] = temp_f

    # Recalculate beta_cell
    for i_p in range(p.shape[0]):
        for i_z in range(mesh_z):
            beta_cell[i_p, i_z] = crystal['tfluo'] * (dndt_pump[i_p, i_z] - dndt_ASE[i_p, i_z]) * (1 - np.exp(-time_t / crystal['tfluo'])) + beta_cell[i_p, i_z] * np.exp(-time_t / crystal['tfluo'])

# Final output
file_b = f'beta_cell_{timeslice_tot}.vtk'
file_p = f'dndt_pump_{timeslice_tot}.vtk'
file_A = f'dndt_ASE_{timeslice_tot}.vtk'
file_C = f'flux_clad_{timeslice_tot}.vtk'

vtk_wedge(file_b, beta_cell, p, t, mesh_z, z_mesh)
vtk_wedge(file_p, dndt_pump, p, t, mesh_z, z_mesh)
vtk_wedge(file_A, dndt_ASE, p, t, mesh_z, z_mesh)
vtk_wedge(file_C, flux_clad, p, t, mesh_z, z_mesh)

print('Calculations finished')
end_time = time.time()
print(f'Elapsed time: {end_time - start_time} seconds')
