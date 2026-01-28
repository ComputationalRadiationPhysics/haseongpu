# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Copyright 2013 Daniel Albach, Erik Zenker, Carlchristian Eckert
#
#  This file is part of HASEonGPU
#
#  HASEonGPU is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  HASEonGPU is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with HASEonGPU.
#  If not, see <http://www.gnu.org/licenses/>.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
import scipy.io
import time
from scipy.interpolate import griddata

# Start timer
start_time = time.time()

timeslices_tot = 150

# Main loop over time slices
for i_s in range(1, timeslices_tot + 1):
    print(f"TimeSlice {i_s}")

    i_ss = i_s - 1
    filename = f'save_{i_ss}.mat'

    # Load MATLAB file
    mat_data = scipy.io.loadmat(filename)
    
    # Extract relevant variables from mat_data
    pump = mat_data['pump']
    steps = mat_data['steps']
    crystal = mat_data['crystal']
    p = mat_data['p']
    beta_cell = mat_data['beta_cell']
    laser = mat_data['laser']
    
    time_step = pump['T'][0][0] / (steps['time'][0][0] - 1)
    crystal_step = crystal['length'][0][0] / (steps['crys'][0][0] - 1)

    gainmap = np.zeros((p.shape[0], 1))

    for i_p in range(p.shape[0]):

        beta_crystal = np.zeros_like(beta_cell[i_p, :])
        beta_crystal[0, :] = beta_cell[i_p, :]

        laser['I'] = 1

        # Grids for plotting the lines
        grid_t = np.linspace(0, pump['T'][0][0], steps['time'][0][0])
        grid_z = np.linspace(0, crystal['length'][0][0], steps['crys'][0][0])

        gain_local = np.zeros((steps['crys'][0][0], 1))

        pulse = np.ones((steps['time'][0][0], 1)) * laser['I']

        # Extraction mode
        mode = {'BRM': 1, 'R': 1, 'extr': 1}

        # Laser wavelength
        laser['max_ems'] = np.max(laser['s_ems'])
        i = np.argmax(laser['s_ems'])
        laser['max_abs'] = laser['s_abs'][i]

        sigma_abs_L = laser['max_abs']
        sigma_ems_L = laser['max_ems']
        beta_min = np.full((25,), sigma_abs_L / (sigma_abs_L + sigma_ems_L))
        grid_z_beta_min = np.linspace(0, crystal['length'][0][0], 25)

        # Time constants for the pump
        time_step_ex = laser['T'][0][0] / (steps['time'][0][0] - 1)
        grid_t_ex = np.linspace(0, laser['T'][0][0], steps['time'][0][0])

        energy_pulse = np.trapezoid(pulse, grid_t_ex)

        # Energy dump for the plot
        energy_pulse_round = np.zeros((1, 1))
        energy_pulse_round[0, 0] = energy_pulse

        # Call to the function for extraction
        # Assuming beta_int is a function defined elsewhere, replace with its equivalent Python implementation
        beta_crystal, beta_store, pulse = beta_int(beta_crystal, pulse, const, crystal, steps, laser, mode)

        gainmap[i_p, 0] = pulse[0] / laser['I']

    # Create meshgrid for interpolation
    x_grid, y_grid = np.meshgrid(np.arange(-3, 3.5, 0.5), np.arange(-3, 3.5, 0.5))
    gainmap_Interp = griddata(p[:, :2], gainmap[:, 0], (x_grid, y_grid), method='cubic')

# Extract data
gain_line = np.zeros((timeslices_tot, 1))
for i_t in range(timeslices_tot):
    gain_line[i_t, 0] = gainmap_Interp[7, 7, i_t]

# Save gain_line to a text file
with open('gain_line.txt', 'w') as f:
    np.savetxt(f, gain_line, fmt='%.50f')

# End timer
end_time = time.time()
print(f"Elapsed time: {end_time - start_time} seconds")
