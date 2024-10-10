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
 # but WITHOUT ANY WARRANTY without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with HASEonGPU.
 # If not, see <http://www.gnu.org/licenses/>.
 #################################################################################

# This will calculate the gain distribution within an Yb:YAG-Crystal
# Monochromatic approach                Daniel Albach
#                   last edit:  2024/10/09







import numpy as np
import matplotlib.pyplot as plt
from beta_int import beta_int  # Importing beta_int from the same directory


# Start timer
import time
start_time = time.time()

# Constants
const = {
    'N1per': 1.38e20,
    'c': 3e8,
    'h': 6.626e-34
}

# Crystal parameters
crystal = {
    'doping': 2,
    'length': 0.7,  # [cm]
    'tfluo': 9.5e-4,
    'nlexp': 1
}

# Steps
steps = {
    'time': 1000,
    'crys': 300
}

# Pump parameters
pump = {
    's_abs': np.array([0.76e-20]),  # absorption cross section in cm²
    's_ems': np.array([0.220e-20]),  # emission cross section in cm²
    'I': 10e3,  # Pump intensity in W/cm²
    'T': 1e-3,  # pump duration in seconds
    'wavelength': 940e-9  # pump wavelength in meters
}

# Laser parameters
laser = {
    's_abs': np.array([1.10e-21]),  # absorption cross section in cm²
    's_ems': np.array([2.40e-20]),  # emission cross section in cm²
    'I': 1e6,  # laser intensity in W/cm²
    'T': 1e-8,  # laser duration
    'wavelength': 1030e-9  # lasing wavelength in meters
}

# Mode definitions
mode = {
    'BRM': 1,  # 1 Backreflection, 0 Transmission mode
    'R': 1  # reflectivity of the coating
}

# Prepare beta_crystal as a vector
beta_crystal = np.zeros(steps['crys'])
pulse = np.zeros(steps['time'])

time_step = pump['T'] / (steps['time'] - 1)
crystal_step = crystal['length'] / (steps['crys'] - 1)

gain_local = np.zeros(steps['crys'])

# Call beta_int for the pump
beta_crystal, beta_store, pulse, Ntot_gradient = beta_int(
    beta_crystal, pulse, const, crystal, steps, pump, mode, Ntot_gradient=None)

# Grids for plotting
grid_t = np.linspace(0, pump['T'], steps['time'])
grid_z = np.linspace(0, crystal['length'], steps['crys'])

# Laser wavelength parameters
sigma_abs_L = 1.16e-21
sigma_ems_L = 2.48e-20
beta_min = np.ones(25) * sigma_abs_L / (sigma_abs_L + sigma_ems_L)
grid_z_beta_min = np.linspace(0, crystal['length'], 25)

# Energy density calculation
energy_density = {
    'before': np.trapezoid(beta_crystal * Ntot_gradient, grid_z) * (const['h'] * const['c'] / pump['wavelength'])
}

# Overwrite pulse with zeros for laser extraction
pulse = np.zeros(steps['time'])

# Call beta_int for the laser extraction
beta_crystal_l, beta_store_l, pulse, Ntot_gradient = beta_int(
    beta_crystal, pulse, const, crystal, steps, laser, mode, Ntot_gradient=Ntot_gradient)

# Energy density after laser extraction
energy_density['after'] = np.trapezoid(beta_crystal_l * Ntot_gradient, grid_z) * (const['h'] * const['c'] / pump['wavelength'])

# Fit and find intersection of beta with beta_min
beta_fit_comp = np.polyfit(grid_z, beta_crystal, 4)
beta_fit = np.polyval(beta_fit_comp, grid_z)

intersect_coeffs = beta_fit_comp.copy()
intersect_coeffs[-1] -= beta_min[0]  # Modify constant term
roots_beta = np.roots(intersect_coeffs)

intersection = 0
for root in roots_beta:
    if np.isreal(root) and root > 0:
        intersection = np.real(root)
        break

# End timer
elapsed_time = time.time() - start_time

# Calculate local gain
gain_local = -(sigma_abs_L - beta_crystal * (sigma_abs_L + sigma_ems_L)) * Ntot_gradient

# Plot beta distribution in absence of ASE
plt.figure(1)
plt.plot(grid_z, beta_crystal, '-', grid_z_beta_min, beta_min, '--', linewidth=2)
plt.title(r'$\beta$ distribution in absence of ASE', fontsize=16)
plt.xlabel('crystal thickness [cm]', fontsize=16)
plt.ylabel(r'$\beta$', fontsize=16)
plt.xlim([0, crystal['length']])
plt.text(0.02, 0.025, r'$\uparrow \beta_{min}$ for 1030nm', fontsize=14)
plt.text(0.2, 0.3, f"{crystal['doping']}% doping", fontsize=14)
plt.grid(True)
plt.show()

# Plot beta distribution after one pass
plt.figure(4)
plt.plot(grid_z, beta_crystal_l, '-', grid_z_beta_min, beta_min, '--', linewidth=2)
plt.title(r'$\beta$ distribution after one pass')
plt.xlabel('crystal length [cm]')
plt.ylabel(r'$\beta$')
plt.xlim([0, crystal['length']])
plt.text(0.02, 0.025, r'$\uparrow \beta_{min}$ for 1030nm')
plt.text(0.2, 0.2, f"{crystal['doping']}% doping")
plt.grid(True)
plt.show()

# Plot beta distribution and gain
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(grid_z, beta_crystal, 'g-', linewidth=2)
ax2.plot(grid_z, gain_local, 'b-', linewidth=2)

ax1.set_xlabel('crystal thickness [cm]')
ax1.set_ylabel(r'$\beta (z)$', color='g')
ax2.set_ylabel('gain (z)', color='b')
plt.title(r'$\beta$ distribution in absence of ASE', fontsize=16)
ax1.set_xlim([0, crystal['length']])
ax2.set_xlim([0, crystal['length']])
plt.show()
