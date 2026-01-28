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

import time
import numpy as np
from numba import njit

"""
    Just in time compilation of the hot loop (over time and pumping steps) using scratch memory arrays
"""
@njit(cache=True, fastmath=True)
def beta_int3Kernel_jit(
          beta_crystal, pulse, beta_store, Ntot_gradient,
          pump, pump_BRM, pump_l,
          I_pump,
          time_step, crystal_step,
          sigma_abs, sigma_ems,
          inv_photon_energy, inv_tau_fluo,
          extr, brm, R
):
      # discretization
      steps_time = pulse.shape[0]
      steps_crystal = beta_crystal.shape[0]

      # declarations
      # sigma_abs = ...  # cm^2
      # sigma_ems = ...  # cm^2

      for itime in range(steps_time):

          if extr == 0:
              pump[0] = I_pump
          else:
              pump[0] = pulse[itime, 0]

          # this is the positive direction
          for icrys in range(steps_crystal - 1):
              # step one is from point one to two for I_pump
              beta_average = (beta_crystal[icrys] + beta_crystal[icrys + 1]) / 2.0
              pump[icrys + 1] = pump[icrys] * np.exp(
                  -(sigma_abs - beta_average * (sigma_abs + sigma_ems))
                  * Ntot_gradient[icrys]
                  * crystal_step
              )

          if brm == 1:
              # (semantic equivalent to flipud/rot90, but without copies)
              #   this is the negative direction
              #   full pump intensity is I+ + I-
              pump_BRM[steps_crystal - 1] = pump[steps_crystal - 1] * R

              for jcrys in range(steps_crystal - 2, -1, -1):
                  # step one is from point one to two for I_pump
                  beta_average = (beta_crystal[jcrys] + beta_crystal[jcrys + 1]) / 2.0
                  pump_BRM[jcrys] = pump_BRM[jcrys + 1] * np.exp(
                      -(sigma_abs - beta_average * (sigma_abs + sigma_ems))
                      * Ntot_gradient[jcrys]
                      * crystal_step
                  )

              # in the case of BRM the return value has to be I-!
              pulse[itime, 0] = pump_BRM[0]
          else:
              # no BRM
              for k in range(steps_crystal):
                  pump_BRM[k] = 0.0
              pulse[itime, 0] = pump[steps_crystal - 1]

          #         full pump intensity is I+ + I-
          for k in range(steps_crystal):
              pump_l[k] = pump[k] + pump_BRM[k]

          #   now calculate the local beta
          for ibeta in range(steps_crystal):
              A1 = sigma_abs * pump_l[ibeta] * inv_photon_energy
              C1 = (sigma_abs + sigma_ems) * pump_l[ibeta] * inv_photon_energy + inv_tau_fluo
              beta_crystal[ibeta] = (A1 / C1) * (1.0 - np.exp(-C1 * time_step)) + beta_crystal[ibeta] * np.exp(-C1 * time_step)

          #     if icrys or jcrys makes no difference
          for k in range(steps_crystal):
              beta_store[k, itime] = beta_crystal[k]
"""
    Runs the outer loop over all sample points.
"""
def beta_int3LoopOverPoints(
          p, beta_cell, pump_dict, intensity, beta_c_2,
          beta_crystal_scratch, pulse_scratch, beta_store_scratch,
          Ntot_gradient,
          # precomputed scalars:
          time_step, crystal_step,
          sigma_abs, sigma_ems,
          inv_photon_energy, inv_tau_fluo,
          extr, brm, R,
          # scratch arrays for kernel:
          pump_scratch, pump_brm_scratch, pump_l_scratch,
          debug=False
):
    rx = float(pump_dict['rx'])
    ry = float(pump_dict['ry'])
    exp_ = float(pump_dict['exp'])

    for i_p in range(p.shape[0]):

          # beta_crystal := beta_cell row (copy into reusable scratch)
          np.copyto(beta_crystal_scratch, beta_cell[i_p, :])

          # pulse := zeros (reuse)
          pulse_scratch.fill(0.0)

          # compute pump I for this point
          r_term = np.sqrt((p[i_p, 0] ** 2) / (ry ** 2) + (p[i_p, 1] ** 2) / (rx ** 2))
          pump_factor = np.exp(-(r_term ** exp_))
          I_pump_point = intensity * pump_factor
          beta_int3Kernel_jit(
              beta_crystal_scratch,
              pulse_scratch,
              beta_store_scratch,
              Ntot_gradient,
              pump_scratch, pump_brm_scratch, pump_l_scratch,
              I_pump_point,
              time_step, crystal_step,
              sigma_abs, sigma_ems,
              inv_photon_energy,
              inv_tau_fluo,
              extr, brm, R
          )

          # write back result for this point
          beta_c_2[i_p, :] = beta_crystal_scratch
"""
    Runs the beta-pumping step for each point.
    This function acts as a hub that initializes constants and fields (once) and setup scratch pad memory arrays.
"""
def beta_int3Main(p, beta_cell, pump, beta_c_2, intensity, const, crystal, steps, int_field, mode, Ntot_gradient):
     # discretization
     steps_time = int(steps['time'])
     steps_crystal = int(steps['crys'])

     # --- extract constants once ---
     h = float(const['h'])
     c = float(const['c'])

     # pump/intensity field constants
     tau_pump = float(int_field['T'])
     wavelength = float(int_field['wavelength'])

     sigma_abs = float(np.max(int_field['s_abs']))
     sigma_ems = float(np.max(int_field['s_ems']))

     # crystal constants
     tau_fluo = float(crystal['tfluo'])
     crystal_length = float(crystal['length'])

     # modes
     extr = int(mode['extr'])
     brm = int(mode['BRM'])
     R = float(mode.get('R', 0.0))

     # --- precompute scalars used in inner loops ---
     time_step = tau_pump / (steps_time - 1)
     crystal_step = crystal_length / (steps_crystal - 1)

     inv_photon_energy = wavelength / (h * c)   # == wavelength/(h*c)
     inv_tau_fluo = 1.0 / tau_fluo

     # Ensure gradient is 1D contiguous float64
     Ntot_gradient = np.ascontiguousarray(np.asarray(Ntot_gradient, dtype=np.float64).reshape(-1))
     if Ntot_gradient.shape[0] < steps_crystal - 1:
         raise ValueError(f"Ntot_gradient too short: got {Ntot_gradient.shape[0]}, need at least {steps_crystal-1}")

     # --- allocate reusable scratch ---

     beta_crystal_scratch = np.empty((steps_crystal,), dtype=np.float64)
     pulse_scratch = np.zeros((steps_time, 1), dtype=np.float64)
     beta_store_scratch = np.empty((steps_crystal, steps_time), dtype=np.float64)

     pump_scratch = np.empty((steps_crystal,), dtype=np.float64)
     pump_brm_scratch = np.empty((steps_crystal,), dtype=np.float64)
     pump_l_scratch = np.empty((steps_crystal,), dtype=np.float64)

     # run outer loop
     beta_int3LoopOverPoints(
         p=p,
         beta_cell=beta_cell,
         pump_dict=pump,
         intensity=float(intensity),
         beta_c_2=beta_c_2,
         beta_crystal_scratch=beta_crystal_scratch,
         pulse_scratch=pulse_scratch,
         beta_store_scratch=beta_store_scratch,
         Ntot_gradient=Ntot_gradient,
         time_step=time_step,
         crystal_step=crystal_step,
         sigma_abs=sigma_abs,
         sigma_ems=sigma_ems,
         inv_photon_energy=inv_photon_energy,
         inv_tau_fluo=inv_tau_fluo,
         extr=extr,
         brm=brm,
         R=R,
         pump_scratch=pump_scratch,
         pump_brm_scratch=pump_brm_scratch,
         pump_l_scratch=pump_l_scratch,
     )