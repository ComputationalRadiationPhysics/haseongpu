# Copyright 2013 Daniel Albach, Erik Zenker, Carlchristian Eckert
# Copyright 2026 Tim Hanel
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

"""Pump integration routines for evolving the excited-state fraction beta.

The public entry point is ``integrateLaserPump``. The lower functions are the
legacy beta-int3 algorithm, kept as implementation detail and accelerated with
numba when available.
"""

import numpy as np
try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        def decorator(fn):
            return fn
        return decorator


CONSTANT_SPECS = {
     "speedOfLight": {
          "description": "Speed of light in vacuum.",
          "unit": "m s^-1",
          "default": np.float64(299792458),
     },
     "planckConstant": {
          "description": "Planck constant.",
          "unit": "J s",
          "default": np.float64(6.62607015e-34),
     },
}

CONSTANT_ALIASES = {
     "c": "speedOfLight",
     "h": "planckConstant",
}


class Constants:
     """Physical constants used by the pump energy-to-photon conversion."""

     def __init__(
          self,
          speedOfLight=CONSTANT_SPECS["speedOfLight"]["default"],
          planckConstant=CONSTANT_SPECS["planckConstant"]["default"],
          *,
          c=None,
          h=None,
     ):
          """Create constants, accepting either long names or ``c``/``h`` aliases."""
          if c is not None:
               speedOfLight = c
          if h is not None:
               planckConstant = h
          self.speedOfLight = float(speedOfLight)
          self.planckConstant = float(planckConstant)

     @property
     def c(self):
          return self.speedOfLight

     @c.setter
     def c(self, value):
          self.speedOfLight = float(value)

     @property
     def h(self):
          return self.planckConstant

     @h.setter
     def h(self, value):
          self.planckConstant = float(value)

     def describeConstant(self, name):
          """Return description, unit, and default for a known constant."""
          canonical = CONSTANT_ALIASES.get(name, name)
          try:
               return dict(CONSTANT_SPECS[canonical])
          except KeyError as exc:
               known = ", ".join([*CONSTANT_SPECS, *CONSTANT_ALIASES])
               raise KeyError(f"unknown physical constant '{name}'. Known constants: {known}") from exc

     def toDict(self):
          """Return the compact dictionary consumed by the legacy pump kernel."""
          return {
              "c": float(self.speedOfLight),
              "h": float(self.planckConstant),
          }


def _constantsDict(constants):
     """Normalize ``None``, ``Constants``, or mapping constants to a dictionary."""
     if constants is None:
          return Constants().toDict()
     if isinstance(constants, Constants):
          return constants.toDict()
     return constants

# Numba-compiled hot loop over time samples and crystal z steps.
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
      """Propagate pump intensity and update beta along one crystal line."""
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
    """Run the one-dimensional pump kernel for every transverse mesh point.

    The super-Gaussian factor turns the global pump intensity into a local
    intensity for each ``(x, y)`` sample before z propagation starts.
    """
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
def beta_int3Main(p, beta_cell, pump, beta_c_2, intensity, const, crystal, steps, int_field, mode, Ntot_gradient):
     """Run one pump update over all transverse points and z levels.

     ``beta_cell`` is the incoming point-level beta array. ``beta_c_2`` receives
     the beta state after pumping. Scratch arrays are allocated once here and
     reused by the inner loops for performance.
     """
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


def integrateLaserPump(
     points,
     betaCells,
     pumpProperties,
     gainMedium,
     timeStep,
     constants=None,
     pumpSubsteps=100,
):
     """Return ``betaCells`` after applying the configured pump.

     ``points`` are transverse ``(x, y)`` topology coordinates. ``betaCells``
     has shape ``(numberOfPoints, numberOfLevels)`` and stores the excited-state
     fraction before pumping. The returned array has the same shape.
     """
     constants = _constantsDict(constants)
     topology = gainMedium.topology
     beta_cells = np.ascontiguousarray(np.asarray(betaCells, dtype=np.float64))
     beta_after_pump = np.empty_like(beta_cells)
     n_tot = float(gainMedium.get("nTot").value)
     pump = pumpProperties.toDict(timeFrame=timeStep)
     tau = (
         float(gainMedium.get("crystalTFluo").value)
         if pumpProperties.temporaryFluorescence is None
         else float(pumpProperties.temporaryFluorescence)
     )
     crystal = {
         "tfluo": tau,
         "length": float(topology.thickness * (topology.levels - 1)),
     }
     steps = {"time": int(pumpProperties.pumpSubsteps if pumpSubsteps is None else pumpSubsteps), "crys": int(topology.levels)}
     mode = pumpProperties.modeDict()
     n_tot_gradient = np.full(topology.levels, n_tot, dtype=np.float64)
     beta_int3Main(
         p=np.asarray(points, dtype=np.float64),
         beta_cell=beta_cells,
         pump=pump,
         beta_c_2=beta_after_pump,
         intensity=float(pumpProperties.intensity),
         const=constants,
         crystal=crystal,
         steps=steps,
         int_field=pump,
         mode=mode,
         Ntot_gradient=n_tot_gradient,
     )
     return beta_after_pump


runLaserPumpStep = integrateLaserPump


def _pumpSpectralSamples(pumpProperties):
     """Return wavelength, sigma_abs, sigma_ems, and raw spectral weights arrays."""
     spectra = pumpProperties.crossSections
     explicit_wavelengths = pumpProperties.getProperty("wavelengths")
     multichromatic = bool(pumpProperties.getProperty("multichromatic", False))
     explicit_wavelength = bool(getattr(pumpProperties, "_explicitWavelength", True))

     if explicit_wavelengths is not None:
          wavelengths = np.asarray(explicit_wavelengths, dtype=np.float64).reshape(-1)
          sigma_abs = np.asarray([spectra.absorptionAt(wavelength) for wavelength in wavelengths], dtype=np.float64)
          sigma_ems = np.asarray([spectra.emissionAt(wavelength) for wavelength in wavelengths], dtype=np.float64)
     elif multichromatic or not explicit_wavelength:
          wavelengths = np.asarray(spectra.wavelengthsAbsorption, dtype=np.float64).reshape(-1)
          sigma_abs = np.asarray(spectra.crossSectionAbsorption, dtype=np.float64).reshape(-1)
          sigma_ems = np.asarray([spectra.emissionAt(wavelength) for wavelength in wavelengths], dtype=np.float64)
     else:
          wavelengths = np.asarray([float(pumpProperties.wavelength)], dtype=np.float64)
          sigma_abs = np.asarray([spectra.absorptionAt(float(pumpProperties.wavelength))], dtype=np.float64)
          sigma_ems = np.asarray([spectra.emissionAt(float(pumpProperties.wavelength))], dtype=np.float64)

     if wavelengths.size == 0:
          raise ValueError("one-dimensional z-traversal pump requires at least one wavelength sample")

     weights = pumpProperties.getProperty("spectralWeights")
     if weights is None:
          weights = np.ones(wavelengths.shape, dtype=np.float64)
     else:
          weights = np.asarray(weights, dtype=np.float64).reshape(-1)
          if weights.size != wavelengths.size:
               raise ValueError("spectralWeights must have the same length as pump wavelength samples")
          if not np.all(np.isfinite(weights)):
               raise ValueError("spectralWeights must be finite")
     return wavelengths, sigma_abs, sigma_ems, weights


def _directionSign(direction):
     if direction is None:
          raise ValueError("OneDimensionalZTraversal requires propagationDirection")
     if isinstance(direction, str) or np.isscalar(direction):
          raise TypeError("OneDimensionalZTraversal propagationDirection must be a 3-vector")
     values = np.asarray(direction, dtype=np.float64).reshape(-1)
     if values.size != 3:
          raise ValueError("OneDimensionalZTraversal propagationDirection must be a 3-vector")
     if not np.all(np.isfinite(values)):
          raise ValueError("OneDimensionalZTraversal propagationDirection must be finite")
     norm = float(np.linalg.norm(values))
     if norm == 0.0:
          raise ValueError("OneDimensionalZTraversal propagationDirection must not be zero")
     z_component = float(values[2] / norm)
     if z_component > 0.0:
          return 1
     if z_component < 0.0:
          return -1
     raise ValueError("OneDimensionalZTraversal currently requires propagationDirection with a non-zero z component")


def oneDimensionalZTraversalPumpRate(points, betaCells, pumpProperties, gainMedium, constants=None):
     """Return frozen-state pump contribution using a 1D z traversal model."""
     constants = _constantsDict(constants)
     topology = gainMedium.topology
     topology._require_levels()
     beta_cells = np.asarray(betaCells, dtype=np.float64).reshape(
          (topology.numberOfPoints, topology.levels),
          order="F",
     )
     points = np.asarray(points, dtype=np.float64)
     wavelengths, sigma_abs, sigma_ems, weights = _pumpSpectralSamples(pumpProperties)
     direction = _directionSign(pumpProperties.propagationDirection)

     transverse_profile = str(pumpProperties.getProperty("transverseProfile", "superGaussian"))
     if transverse_profile != "superGaussian":
          raise ValueError("OneDimensionalZTraversal currently supports transverseProfile='superGaussian'")
     rx = float(pumpProperties.radiusX)
     ry = float(pumpProperties.radiusY)
     center = np.asarray(pumpProperties.getProperty("center", (0.0, 0.0)), dtype=np.float64).reshape(-1)
     if center.size not in {2, 3}:
          raise ValueError("OneDimensionalZTraversal center must be a 2D or 3D point")
     x_rel = points[:, 0] - center[0]
     y_rel = points[:, 1] - center[1]
     super_gaussian_order = float(pumpProperties.superGaussianOrder)
     r = np.sqrt((x_rel ** 2) / (ry ** 2) + (y_rel ** 2) / (rx ** 2))
     pump_input = float(pumpProperties.intensity) * np.exp(-(r ** super_gaussian_order))
     spectral_input = pump_input[:, None] * weights[None, :]

     sigma_sum = sigma_abs + sigma_ems
     n_tot = float(gainMedium.get("nTot").value)
     crystal_step = float(topology.thickness)
     beta_average = 0.5 * (beta_cells[:, :-1] + beta_cells[:, 1:])
     propagation = np.exp(
          -(
              sigma_abs[None, :, None]
              - beta_average[:, None, :] * sigma_sum[None, :, None]
          )
          * n_tot
          * crystal_step
     )

     pump_primary = np.empty((beta_cells.shape[0], wavelengths.size, beta_cells.shape[1]), dtype=np.float64)
     if direction == 1:
          pump_primary[:, :, 0] = spectral_input
          pump_primary[:, :, 1:] = spectral_input[:, :, None] * np.cumprod(propagation, axis=2)
     else:
          pump_primary[:, :, -1] = spectral_input
          pump_primary[:, :, :-1] = spectral_input[:, :, None] * np.cumprod(propagation[:, :, ::-1], axis=2)[:, :, ::-1]

     pump_total = pump_primary
     if bool(pumpProperties.backReflection):
          reflected = np.empty_like(pump_primary)
          reflectivity = float(pumpProperties.reflectivity)
          if direction == 1:
               reflected[:, :, -1] = pump_primary[:, :, -1] * reflectivity
               reflected[:, :, :-1] = reflected[:, :, -1:] * np.cumprod(propagation[:, :, ::-1], axis=2)[:, :, ::-1]
          else:
               reflected[:, :, 0] = pump_primary[:, :, 0] * reflectivity
               reflected[:, :, 1:] = reflected[:, :, :1] * np.cumprod(propagation, axis=2)
          pump_total = pump_total + reflected

     # Φₖ = Iₖ λₖ / (h c), ∂β/∂t|pump = Σₖ [σₐₖ - β(σₐₖ + σₑₖ)] Φₖ.
     photon_flux_scale = wavelengths / (float(constants["h"]) * float(constants["c"]))
     gain_per_flux = sigma_abs[None, :, None] - beta_cells[:, None, :] * sigma_sum[None, :, None]
     rate_by_wavelength = gain_per_flux * pump_total * photon_flux_scale[None, :, None]
     return np.sum(rate_by_wavelength, axis=1)


class OneDimensionalZTraversal:
     """Pump solver using z-aligned 1D traversal."""

     def step(self, input, pump):
          """Return beta advanced by one outer time step using instantaneous pump rate."""
          medium = input["_medium"]
          beta_cells = np.asarray(input["betaCell"], dtype=np.float64)
          rate = oneDimensionalZTraversalPumpRate(
               points=medium.topology.points,
               betaCells=beta_cells,
               pumpProperties=pump,
               gainMedium=medium,
               constants=input.get("_constants"),
          )
          return beta_cells + float(input["_timeStep"]) * rate


class BetaIntegrationGaussianSolver:
     """Default ``PumpProperties`` solver for a super-Gaussian pump beam."""

     def step(self, input, pump):
          """Return beta after pumping using the Simulation pump-solver protocol."""
          medium = input["_medium"]
          timeStep = input["_timeStep"]
          return integrateLaserPump(
              points=medium.topology.points,
              betaCells=input["betaCell"],
              pumpProperties=pump,
              gainMedium=medium,
              timeStep=timeStep,
              constants=input.get("_constants"),
              pumpSubsteps=input.get("_substeps"),
          )


BetaIntegrationSolver = BetaIntegrationGaussianSolver
BetaInt3PumpSolver = BetaIntegrationGaussianSolver
