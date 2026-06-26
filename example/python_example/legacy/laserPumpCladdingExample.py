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
import argparse
import time
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat
from scipy.interpolate import griddata
from beta_int3 import beta_int3Main
from set_variables import set_variables
from vtk_wedge import vtk_wedge
import csv


Legacy_DIR = Path(__file__).resolve().parent
SCRIPT_DIR= Legacy_DIR
from HASEonGPU import calcPhiASE


def interpolate_beta_volume(beta_cell, p, x_center, y_center, number_of_triangles, mesh_z, z_mesh):
    beta_vol = np.zeros((number_of_triangles, mesh_z - 1))
    x_1 = p[:, 0]
    y_1 = p[:, 1]
    x_2 = p[:, 0]
    y_2 = p[:, 1]
    x = np.concatenate((x_1, x_2))
    y = np.concatenate((y_1, y_2))
    z = np.concatenate((np.zeros(x_1.shape[0]), np.zeros(x_2.shape[0]) + z_mesh))
    xi = x_center
    yi = y_center
    zi = np.zeros(xi.shape) + z_mesh / 2

    for i_z in range(mesh_z - 1):
        v = np.concatenate((beta_cell[:, i_z], beta_cell[:, i_z + 1]))
        beta_vol[:, i_z] = griddata((x, y, z), v, (xi, yi, zi), method='linear')
        z = z + z_mesh
        zi = zi + z_mesh
    return beta_vol


def main(
    material_path: Path = SCRIPT_DIR / "pt.mat",
    gpus: int = 1,
    parallel_mode: str = "single",
    backend_override: str | None = None,
    min_rays: int | None = None,
    max_rays: int | None = None,
    use_reflections_override: bool | None = None,
    repetitions_override: int | None = None,
    adaptive_steps_override: int | None = None,
    timeslice_override: int | None = None,
    timeslice_total_override: int | None = None,
    min_sample_range: int = 0,
    max_sample_range: int | None = None,
    rng_seed: int | None = None,
    mse_threshold_override: float | None = None,
    snapshot_path: Path | None = None,
):
    print(f' root: {SCRIPT_DIR}')
    #precheck the file path and throw in case this is not present
    # Crystal parameter
    crystal = {
        'doping': 2,
        'length': 0.7,    # [cm]
        'tfluo': 9.5e-4,  # 1.0e-3 pour DA
        'nlexp': 1,
        'levels': 10,
    }
    # nonlinearity of the doping - exponential - factor front/end of the
    # crystal

    # Timesteps
    steps = {
        'time': 100,
        'crys': crystal['levels'],
    }

    # Pump parameter
    pump_stretch = 1
    pump_aspect = 1
    pump_diameter = 3
    pump = {
        's_abs': 0.778e-20,  # Absorption cross section in cm^2 (0.778e-20 pour DA)
        's_ems': 0.195e-20,  # Emission cross section in cm^2 (0.195e-20 pour DA)
        'I': 16e3,  # Pump intensity in W/cm^2
        'T': 1e-3,  # Pump duration in s
        'wavelength': 940e-9,  # Pump wavelength in m
        'ry': pump_diameter / 2 * pump_aspect,
        'rx': pump_diameter / 2,  # Pump radius
        'exp': 40,  # Pump exponent
    }
    # Laser parameter
    laser = {
        's_abs' : np.loadtxt(material_path.parent.parent/'input'/'sigma_a.txt'),  # Absorption spectrum cm2(1.16e-21 pour DA)
        's_ems' : np.loadtxt(material_path.parent.parent/'input'/'sigma_e.txt'), # Emission spectrum in cm2(2.48e-20)
        'l_abs' : np.loadtxt(material_path.parent.parent/'input'/'lambda_a.txt'), # Wavelengths absorption spectrum in nm (x values for absorption)
        'l_ems' : np.loadtxt(material_path.parent.parent/'input'/'lambda_e.txt'), # Wavelengths emission spectrum in nm (y values for emission)
        'l_res' : 1000,                      # Resolution of linear interpolated spectrum
        'I' : 1e6,                           # Laser intensity
        'T' : 1e-8,                          # Laser duration
        'wavelength' : 1030e-9,              # Lasing wavelength in m
    }
    laser['max_ems'] = np.max(laser['s_ems'])
    laser['i'] = np.argmax(laser['s_ems'])
    laser['max_abs']= laser['s_abs'][laser['i']]

    # Mode parameter
    mode = {
        'BRM': 1,  # 1 Backreflection, 0 Transmissionmode
        'R': 1,    # reflectivity of the coating
        'extr': 0,  # no extraction!!
    }

    # Constants
    phy_const= {
        'N1per': 1.38e20,
        'c': 3e8,
        'h': 6.626e-34,
    }

    N_tot = phy_const['N1per'] * crystal['doping']
    Ntot_gradient = np.zeros(crystal['levels'])   # doping gradient if not given by load!
    Ntot_gradient[:] = crystal['doping'] * phy_const['N1per']
    mesh_z = crystal['levels']
    z_mesh = crystal['length'] / (mesh_z-1)
    timeslice = 50 if timeslice_override is None else int(timeslice_override)
    timeslice_tot = 150 if timeslice_total_override is None else int(timeslice_total_override)
    timetotal = 1e-3  # [s]
    time_t = timetotal / timeslice

    # ASE application
    numDevices = gpus
    nPerNode = 1
    backend = 'gpu' if backend_override is None else str(backend_override)
    # parallelMode = 'mpi'
    # parallelMode = 'single'
    parallelMode = parallel_mode
    useReflections = True if use_reflections_override is None else bool(use_reflections_override)
    refractiveIndices = [1.83, 1, 1.83, 1]
    repetitions = 4 if repetitions_override is None else int(repetitions_override)
    adaptiveSteps = 5 if adaptive_steps_override is None else int(adaptive_steps_override)
    minRaysPerSample = 1e5 if min_rays is None else int(min_rays)
    maxRaysPerSample = minRaysPerSample * 100 if max_rays is None else int(max_rays)
    # this is the volume weighting adjusted threshold (scaled the original 0.005 by thickness * mesh.surface )
    mseThreshold = 0.01087 if mse_threshold_override is None else float(mse_threshold_override)

    ## this is our starting counter
    tic= time.perf_counter()

    # Constants for short use
    c = phy_const['c']  # m/s
    h = phy_const['h']  # Js

    ################################ Create mesh ######################################%
    # Load the grid from file
    # the grid is done using the distmesh routines
    # use e.g. cr_60mm_30mm.m in meshing folder
    print("loading points data")
    # load data from .mat files
    pt_data = loadmat(str(material_path))
    print("finished loading points data")
    # extract necessary variables
    p = np.array(pt_data['p'])
    print(f'points shape from loadmat {p.shape[0]}')
    t = np.array([[x -1 for x in y] for y in  pt_data['t']]) # indes shift due to earlier matlab konvention
    print("t min/max:", t.min(), t.max(), "t shape:", t.shape)

    set_variables(p, t)

    variable_data = loadmat('variable.mat')

    normals_x= variable_data['normals_x']
    normals_y= variable_data['normals_y']
    normals_z= variable_data['normals_z']
    ordered_int= variable_data['sorted_int']
    t_int= variable_data['t_int']
    x_center= variable_data['x_center'][0]
    y_center= variable_data['y_center'][0]
    surface= variable_data['surface']
    normals_p= variable_data['normals_p']
    forbidden=variable_data['forbidden']

    beta_cell = np.zeros((p.shape[0], mesh_z))
    beta_vol = np.zeros((t.shape[0], mesh_z-1))
    nT = ordered_int.shape[0]
    reflectivities = np.zeros((1, nT*2))

    # initialize arrays
    phi_ASE = np.zeros((p.shape[0], mesh_z))
    dndt_ASE = np.copy(phi_ASE)
    flux_clad = np.copy(phi_ASE)
    dndt_pump = np.copy(beta_cell)

    # save initial values to file
    vtk_wedge('beta_cell_0.vtk', beta_cell, p, t_int, mesh_z, z_mesh)
    vtk_wedge('dndt_pump_0.vtk', dndt_pump, p, t_int, mesh_z, z_mesh)
    vtk_wedge('dndt_ASE_0.vtk', dndt_ASE, p, t_int, mesh_z, z_mesh)
    savemat('save_0.mat', {'pys_const': phy_const, 'mode': mode, 'steps': steps, 'p': p })

    # simulate pumping and ASE for each point in the mesh
    temp = pump['T']
    time_beta = 1e-6
    pump['T'] = time_beta
    temp_f = crystal['tfluo']
    crystal['tfluo'] = 1
    beta_c_2 = np.zeros((p.shape[0], mesh_z))
    intensity = pump['I']
    # first of beta int
    beta_int3Main(
        beta_c_2=beta_c_2,
        p=p,
        beta_cell=beta_cell,
        pump=pump,
        intensity=intensity,
        const=phy_const,
        crystal=crystal,
        steps=steps,
        int_field=pump,
        mode=mode,
        Ntot_gradient=Ntot_gradient
    )

    # calculate change in beta_cell due to pumping
    dndt_pump = (beta_c_2 - beta_cell) / time_beta

    # reset parameters to initial values
    pump['I'] = intensity
    pump['T'] = temp
    crystal['tfluo'] = temp_f
    time_int = 1e-6

    for i_p in range(p.shape[0]):
        for i_z in range(mesh_z):
            beta_cell[i_p, i_z] = crystal['tfluo'] * (dndt_pump[i_p,i_z] - dndt_ASE[i_p,i_z]) * (1-np.exp(-time_t/crystal['tfluo'])) + beta_cell[i_p,i_z]*np.exp(-time_t/crystal['tfluo'])

    vtk_wedge('beta_cell_1.vtk', beta_cell, p, t_int, mesh_z, z_mesh)
    initial_beta_cell = beta_cell.copy()

    ################################ cladding ###########################

    clad_number = 1
    clad_number = np.int32(clad_number)
    numberOfTriangles, b = t_int.shape
    clad = np.zeros((numberOfTriangles,1))

    # Absorption of the cladding
    clad_abs = 5.5


    clad_index = np.where(clad==clad_number)[0]
    Yb_index = np.where(clad!=clad_number)[0]

    print("phi_ASE.shape[0]:", phi_ASE.shape[0])

    Yb_points_t = t.copy()
    Yb_points_t = np.delete(Yb_points_t, clad_index, axis=0)
    print("Yb_points_t min/max:", Yb_points_t.min(), Yb_points_t.max(), "len(unique):", len(Yb_points_t))

    Yb_points_t2 = np.reshape(Yb_points_t,(-1,1))
    del Yb_points_t
    Yb_points_t = np.unique(Yb_points_t2)
    del Yb_points_t2

    clad_points_t = t.copy()
    clad_points_t = np.delete(clad_points_t, Yb_index, axis=0)
    clad_points_t2 = np.reshape(clad_points_t,(-1,1))
    del clad_points_t
    clad_points_t = np.unique(clad_points_t2)
    del clad_points_t2

    #Now get the lengths of the
    length_clad = clad_points_t.shape[0]
    length_Yb = Yb_points_t.shape[0]

    #Make sure that clad is integer variable
    clad_int = np.int32(clad)
    slice_times_series=[]
    ################################ Main pump loop ###################################

    for i_slice in range(0, timeslice_tot):
        start_slice_t=time.perf_counter()
        print('TimeSlice', i_slice, 'of', timeslice_tot-1, 'started')
        # ******************* BETA PUMP TEST ******************************
        # make a test with the gain routine "gain.m" for each of the points
        # define the intensity at the nominal points of the grid and make the
        # pumping to get the temporal evolution of the beta
        # make a rectangular field, make the interpolation of the points on it
        # if we don't have any formula discribing the intensity distribution
        # as first test do a supergaussian distribution for 1 s and then estimate
        # the dndt|pump

        # Write output for this timeslice to file
        file_b = 'beta_cell_' + str(i_slice) + '.vtk'
        file_p = 'dndt_pump_' + str(i_slice) + '.vtk'
        file_A = 'dndt_ASE_' + str(i_slice) + '.vtk'
        file_C = 'flux_clad_' + str(i_slice) + '.vtk'
        vtk_wedge_in_t=time.perf_counter()
        vtk_wedge(file_b, beta_cell, p, t_int, mesh_z, z_mesh)
        vtk_wedge(file_p, dndt_pump, p, t_int, mesh_z, z_mesh)
        vtk_wedge_out_t=time.perf_counter()
        print(f'[TIME] writing vtk files: {vtk_wedge_out_t-vtk_wedge_in_t}')
        # Interpolate beta_vol from beta_cell
        x_1 = p[:,0]
        y_1 = p[:,1]
        x_2 = p[:,0]
        y_2 = p[:,1]

        x = np.concatenate((x_1,x_2))
        y = np.concatenate((y_1,y_2))

        z_1 = np.zeros(x_1.shape[0])
        z_2 = np.zeros(x_2.shape[0])
        z_2 = z_2 + z_mesh

        z = np.concatenate((z_1,z_2))

        xi = x_center
        yi = y_center
        zi = np.zeros(xi.shape)
        zi = zi + z_mesh/2

        for i_z in range(mesh_z-1):
            v_1 = beta_cell[:,i_z]
            v_2 = beta_cell[:,i_z+1]

            v = np.concatenate((v_1,v_2))
            pts= np.array([x, y, z])
            ges= np.array([xi,yi,zi])

            beta_vol[:,i_z] = griddata((x,y,z), v, (xi,yi,zi), method= 'linear')
            #beta_vol[:,i_z] = griddata(pts, v, ges)

            z = z + z_mesh
            zi = zi + z_mesh
        phi_ASECall_in_t=time.perf_counter()
        ############################ Call external ASE application ####################
        phi_ASE, mse_values, N_rays = calcPhiASE(
            p,
            t_int,
            beta_cell,
            beta_vol,
            clad_int,
            clad_number,
            clad_abs,
            useReflections,
            refractiveIndices,
            reflectivities,
            normals_x,
            normals_y,
            ordered_int,
            surface,
            x_center,
            y_center,
            normals_p,
            forbidden,
            minRaysPerSample,
            maxRaysPerSample,
            mseThreshold,
            repetitions,
            N_tot,
            z_mesh,
            laser,
            crystal,
            mesh_z,
            backend,
            parallelMode,
            numDevices,
            adaptiveSteps,
            nPerNode,
            minSampleRange=min_sample_range,
            maxSampleRange=max_sample_range,
            rngSeed=rng_seed,
        )
        phi_ASECall_out_t=time.perf_counter()
        print(f'[TIME] Total extern phi call time: {phi_ASECall_out_t-phi_ASECall_in_t}')

        dndt_ASE.fill(0.0)
        flux_clad.fill(0.0)

        # Calc dn/dt ASE, change it to only the Yb:YAG points
        for i_p in range(length_Yb):
            for i_z in range(mesh_z):
                pos_Yb = Yb_points_t[i_p]
                # calcPhiASE returns phi_ASE already scaled by N_tot / crystal.tfluo.
                g_l = beta_cell[pos_Yb,i_z] * (laser['max_ems'] + laser['max_abs']) - laser['max_abs']
                dndt_ASE[pos_Yb,i_z] = g_l * phi_ASE[pos_Yb,i_z]

        # Now for the cladding points
        for i_p in range(length_clad):
            for i_z in range(mesh_z):
                # Calc local gain (g_l)
                pos_clad = clad_points_t[i_p]
                flux_clad[pos_clad,i_z] = clad_abs*phi_ASE[pos_clad,i_z]

        ############################ Prepare next timeslice ###########################
        temp = pump['T']
        time_beta = 1e-6
        pump['T'] = time_beta
        temp_f = crystal['tfluo']
        crystal['tfluo'] = 1

        beta_c_2 = np.zeros((p.shape[0], mesh_z))
        intensity = pump['I']
        beta_int3Main(
            beta_c_2=beta_c_2,
            p=p,
            beta_cell=beta_cell,
            pump=pump,
            intensity=intensity,
            const=phy_const,
            crystal=crystal,
            steps=steps,
            int_field=pump,
            mode=mode,
            Ntot_gradient=Ntot_gradient
        )

        dndt_in_t = time.perf_counter()
        if i_slice < timeslice:
            dndt_pump = np.divide(np.subtract(beta_c_2, beta_cell), time_beta)
        else:
            dndt_pump = np.zeros((p.shape[0], mesh_z))

        # optional NaN/inf checks
        if not np.isfinite(beta_c_2).all():
            print('[WARN] beta_c_2 contains non-finite values (NaN/Inf)')
        if not np.isfinite(dndt_pump).all():
            print('[WARN] dndt_pump contains non-finite values (NaN/Inf)')

        # restore pump/crystal state
        pump['I'] = intensity
        pump['T'] = temp
        crystal['tfluo'] = temp_f

        # update beta_cell
        tfluo = crystal['tfluo']
        exp_fac = np.exp(-time_t / tfluo)
        one_minus = (1.0 - exp_fac)

        for i_p in range(p.shape[0]):
            for i_z in range(mesh_z):
                beta_cell[i_p, i_z] = (
                        tfluo * (dndt_pump[i_p, i_z] - dndt_ASE[i_p, i_z]) * one_minus
                        + beta_cell[i_p, i_z] * exp_fac
                )
        beta_cell[:] = np.clip(beta_cell, 0.0, 1.0)

        end_slice_t = time.perf_counter()
        sliceTime=end_slice_t-start_slice_t
        print(f'[TIME] Total time for the current time slice: {sliceTime}')
        slice_times_series.append(sliceTime)

    csv_name = "slice_times.csv"
    with open(csv_name, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["slice", "seconds"])
        for i, dt in enumerate(slice_times_series, start=1):
            w.writerow([i, dt])

    slice_times = np.asarray(slice_times_series, dtype=np.float64)
    print("\n=== Performance Measurements (Timings) ===")
    print(f"slices measured: {slice_times.size}")
    print(f"mean: {slice_times.mean():.3f} s, std: {slice_times.std():.3f} s")
    print(f"min:  {slice_times.min():.3f} s, max: {slice_times.max():.3f} s")
    print(f"total (measured slices): {slice_times.sum():.3f} s")

    print(f"[INFO] Wrote slice time series to {csv_name}")
    file_b = 'beta_cell_' + str(timeslice_tot) + '.vtk'
    file_p = 'dndt_pump_' + str(timeslice_tot) + '.vtk'
    file_A = 'dndt_ASE_' + str(timeslice_tot) + '.vtk'
    file_C = 'flux_clad_' + str(timeslice_tot) + '.vtk'

    vtk_wedge(file_b, beta_cell, p, t_int, mesh_z, z_mesh)
    vtk_wedge(file_p, dndt_pump, p, t_int, mesh_z, z_mesh)
    vtk_wedge(file_A, dndt_ASE, p, t_int, mesh_z, z_mesh)
    vtk_wedge(file_C, flux_clad, p, t_int, mesh_z, z_mesh)
    print('Calculations finished')
    toc = time.perf_counter()
    elapsed_time= toc - tic

    print(f'Total time taken: {elapsed_time} seconds')
if __name__ == "__main__":
    print(SCRIPT_DIR)
    parser = argparse.ArgumentParser(description="HASEonGPU pump/ASE simulation")
    parser.add_argument(
        "--material",
        "-m",
        default=str(SCRIPT_DIR / "pt.mat"),
        help="Path to material/mesh MAT file (default: pt.mat)",
    )
    parser.add_argument(
        "--number_of_gpus",
        "-nr_gpus",
        default=1,
        help="The number of gpus used for phiASE calculation (default: 1)",
    )
    parser.add_argument(
        "--parallel_mode",
        "-par",
        default="single",
        help="The Parallel Mode used for phiASE [single|mpi] (default: single)",
    )
    parser.add_argument(
        "--backend",
        default=None,
        help="Override the phiASE backend for this run (default: legacy script default)",
    )
    parser.add_argument(
        "--min-rays",
        type=int,
        default=None,
        help="Override minRaysPerSample for this run",
    )
    parser.add_argument(
        "--max-rays",
        type=int,
        default=None,
        help="Override maxRaysPerSample for this run",
    )
    parser.add_argument(
        "--repetitions",
        type=int,
        default=None,
        help="Override the number of phiASE repetitions for this run",
    )
    parser.add_argument(
        "--adaptive-steps",
        type=int,
        default=None,
        help="Override the number of phiASE adaptive ray-count steps for this run",
    )
    parser.add_argument(
        "--mse-threshold",
        type=float,
        default=None,
        help="Override the phiASE MSE threshold for this run",
    )
    parser.add_argument(
        "--reflection",
        type=int,
        choices=(0, 1),
        default=None,
        help="Override useReflections for this run (0 or 1)",
    )
    parser.add_argument(
        "--timeslice",
        type=int,
        default=None,
        help="Override the pump-on timeslice count for this run",
    )
    parser.add_argument(
        "--timeslice-total",
        type=int,
        default=None,
        help="Override the total legacy loop count for this run",
    )
    parser.add_argument(
        "--min-sample-range",
        type=int,
        default=0,
        help="First flattened sample index processed by phiASE (default: 0)",
    )
    parser.add_argument(
        "--max-sample-range",
        type=int,
        default=None,
        help="Last flattened sample index processed by phiASE (default: full mesh)",
    )
    parser.add_argument(
        "--rng-seed",
        type=int,
        default=None,
        help="Optional phiASE RNG seed for reproducible comparison runs",
    )
    args = parser.parse_args()
    material_path = Path(args.material)

    # check file path
    if not material_path.exists():
        raise FileNotFoundError(
            f"Material file '{material_path}' does not exist."
        )

    if not material_path.is_file():
        raise ValueError(
            f"Material path '{material_path}' is not a file."
        )
    main(
        material_path=material_path,
        gpus=args.number_of_gpus,
        parallel_mode=args.parallel_mode,
        backend_override=args.backend,
        min_rays=args.min_rays,
        max_rays=args.max_rays,
        use_reflections_override=None if args.reflection is None else bool(args.reflection),
        repetitions_override=args.repetitions,
        adaptive_steps_override=args.adaptive_steps,
        mse_threshold_override=args.mse_threshold,
        timeslice_override=args.timeslice,
        timeslice_total_override=args.timeslice_total,
        min_sample_range=args.min_sample_range,
        max_sample_range=args.max_sample_range,
        rng_seed=args.rng_seed,
    )
