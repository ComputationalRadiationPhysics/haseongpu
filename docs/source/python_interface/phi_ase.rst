PhiASE
======

``PhiASE`` configures and runs the ASE calculation.  It behaves like a Python
configuration class for the lower-level HASEonGPU host mesh, experiment
parameters, compute parameters, backend selection, and execution mode.

.. code-block:: python

   from HASEonGPU import PhiASE

   phi_ase = PhiASE(
       spectralProperties=spectra,
       minRaysPerSample=1000,
       maxRaysPerSample=1000,
       repetitions=1,
       adaptiveSteps=1,
       mseThreshold=0.005,
       useReflections=True,
       backend="Host_Cpu_CpuSerial",
       parallelMode="single",
       numDevices=1,
       rngSeed=1234,
   )

Run One PhiASE Step
-------------------

ASE can be run once without a ``Simulation`` time loop:

.. code-block:: python

   phi_ase.run(gainMedium=medium, crossSections=spectra)
   result = phi_ase.getResults()

   phi = np.asarray(result.phiAse).reshape(
       medium.get("betaCells").expectedShape,
       order="F",
   )

``run(...)`` constructs the low-level host mesh from ``GainMedium``, creates
the experiment and compute parameter objects, calls HASEonGPU, stores the raw
result, and returns ``self``.  The returned ``result.phiAse`` values correspond
to the ASE flux :math:`\Phi_i` described in the scientific background.

Sampling and Physics Settings
-----------------------------

``minRaysPerSample``
   Minimum number of Monte Carlo rays :math:`N` used for each sample point.

``maxRaysPerSample``
   Maximum number of Monte Carlo rays used for adaptive sampling.

``mseThreshold``
   Target mean squared error threshold for the ASE estimate.

``repetitions``
   Maximum number of repeated phiASE runs using the same number of rays if the MSE target is not reached.
   Since the importance-sampling distribution assigns rays to prisms stochastically, repeating the phiASE
   calculation can improve the estimate without increasing the number of rays per run.

``adaptiveSteps``
   Given the mse threshold is not reached for one phiASE backend run - HASEonGPU will linearly increase the rays per sample.
   This parameter controls how many steps there should be between minRaysPerSample and maxRaysPerSample.

``useReflections``
   Enables or disables reflections at the top and bottom surfaces.

``monochromatic``
   This will enforce the phiASE computation to only consider the first cross-section emission and absorption for its calculation.

``rngSeed``
   Optional unsigned RNG seed for reproducible Monte Carlo ray sampling.  Set
   this explicitly for reproducible runs.  If omitted, the Python wrapper
   initializes a process-local NumPy seed stream from ``np.random.SeedSequence()``
   and draws one unsigned 32-bit backend seed for each ASE invocation.

Backend and Parallel Settings
-----------------------------

``backend``
   Alpaka backend name.  The minimal example uses ``"Host_Cpu_CpuSerial"``
   because it is available in a plain CPU build.  Query the installed build
   with ``AlpakaBackends.all()`` and pass one of the returned strings here.
   See :doc:`../backendSelection`.

``parallelMode``
   Execution mode.  ``"single"`` uses the direct pybind binding.
   ``"mpi"`` uses the internal MPI launcher path, writes the temporary binary
   input files, launches ``calcPhiASE`` with ``mpiexec``, and reads the result
   files back into the ``PhiASE`` result object.

``numDevices``
   Maximum number of devices made available on each node for the compute run.
   In MPI execution, HASEonGPU distributes those devices across the MPI ranks
   that are active on the same node.

``nPerNode``
   MPI launcher setting used when ``parallelMode="mpi"``.  The Python
   interface launches the ``calcPhiASE`` binary with
   ``mpiexec -npernode nPerNode``.  This option controls how many MPI ranks
   are started per node by the launcher; it is not a device count inside one
   rank.  See :doc:`../mpi` for the interaction between ``parallelMode``,
   ``numDevices``, and ``nPerNode``.

``minSampleRange`` and ``maxSampleRange``
   Optional inclusive sample-index range.  When omitted, all beta samples
   :math:`\beta_i` are processed.

``writeVtk``
   Requests VTK output from the lower-level compute path when supported.

Configuration Helpers
---------------------

``PhiASE`` can read settings from a dictionary or YAML file.  This is intended
for run-control values: sampling, convergence, reflection flags, backend
selection, MPI launcher settings, optional sample ranges, and VTK output.
Objects such as ``GainMedium``, ``SpectralDecomposition``, and pump solvers are
still passed from Python.

.. code-block:: python

   phi_ase = PhiASE({"minRaysPerSample": 1000, "backend": "Host_Cpu_CpuSerial"})
   phi_ase = PhiASE.fromYaml(
       "phi_ase.yaml",
       spectralProperties=spectra,
       gainMedium=medium,
   )

A YAML file can keep experiment and compute settings together:

.. code-block:: yaml

   experiment:
     min_rays_per_sample: 100000
     max_rays_per_sample: 1000000
     mse_threshold: 0.05
     repetitions: 2
     adaptive_steps: 4
     use_reflections: true
     monochromatic: false

   compute:
     backend: Host_Cpu_CpuSerial
     parallel_mode: single
     numDevices: 1
     n_per_node: 1
     write_vtk: false
     min_sample_range: 0
     max_sample_range: 999
     rng_seed: 1234

YAML keys may be placed at the top level or under ``phiASE``, ``phi_ase``,
``experiment``, or ``compute``.  If the same setting appears more than once,
``PhiASE`` applies sections in this order: ``phiASE``, ``phi_ase``,
``experiment``, ``compute``, then the top-level mapping.  Explicit keyword
overrides passed to ``fromYaml(...)`` are applied after the file is read.

Accepted setting names are the ``PhiASE`` attribute names plus these aliases:
``minRays`` -> ``minRaysPerSample``, ``maxRays`` ->
``maxRaysPerSample``, ``min_rays_per_sample``, ``max_rays_per_sample``,
``mse_threshold``, ``adaptive_steps``, ``use_reflections``,
``parallel_mode``, ``max_gpus`` -> ``numDevices``, ``n_per_node``,
``write_vtk``, ``min_sample_range``, ``max_sample_range``, and
``rng_seed``.

Loading YAML requires ``PyYAML``.  The editable package installation installs
this dependency from ``pyproject.toml``; source-tree usage must provide it in
the Python environment.

For command-line tools:

.. code-block:: python

   parser = PhiASE.addArguments(parser)
   args = parser.parse_args()
   phi_ase = PhiASE.fromArgs(args, spectralProperties=spectra)

The command-line helper accepts ``--phi-ase-config`` first and then applies
explicit command-line options such as ``--backend`` or
``--min-rays-per-sample`` as overrides.  It also accepts ``--rng-seed`` for
reproducible Monte Carlo sampling.

Inspection After a Run
----------------------

After ``run(...)``:

.. code-block:: python

   phi_ase.hostMesh
   phi_ase.experimentParameters
   phi_ase.computeParameters
   phi_ase.getResults()

``getResults()`` raises ``RuntimeError`` if the object has not been run yet.
