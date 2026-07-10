Legacy Low-Level Python Interface
=================================

This page documents the compatibility status of the historical low-level
``calcPhiASE(...)`` Python entry point. New workflows should use the high-level
:doc:`Python Interface Guide <pythonInterface>` and configure ASE through
``PhiASE``.

Current Compatibility
---------------------

``from HASEonGPU import calcPhiASE`` remains available for direct in-process
Python binding calls. It accepts the historical many-argument signature and
normalizes flat, matrix, NumPy, and list inputs before constructing the C++
``HostMesh``, ``ExperimentParameters``, and ``ComputeParameters`` objects.

The removed legacy text-file backend parser is not restored. In particular:

* ``parallelMode="single"`` uses the direct Python binding path.
* ``parallelMode="debugFileIOPath"`` is accepted for compatibility but is
  routed to the direct binding path with a warning.
* ``parallelMode="mpi"`` no longer launches the retired text-file executable
  path; use the openPMD transport path instead.

Modern Execution Path
---------------------

``PhiASE.run(...)`` serializes ``GainMedium``, spectra, and compute settings to
the HASE openPMD transport, launches the compiled ``calcPhiASE`` backend, and
reads ``core_result_*`` records back into Python. The standalone binary accepts
only ``--input-path`` and ``--output-path`` and reads simulation settings from
the openPMD input series.

See :doc:`openpmdTransport` for the openPMD record layout, storage backend
options, dynamic-iteration behavior, and MPI launch examples.

Migration Guide
---------------

Prefer replacing direct argument lists with domain objects:

* ``MeshTopology`` owns points, triangle connectivity, levels, and thickness.
* ``GainMedium`` owns ``betaCells``, ``betaVolume``, cladding data, refractive
  indices, reflectivities, ``nTot``, and ``crystalTFluo``.
* ``CrossSectionData`` or ``SpectralDecomposition`` owns wavelength and
  cross-section tables.
* ``PhiASE`` owns sampling, convergence, reflection, backend, sample-range, and
  optional RNG seed settings.

For workflow examples, start with :doc:`pythonInterface`. For generated
signatures, use :doc:`pythonAPI`.
