Interpret forward ASE uncertainty
=================================

Forward ASE results include a relative standard error (RSE) for every cell.
Use this value to decide whether the sampled ray histories provide sufficient
precision for your simulation.

Set an RSE target
-----------------

Configure the target and the available ray range with ``PhiASE``:

.. code-block:: python

   phi_ase = PhiASE(
       minRays=100_000,
       maxRays=1_000_000,
       adaptiveSteps=4,
       relativeStandardErrorThreshold=0.05,
   )

The solver adds ray histories until every cell reaches the target or
``maxRays`` is reached. A threshold of ``0.05`` requests an estimated
one-standard-error uncertainty of 5% relative to the cell value.

Read the result
---------------

After the run, inspect ``result.relativeStandardError`` together with
``result.phiAse`` and ``result.totalRays``. Cells that remain above the target
at ``maxRays`` need a larger ray budget if lower sampling uncertainty is
required.

How the RSE is estimated
------------------------

HASEonGPU divides the ray histories into independent statistical batches. Each
batch samples the complete source and wavelength distributions and produces a
complete cell-field estimate. The variation between these estimates supplies
the reported RSE.

Batching does not change the physical ASE field. The batch contributions are
combined and normalized by the total number of ray histories.

Limitations
-----------

RSE measures Monte Carlo sampling uncertainty. It does not measure mesh
discretization error, uncertainty in material data, or model error. For
reflected ASE, the batches share the surface reflection reservoir, so the
reported RSE does not independently resample uncertainty introduced while
constructing that reservoir.
