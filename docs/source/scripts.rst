Scripts
========

The ``scripts`` directory contains checkout-only helpers for inspecting example
and simulation output.  They are not part of the Python package API.

Small-Signal Gain Plot
----------------------

``scripts/plot_ssg.py`` reads ``laserPumpCladding_*.vtk`` files, integrates a
point-data scalar field along one coordinate direction, and writes a CSV plus a
PNG plot of the resulting net gain factor over time.

For the default ``laserPumpCladding.py`` output, run for example:

.. code-block:: bash

   python3 scripts/plot_ssg.py \
       --input-dir example/python_example \
       --direction z \
       --x 0 \
       --y 0 \
       --output-prefix scripts/gain_z_origin

By default the script integrates the ``gain`` field and applies the round-trip
reflection convention used by the legacy gain workflow:
``net_gain = reflectivity * single_pass_gain**2``.  Use ``--field`` to select a
different VTK scalar, ``--no-back-reflection`` for single-pass gain, and
``--show`` to display the figure interactively after writing it.
