Performance and scaling
=======================

.. _benchmark_examples:

Choose the API based on the shape of the problem:

* Use ``Simulation`` for one or a few configurations and for interactive
  exploration.
* Use ``Experiment`` for Cartesian parameter grids and repeated measurements.
* Use ``get(...).as_numpy()`` when plotting and unit metadata are not
  needed; this avoids DataFrame construction.

Grid size and memory
--------------------

An experiment grid contains the Cartesian product of all non-singleton
parameters. A sweep with dimensions ``(100, 50, 20)`` has 100,000 simulation
configurations, even if each individual calculation is inexpensive. Check the
shape before running::

   print(experiment.array_shape)
   print(experiment.total_iterations)

Start with a small grid, verify the physics, then increase resolution. Split
very large sweeps into chunks when their complete DataFrame or NumPy array does
not fit comfortably in memory.

Timing fairly
-------------

For reproducible timing:

* run several repetitions and report the median;
* exclude imports, plotting, and notebook display from the timed region;
* record Python, PyMieSim, NumPy, compiler, operating-system, and CPU details;
* use the same parameter values and requested measures for every comparison;
* warm up the first call before collecting timings.

The runnable :ref:`benchmark examples <benchmark_examples>` follow this
pattern. Sphinx-Gallery also exports them as notebooks.

Output choices
--------------

``.as_numpy()`` is usually fastest for numeric post-processing. DataFrame
output is available through ``.as_dataframe()`` when you need tabular
interoperability. Request only the measures needed for the study, especially
for large grids.

Far fields and near fields
--------------------------

Field representations can dominate runtime and memory because they add an
angular or spatial sampling dimension. Begin with a coarse mesh, inspect
convergence, and increase sampling only after the desired features are
resolved. Do not include plotting in solver benchmarks.

Reproducible examples
---------------------

* :download:`parameter-grid benchmark <../examples/benchmarks/parameter_grid.py>`
* :download:`sweep-runtime benchmark <../examples/benchmarks/sweep_runtime.py>`
* :download:`reproducible sweep benchmark <../examples/benchmarks/reproducible_parameter_sweep.py>`
