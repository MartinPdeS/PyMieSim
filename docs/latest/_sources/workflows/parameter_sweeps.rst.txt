.. _workflow_parameter_sweeps:

Parameter sweeps
================

Use :class:`PyMieSim.Experiment` and set classes when several values should be
evaluated as a parameter grid.

.. code-block:: python

   import numpy as np

   from PyMieSim import (
       Measure,
       Experiment,
       GaussianSet,
       PolarizationSet,
       SphereSet,
       ureg,
   )

   source = GaussianSet(
       wavelength=np.linspace(400, 800, 5) * ureg.nanometer,
       polarization=PolarizationSet(angles=[0] * ureg.degree),
       optical_power=[1e-3] * ureg.watt,
       numerical_aperture=[0.2],
   )

   scatterer = SphereSet(
       diameter=np.linspace(100, 500, 10) * ureg.nanometer,
       material=[1.5],
       medium=[1.0],
   )

   experiment = Experiment(scatterer_set=scatterer, source_set=source)
   result = experiment.run(Measure.QSCA, Measure.QEXT)
   result["Qsca"].plot(x="source:wavelength", y="Qsca")

This grid contains five wavelengths and ten diameters. The result maps each requested measure to a result whose parameter axes
have lengths five and ten. ``Qsca`` and ``Qext`` are dimensionless. Curves at different diameters
can move spectral resonances even when the refractive index is fixed.

Each entry is a unit-aware ``SimulationResult``.  Use ``.as_numpy()`` when a
raw NumPy array is preferable, ``.as_dataframe()`` for tabular interoperability,
or ``get_sequential`` for aligned sequential configurations.

Extracting NumPy and pandas data
---------------------------------

The experiment result remains labeled by default. Convert it explicitly when
using libraries that expect a raw array or a pandas table:

.. code-block:: python

   result = experiment.run(Measure.QSCA)

   values = result.as_numpy()
   print(values.shape)

   dataframe = result.as_dataframe()
   print(dataframe[["source:wavelength", "scatterer:diameter", "Qsca"]])

For multiple measures, ``as_dataframe()`` creates one output column per measure:

.. code-block:: python

   result = experiment.run(Measure.QEXT, Measure.QSCA)
   dataframe = result.as_dataframe()
   # Columns: source:wavelength, scatterer:diameter, Qext, Qsca

See the :ref:`experiment gallery <experiment_index>` for larger
examples.

Singleton parameter axes are dropped by default. Pass ``drop_unique_level=False``
to ``run`` to retain them. Use ``result["Qsca"].isel({"source:wavelength": 0})``
to select a position, and inspect ``dims``, ``coords``, ``coordinate_units`` and
``units`` before converting or plotting.
