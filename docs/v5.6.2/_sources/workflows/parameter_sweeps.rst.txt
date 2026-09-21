.. _workflow_parameter_sweeps:

Parameter sweeps
================

Use :class:`PyMieSim.Experiment` and set classes when several values should be
evaluated as a parameter grid.

.. code-block:: python

   import numpy as np

   from PyMieSim import (
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
   result = experiment.get("Qsca", "Qext")
   result.isel({"measure": 0}).plot(x="source:wavelength", y="Qsca")

The result is a unit-aware PyMieSim LabeledArray.  Use ``.as_numpy()`` when a
raw NumPy array is preferable, ``.as_dataframe()`` for tabular interoperability,
or ``get_sequential`` for aligned sequential configurations.

Extracting NumPy and pandas data
---------------------------------

The experiment result remains labeled by default. Convert it explicitly when
using libraries that expect a raw array or a pandas table:

.. code-block:: python

   result = experiment.get("Qsca")

   values = result.as_numpy()
   print(values.shape)

   dataframe = result.as_dataframe()
   print(dataframe[["source:wavelength", "scatterer:diameter", "Qsca"]])

For multiple measures, the labeled array contains a ``measure`` dimension and
``as_dataframe()`` creates one output column per measure:

.. code-block:: python

   result = experiment.get("Qext", "Qsca")
   dataframe = result.as_dataframe()
   # Columns: source:wavelength, scatterer:diameter, Qext, Qsca

See the :ref:`experiment gallery <sphx_glr_gallery_experiment>` for larger
examples.
