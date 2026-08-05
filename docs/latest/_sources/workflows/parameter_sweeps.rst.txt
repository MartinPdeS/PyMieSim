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
   dataframe = experiment.get("Qsca", "Qext")
   dataframe.plot(x="source:wavelength")

The result is a unit-aware PyMieSim DataFrame.  Use ``as_numpy=True`` when a
raw NumPy array is preferable, or ``get_sequential`` for aligned sequential
configurations.

See the :ref:`experiment gallery <sphx_glr_gallery_experiment>` for larger
examples.
