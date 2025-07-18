.. _experiment_index:

Experiment Module
=================

The ``experiment`` package is used to build large parameter sweeps.  It lets you
combine sources, scatterers and detectors to explore how scattering quantities
change across many configurations.

Key components
--------------
- **Sources** – plane waves and Gaussian beams with user-defined wavelength and
  polarization.
- **Scatterers** – spherical, cylindrical and core--shell geometries.
- **Detectors** – photodiodes, coherent modes and integrating spheres.

Example
-------
The snippet below retrieves the scattering efficiency for a range of refractive
indices.

.. code-block:: python

   from PyMieSim.experiment.scatterer import Sphere
   from PyMieSim.experiment.source import Gaussian
   from PyMieSim.experiment import Setup
   from PyMieSim.units import nanometer, degree, watt, AU, RIU
   import numpy as np

   source = Gaussian(
       wavelength=[500., 1000., 1500.] * nanometer,
       polarization=30. * degree,
       optical_power=1e-3 * watt,
       NA=0.2 * AU,
   )

   scatterer = Sphere(
       diameter=800. * nanometer,
       property=np.linspace(1.3, 1.9, 150) * RIU,
       medium_property=1. * RIU,
       source=source,
   )

   experiment = Setup(scatterer=scatterer, source=source)
   df = experiment.get('Qsca')
   df.plot_data(x='scatterer:property')
