.. _single_index:

Single Module
=============

The ``single`` package focuses on scattering from a single particle.  It exposes
simple classes for defining a source, a scatterer and a detector so that the
interaction can be studied in detail.

Key components
--------------
- **Sources** – plane waves and Gaussian beams.
- **Scatterers** – spheres, cylinders and core--shell particles.
- **Detectors** – photodiodes as well as coherent and incoherent detectors.

Example
-------
This example computes the power coupled into a photodiode from a tiny sphere.

.. code-block:: python

   from PyMieSim.single.scatterer import Sphere
   from PyMieSim.single.source import Gaussian
   from PyMieSim.single.detector import Photodiode
   from PyMieSim.units import nanometer, degree, watt, AU, RIU

   source = Gaussian(
       wavelength=450 * nanometer,
       polarization=0 * degree,
       optical_power=1 * watt,
       NA=0.3 * AU,
   )

   scatterer = Sphere(
       diameter=6 * nanometer,
       source=source,
       medium_property=1.0 * RIU,
       property=1.4 * RIU,
   )

   detector = Photodiode(
       NA=0.1 * AU,
       phi_offset=0 * degree,
       gamma_offset=0 * degree,
       sampling=200 * AU,
       polarization_filter=None,
   )

   coupling = detector.coupling(scatterer)
   print(coupling)
