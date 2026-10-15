.. _workflow_detector_coupling:

Detector coupling
==================

Add a detector to a :class:`PyMieSim.Simulation` when the quantity of interest
is collected or coupled power rather than only a scatterer property.

.. code-block:: python

   from PyMieSim import (
       Measure,
       Gaussian,
       Photodiode,
       PolarizationState,
       Simulation,
       Sphere,
       ureg,
   )

   source = Gaussian(
       wavelength=1000 * ureg.nanometer,
       polarization=PolarizationState(angle=30 * ureg.degree),
       optical_power=1 * ureg.watt,
       numerical_aperture=0.3,
   )

   detector = Photodiode(
       sampling=500,
       numerical_aperture=0.2,
       phi_offset=0 * ureg.degree,
       gamma_offset=0 * ureg.degree,
       medium=1.0,
   )

   simulation = Simulation(
       scatterer=Sphere(
           diameter=300 * ureg.nanometer,
           material=1.5,
           medium=1.0,
       ),
       source=source,
       detector=detector,
   )

   coupling = simulation.run(Measure.COUPLING)
   print(coupling)

The output is a scalar quantity in watts. It is power collected within the
detector cone for the specified source power. Repeat with larger ``sampling``
to check angular integration convergence before comparing small differences.

For angularly varying detector parameters, use ``PhotodiodeSet`` and
``Experiment``.  See the :ref:`coupling examples
<sphx_glr_gallery_single_coupling.py>` and experiment gallery.
