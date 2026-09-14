.. _workflow_visualizing_fields:

Visualizing fields
==================

The setup facade exposes the existing plotting and representation helpers.
First create a simulation, then request the representation that matches the
question you want to investigate.

.. code-block:: python

   from PyMieSim import (
       Gaussian,
       PolarizationState,
       Simulation,
       Sphere,
       ureg,
   )

   simulation = Simulation(
       scatterer=Sphere(
           diameter=500 * ureg.nanometer,
           material=1.5,
           medium=1.0,
       ),
       source=Gaussian(
           wavelength=700 * ureg.nanometer,
           polarization=PolarizationState(angle=0 * ureg.degree),
           optical_power=1e-3 * ureg.watt,
           numerical_aperture=0.2,
       ),
   )

   simulation.plot_system()
   farfields = simulation.get_representation("farfields")
   farfields.plot()

For structured near-fields, far-fields, Stokes parameters, and footprints,
see the :ref:`single-particle gallery <sphx_glr_gallery_single>`.
