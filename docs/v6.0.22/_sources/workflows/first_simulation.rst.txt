.. _workflow_first_simulation:

Your first simulation
=====================

Use :class:`PyMieSim.Simulation` for one source, one scatterer, and an
optional detector.  This is the simplest way to calculate scattering
properties.

.. code-block:: python

   from PyMieSim import (
       Measure,
       Gaussian,
       PolarizationState,
       Simulation,
       Sphere,
       ureg,
   )

   source = Gaussian(
       wavelength=750 * ureg.nanometer,
       polarization=PolarizationState(angle=0 * ureg.degree),
       optical_power=1e-3 * ureg.watt,
       numerical_aperture=0.2,
   )

   simulation = Simulation(
       scatterer=Sphere(
           diameter=200 * ureg.nanometer,
           material=1.5,
           medium=1.0,
       ),
       source=source,
   )

   result = simulation.run(Measure.QSCA, Measure.QEXT)
   print(result)

The result is a dictionary with ``"Qsca"`` and ``"Qext"`` keys, each holding
a scalar dimensionless quantity. For the real refractive index used here,
absorption vanishes, so scattering and extinction efficiencies agree.

Next: :ref:`workflow_parameter_sweeps`.
