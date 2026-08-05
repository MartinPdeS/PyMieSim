Getting started in 10 minutes
=============================

This tutorial goes from installation to a parameter sweep and a detector
measurement. All quantities with physical dimensions use PyMieSim's unit
registry, so changing units does not change the calculation.

Install
-------

Install the released package with pip::

   python -m pip install PyMieSim

For development, install the testing and documentation extras as needed::

   python -m pip install "PyMieSim[testing,documentation]"

Single-particle calculation
---------------------------

Create a source, a scatterer, and a simulation. ``Measure`` members are
discoverable in an IDE; the historical string names such as ``"Qsca"`` remain
supported.

.. code-block:: python

   from PyMieSim import (
       Gaussian,
       Measure,
       PolarizationState,
       Simulation,
       Sphere,
       ureg,
   )

   source = Gaussian(
       wavelength=633 * ureg.nanometer,
       polarization=PolarizationState(angle=0 * ureg.degree),
       optical_power=1e-3 * ureg.watt,
       numerical_aperture=0.2,
   )
   sphere = Sphere(
       diameter=200 * ureg.nanometer,
       material=1.5 + 0.01j,
       medium=1.0,
   )

   simulation = Simulation(scatterer=sphere, source=source)
   qsca = simulation.run(Measure.QSCA)
   print(qsca)

Inspect the available measures before running a larger study::

   print(simulation.available_measures)

For explicit result metadata, opt in to the typed result container::

   result = simulation.run(Measure.QSCA, as_result=True)
   print(result.measure, result.quantity, result.units)

Parameter sweep
--------------

Use ``Experiment`` when several parameter dimensions must be evaluated. Scalar
parameters are broadcast over the grid, and the default result is a
unit-aware ``PyMieSimDataFrame``.

.. code-block:: python

   import numpy as np
   from PyMieSim import (
       Experiment,
       GaussianSet,
       PolarizationSet,
       SphereSet,
   )

   sweep_source = GaussianSet(
       wavelength=np.linspace(500, 700, 5) * ureg.nanometer,
       polarization=PolarizationSet(angles=0 * ureg.degree),
       optical_power=1e-3 * ureg.watt,
       numerical_aperture=0.2,
   )
   sweep_sphere = SphereSet(
       diameter=np.linspace(100, 500, 9) * ureg.nanometer,
       material=1.5,
       medium=1.0,
   )

   experiment = Experiment(
       scatterer_set=sweep_sphere,
       source_set=sweep_source,
   )
   dataframe = experiment.get(Measure.QSCA)
   dataframe.plot_standard(x="scatterer:diameter", y="Qsca")

Detector coupling
-----------------

Coupling requires a detector. A photodiode is the simplest detector to add;
coherent modes and integrating spheres are covered in the
:doc:`workflows/detector_coupling` guide.

Next steps
----------

* :doc:`workflows/parameter_sweeps` for larger grids and DataFrame operations.
* :doc:`measures` for definitions, units, and geometry support.
* :doc:`performance` for runtime and memory guidance.
* :doc:`troubleshooting` for installation and numerical issues.
* :doc:`examples` for runnable single-particle, experiment, and validation examples.
