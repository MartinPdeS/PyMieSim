Getting started in 10 minutes
=============================

This tutorial goes from installation to a parameter sweep and a detector
measurement. All quantities with physical dimensions use PyMieSim's unit
registry, so changing units does not change the calculation.

Install
-------

Install the released package with pip::

   python -m pip install PyMieSim

To run tests or build documentation, install the optional dependencies::

   python -m pip install "PyMieSim[testing,documentation]"

Verify the installation with ``python -c "import PyMieSim; print(PyMieSim.__version__)"``.
See :doc:`compatibility` for wheel targets and :doc:`troubleshooting` for source
builds. The examples below use constant refractive indices and need no material
catalog download. Named materials require ``python -m PyOptik setup`` first.

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

``qsca`` is a scalar, dimensionless quantity (about 0.208 for this setup).
It is the scattering cross section divided by the sphere's projected area;
it is not a fraction of the source power.

Inspect the available measures before running a larger study::

   print(simulation.available_measures)

For explicit result metadata, opt in to the typed result container::

   result = simulation.run(Measure.QSCA, as_result=True)
   print(result.measure, result.quantity, result.units)

Parameter sweep
---------------

Use ``Experiment`` when several parameter dimensions must be evaluated. Scalar
parameters are broadcast over the grid, and the default result is a native
``LabeledArray`` with named dimensions, coordinates, and units.

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
   result = experiment.get(Measure.QSCA)
   result.plot(x="scatterer:diameter", y="Qsca")

With five wavelengths and nine diameters, the result has shape ``(5, 9)``
and dimensions ``("source:wavelength", "scatterer:diameter")``. Each plotted
curve holds one wavelength fixed. Values are dimensionless; coordinates retain
their own units. See :ref:`reading_results` for selection and export.

Detector coupling
-----------------

Coupling requires a detector. A photodiode is the simplest detector to add;
coherent modes and integrating spheres are covered in the
:doc:`workflows/detector_coupling` guide. Reuse the source and sphere above:

.. code-block:: python

   from PyMieSim import Photodiode

   detector = Photodiode(
       numerical_aperture=0.2,
       sampling=500,
       phi_offset=0 * ureg.degree,
       gamma_offset=0 * ureg.degree,
       medium=1.0,
   )
   detected_simulation = Simulation(
       source=source,
       scatterer=sphere,
       detector=detector,
   )
   collected_power = detected_simulation.run(Measure.COUPLING)
   print(collected_power.to("nanowatt"))

The output is a positive scalar power, expressed here in nanowatts. It describes
light collected over this detector's acceptance cone. Changing the numerical
aperture changes that cone; increasing ``sampling`` checks angular integration
convergence. This power and ``Qsca`` answer different measurement questions.

Next steps
----------

* :doc:`workflows/parameter_sweeps` for larger grids and DataFrame operations.
* :doc:`measures` for definitions, units, and geometry support.
* :doc:`performance` for runtime and memory guidance.
* :doc:`troubleshooting` for installation and numerical issues.
* :doc:`examples` for runnable single-particle, experiment, and validation examples.
