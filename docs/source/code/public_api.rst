Stable Python API
=================

.. _choosing_an_interface:

Choose an interface
-------------------

.. list-table:: Calculation interfaces
   :header-rows: 1
   :widths: 18 28 28 26

   * - Interface
     - Use it for
     - Inputs
     - Default output
   * - ``Simulation``
     - One particle and illumination
     - Source, scatterer, optional detector
     - ``SimulationResult`` or ``SimulationResults``
   * - ``Experiment``
     - Cartesian grids or aligned sequential configurations
     - Source, scatterer, optional detector sets
     - ``SimulationResult`` or ``SimulationResults``
   * - ``fit_parameters``
     - Estimate parameters from measurements
     - Model callable, ``MeasuredData``, bounded ``FitParameter`` objects
     - ``ParameterFitResult``

These minimal examples share a constant-index sphere and a plane wave:

.. code-block:: python

    from PyMieSim import (
        Experiment, FitParameter, MeasuredData, Measure, PlaneWave,
        PlaneWaveSet, PolarizationSet, PolarizationState, Simulation,
        Sphere, SphereSet, fit_parameters, ureg,
    )

    source = PlaneWave(
        wavelength=532 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        amplitude=1 * ureg.volt / ureg.meter,
    )
    simulation = Simulation(
        source=source,
        scatterer=Sphere(diameter=180 * ureg.nanometer, material=1.59, medium=1.0),
    )
    qsca = simulation.run(Measure.QSCA).quantity  # Scalar dimensionless quantity.

    experiment = Experiment(
        source_set=PlaneWaveSet(
            wavelength=[450, 550, 650] * ureg.nanometer,
            polarization=PolarizationSet(angles=0 * ureg.degree),
            amplitude=1 * ureg.volt / ureg.meter,
        ),
        scatterer_set=SphereSet(
            diameter=[150, 180] * ureg.nanometer,
            material=1.59,
            medium=1.0,
        ),
    )
    grid = experiment.run(Measure.QSCA)  # Shape (3, 2), dimensionless values.

    def model(parameters):
        return Simulation(
            source=source,
            scatterer=Sphere(
                diameter=parameters["diameter"], material=1.59, medium=1.0,
            ),
        ).run(Measure.QSCA).quantity

    fit = fit_parameters(
        model=model,
        observation=MeasuredData(values=qsca),
        parameters=[FitParameter(
            name="diameter",
            initial=150 * ureg.nanometer,
            bounds=(100 * ureg.nanometer, 220 * ureg.nanometer),
        )],
    )
    print(fit.parameters["diameter"])  # Approximately 180 nanometers.

Use ``available_measures`` to inspect what the configured geometry and detector
support. ``Measure.QSCA`` and the historical string ``"Qsca"`` are equivalent.
Measure arguments are positional; constructor
and fitting arguments use keywords in these examples.

API migration
-------------

Use ``Simulation`` or ``Experiment`` and call ``run``. The former public
``Setup`` entry points, ``get`` alias, ``get_representation`` string dispatcher,
and ``as_result`` option on ``run`` have been removed. For example, replace
``simulation.get("Qsca")`` with ``simulation.run("Qsca").quantity`` when a
Pint quantity is needed, and replace ``simulation.get_representation("spf")``
with ``simulation.spf()``. For a sweep, select a result by measure name before
plotting it, or call ``as_labeled_array()`` for the combined measure axis.

.. _reading_results:

Read, select, and export results
--------------------------------

Both ``Simulation.run`` and ``Experiment.run`` return a ``SimulationResult``
for one measure, or a ``SimulationResults`` mapping for several measures.
A result exposes its Pint ``quantity``, ``units``, ``dims``, ``coords`` and
``coordinate_units``. Experiment values carry a Cartesian parameter grid.

.. code-block:: python

    scalar = simulation.run(Measure.CSCA)
    print(scalar.to("nanometer ** 2").quantity)

    results = experiment.run(Measure.QSCA, Measure.CSCA)
    grid = results["Qsca"]
    print(grid.dims, grid.shape, grid.units, grid.coordinate_units)
    first_wavelength = grid.isel(indexers={"source:wavelength": 0})
    first_wavelength.plot(x="scatterer:diameter")
    values = grid.as_numpy()
    table = results.as_dataframe()
    scattering = results["Csca"].to("nanometer ** 2")
    labeled = scattering.as_labeled_array()

Singleton dimensions are removed by default; pass ``drop_unique_level=False``
to ``Experiment.run`` to retain them. Multiple measures are selected by name;
``as_labeled_array()`` on the collection combines them along a ``measure``
axis. ``as_numpy()`` returns magnitudes in each measure's current units.
CSV does not preserve Pint units or dataframe attributes, so record units and
coordinate units separately. Single simulations have no parameter dimensions.

Empty requests, unsupported names, and duplicate measures raise ``ValueError``
before computation. For the meaning of each output, see :doc:`../measures`.

Typing
------

The simulation interface provides overloads for single and multiple measures,
typed results, and structured versus angle-sampled fields. Native setup and
labeled-array interfaces have accompanying ``.pyi`` files. Run ``make typecheck``
with development dependencies installed to check the Python API. After building
the native extensions, run ``make stubcheck`` to compare the stubs with the
installed bindings. Both checks run in CI.

Particle-size distributions
---------------------------

``ParticleSizeDistribution`` supports discrete, monodisperse, uniform, lognormal,
truncated-normal, triangular, and mixed number-based populations.
``Experiment.average_size_distribution`` averages homogeneous-sphere sweeps
using the non-interacting approximation, as specified in its docstring.
It excludes multiple scattering and interparticle interference. For a better
approximation when particle correlations matter, refer to
`PackLab <https://martinpdes.github.io/PackLab/docs/latest/scattering.html>`_.
See :ref:`particle_size_distributions` for physical weighting rules,
limitations, and an example.

Inverse fitting
---------------

The NumPy-only ``fit_parameters`` API fits bounded physical parameters against
unit-aware observations without adding SciPy as a dependency. See
:ref:`inverse_fitting` for the forward-model adapter pattern and limitations.
