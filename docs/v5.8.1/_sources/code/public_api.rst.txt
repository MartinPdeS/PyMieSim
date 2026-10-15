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
     - Quantity or dictionary of quantities
   * - ``Experiment``
     - Cartesian grids or aligned sequential configurations
     - Source, scatterer, optional detector sets
     - ``LabeledArray``
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
    qsca = simulation.run(Measure.QSCA)  # Scalar dimensionless quantity.

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
        ).run(Measure.QSCA)

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
``run`` and ``get`` are synonyms. Measure arguments are positional; constructor
and fitting arguments use keywords in these examples.

.. _reading_results:

Read, select, and export results
--------------------------------

.. list-table:: Return types
   :header-rows: 1
   :widths: 28 22 25 25

   * - Request
     - Return type
     - Selection
     - Units and export
   * - One simulation measure
     - Pint quantity
     - Scalar magnitude
     - ``.to("nanometer ** 2")`` for area measures
   * - Multiple simulation measures
     - Dictionary of quantities
     - ``values["Csca"]``
     - Convert each quantity separately
   * - Experiment, one or multiple measures
     - ``LabeledArray``
     - ``.isel(indexers={...})``
     - Units in ``attrs``; ``.as_numpy()`` or ``.as_dataframe()``
   * - One measure with ``as_result=True``
     - ``SimulationResult``
     - ``.quantity``, ``.dims``, ``.coords``
     - ``.to(...)`` and ``.as_labeled_array()``
   * - Multiple measures with ``as_result=True``
     - ``SimulationResults``
     - ``results["Csca"]`` gives a ``SimulationResult``
     - Convert or export each named result

Continue with ``simulation`` and ``experiment`` from the example above:

.. code-block:: python

    scalar = simulation.run(Measure.CSCA)
    print(scalar.to("nanometer ** 2"))

    values = simulation.run(Measure.QSCA, Measure.CSCA)
    print(values["Qsca"], values["Csca"])

    grid = experiment.run(Measure.QSCA)
    print(grid.dims, grid.shape)  # Wavelength then diameter: (3, 2).
    first_wavelength = grid.isel(indexers={"source:wavelength": 0})
    print(first_wavelength.shape)  # (2,), one value per diameter.
    print(grid.attrs["units"], grid.attrs["coordinate_units"])
    magnitudes = grid.as_numpy()  # Shape (3, 2); plain numbers carry no units.
    table = grid.as_dataframe()
    table.to_csv(path_or_buf="scattering.csv", index=False)

    typed = experiment.run(Measure.CSCA, as_result=True)
    converted = typed.to("nanometer ** 2")
    print(converted.units, converted.dims, converted.coordinate_units)
    table = converted.as_labeled_array().as_dataframe()

    results = experiment.run(Measure.QSCA, Measure.CSCA, as_result=True)
    scattering = results["Csca"].to("nanometer ** 2")
    print(scattering.quantity)

Dimension names describe the swept parameters. Singleton dimensions are removed
by default; pass ``drop_unique_level=False`` to ``Experiment.run`` to retain
them. Multiple default experiment measures add a ``measure`` dimension;
``as_dataframe()`` makes one output column per measure. Typed collections keep
measures separate, with the same parameter coordinates on each result.

CSV stores numerical columns and labels, but does not preserve Pint objects
or dataframe attributes. Record ``attrs["units"]`` and
``attrs["coordinate_units"]`` alongside a CSV, or include units explicitly in
its column headings. Unit conversion on a typed result preserves coordinates.
Single-simulation typed results have no parameter dimensions.

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
