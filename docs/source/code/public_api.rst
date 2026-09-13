Stable Python API
=================

The stable Python-facing API is centered on :class:`PyMieSim.Simulation` for
single-particle calculations and :class:`PyMieSim.Experiment` for parameter
sweeps. Both accept the historical measure strings and the discoverable
:class:`PyMieSim.Measure` enum::

    from PyMieSim import Measure, Simulation

    simulation.run(Measure.QSCA)

Use ``available_measures`` to inspect what the configured geometry and detector
support. Existing string calls remain supported for compatibility.

Typed results
-------------

Typed result containers are opt-in::

    result = simulation.run(Measure.QSCA, as_result=True)
    print(result.measure, result.quantity, result.units)

Both interfaces provide ``run`` and its synonym ``get``. One requested measure
with ``as_result=True`` returns ``SimulationResult``; multiple measures return
``SimulationResults``, a mapping keyed by measure name::

    results = experiment.run(Measure.QSCA, Measure.CSCA, as_result=True)
    scattering = results["Csca"].to("nanometer ** 2")
    print(scattering.dims, scattering.coords, scattering.coordinate_units)
    labeled = scattering.as_labeled_array()
    dataframe = labeled.as_dataframe()

Experiment results retain parameter dimensions, coordinates, and coordinate
units through unit conversion. Use ``drop_unique_level=False`` to retain axes
containing only one value. Single-simulation results convert to scalar labeled
arrays with no parameter dimensions.

Default return values remain lightweight: a single simulation returns a quantity
for one measure or a dictionary of quantities for multiple measures. Experiments
return ``LabeledArray``, with ``as_numpy()`` and ``as_dataframe()`` conversions.
Measure order is preserved. Empty requests, unsupported names, and duplicate
measures raise ``ValueError`` consistently before computation.

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
