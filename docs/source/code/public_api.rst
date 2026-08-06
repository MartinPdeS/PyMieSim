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

For multiple single-particle measures, ``as_result=True`` returns a mapping of
measure names to ``SimulationResult`` objects. For experiments it returns an
``ExperimentResult`` whose original DataFrame is available as
``result.result`` or through ``result.to_pandas()``.

Advanced and legacy access
---------------------------

The ``Simulation`` facade is the recommended stable entry point. The compiled
backend remains available through ``simulation.advanced`` (and the historical
``simulation.setup`` alias) for specialized methods and migration of existing
applications. Backend methods reached through attribute forwarding are not
part of the stable API contract.
