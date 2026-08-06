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
measure names to ``SimulationResult`` objects. Experiment sweeps always return
``LabeledArray``; use ``result.as_numpy()`` or ``result.as_dataframe()`` for
explicit conversion.
