.. _experiment_index:

Experiment Module
=================

The ``experiment`` package is used to build large parameter sweeps.  It lets you
combine sources, scatterers and detectors to explore how scattering quantities
change across many configurations.

Reading a sweep plot
--------------------

Examples request observables with ``Measure`` members and print the actual
``result.dims``, ``result.shape``, and measure units before plotting. Each
remaining parameter dimension selects a family of curves; singleton axes are
removed by default. A shaded standard-deviation band summarizes variation
across the named sweep axis, not measurement uncertainty or a confidence interval.

``Q`` measures are dimensionless efficiencies, ``C`` measures carry area units,
and ``coupling`` is collected power in watts. Resonances can move with diameter,
wavelength, refractive index, and the surrounding medium. Compare curves with
one parameter held fixed before attributing a peak shift to a material change.
See :ref:`reading_results` for indexing and export, and :doc:`/measures` for
normalization conventions.
