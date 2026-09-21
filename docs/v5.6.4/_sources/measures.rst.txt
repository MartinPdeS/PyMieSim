Measure reference
=================

Measures are available through ``PyMieSim.Measure`` or their historical string
names. For example, ``Measure.QSCA`` and ``"Qsca"`` are equivalent.

Efficiencies and dimensionless quantities
-----------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 35 18 25

   * - Name
     - Meaning
     - Unit
     - Geometries
   * - ``Qsca``
     - Scattering efficiency
     - dimensionless
     - sphere, cylinder, core--shell
   * - ``Qext``
     - Extinction efficiency
     - dimensionless
     - sphere, cylinder, core--shell
   * - ``Qabs``
     - Absorption efficiency
     - dimensionless
     - sphere, cylinder, core--shell
   * - ``Qback``
     - Back-scattering efficiency
     - dimensionless
     - sphere, cylinder, core--shell
   * - ``Qforward``
     - Forward-scattering efficiency
     - dimensionless
     - sphere, cylinder, core--shell
   * - ``Qratio``
     - Back-scattering-to-scattering ratio
     - dimensionless
     - sphere, cylinder, core--shell
   * - ``Qpr``
     - Radiation-pressure efficiency
     - dimensionless
     - sphere, cylinder, core--shell
   * - ``g``
     - Anisotropy factor
     - dimensionless
     - sphere, cylinder, core--shell
   * - ``g_with_farfields``
     - Anisotropy factor evaluated from sampled far fields
     - dimensionless
     - single-particle API
   * - ``size_parameter``
     - Dimensionless size parameter
     - dimensionless
     - sphere, cylinder, core--shell

Cross sections
--------------

``Csca``, ``Cext``, ``Cabs``, ``Cback``, ``Cforward``, ``Cratio``, ``Cpr``,
and ``cross_section`` are area-valued measures in square metres. Results are
returned as unit-aware quantities in single simulations and as numerical
columns with unit metadata in experiment DataFrames.

Coefficient measures
--------------------

The experiment API also exposes the first few Mie coefficients (``a1``,
``a2``, ``a3``, ``b1``, ``b2``, and ``b3``) where supported by the configured
scatterer. These are dimensionless complex coefficients and are primarily
useful for diagnostics and theory comparisons.

Detector measures
-----------------

``coupling`` is collected optical power and has units of watts. It is only
available when a detector is configured. The detector type, angular sampling,
polarization filter, and coherent/non-coherent mode determine its exact
meaning. See :doc:`workflows/detector_coupling`.

Discoverability
---------------

Use the configured object rather than guessing which measures are valid::

   simulation.available_measures
   experiment.available_measures

An invalid request produces an error listing the available names for the
current geometry and detector configuration.
