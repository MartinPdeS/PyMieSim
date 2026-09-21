.. _particle_size_distributions:

Particle-size distributions: independent scattering
===================================================

PyMieSim can average the response of homogeneous spheres over a **number-based
particle-size distribution**, using discrete weights or quadrature for continuous distributions.

.. warning::

   This feature uses the **non-interacting (independent-scattering)
   approximation only**. Each particle responds to the prescribed incident
   illumination in isolation. It excludes electromagnetic coupling between
   particles, multiple scattering, positional correlations, and interference
   between fields scattered by different particles.

   It does not predict the response of a dense suspension, an aggregate, or an
   optically thick sample. Low concentration alone does not guarantee that
   multiple scattering is negligible: propagation distance also matters.
   PyMieSim does not infer or verify whether your sample meets these assumptions.

For a better approximation when particle correlations matter, refer to
`PackLab's correlation-based dependent-scattering calculations
<https://martinpdes.github.io/PackLab/docs/latest/scattering.html>`_. These
incorporate interparticle correlations; they are not a general solver for full
electromagnetic multiple scattering.

Averaging uses the independent-scattering approximation without a mode argument.
Results record this assumption in their metadata, including when converted to
typed results or dataframes.

Discrete sizes
--------------

Supply particle counts or number fractions at each diameter::

    from PyMieSim import (
        Experiment, ParticleSizeDistribution, PlaneWaveSet,
        PolarizationSet, SphereSet, ureg,
    )

    distribution = ParticleSizeDistribution(
        diameters=[100, 200, 400] * ureg.nanometer,
        number_weights=[1, 2, 1],
    )
    source = PlaneWaveSet(
        wavelength=[450, 550, 650] * ureg.nanometer,
        polarization=PolarizationSet(angles=0 * ureg.degree),
        amplitude=[1] * ureg.volt / ureg.meter,
    )
    spheres = SphereSet(
        diameter=distribution.diameters,
        material=[1.5],
        medium=[1.0],
    )
    experiment = Experiment(scatterer_set=spheres, source_set=source)
    result = experiment.average_size_distribution(
        distribution, "Csca", "Cext", "g",
        as_result=True,
    )
    scattering = result["Csca"].to("nanometer ** 2")
    table = scattering.as_labeled_array().as_dataframe()

Weights are normalized to sum to one. They must be **number weights**, not
volume fractions, mass fractions, or intensity-weighted size estimates.
For histogram measurements, use the number of particles in each bin. Samples
of a probability density require integration weights (including bin widths);
supplying only density values generally gives the wrong average.

The experiment diameter nodes must match ``distribution.diameters`` in the
same order. The method removes only the diameter axis and optionally singleton
axes; wavelength, material, medium, and detector sweeps remain separate.
It applies the same size distribution at every remaining grid point. It does
not represent a joint distribution of size and material. Sequential parameter
sets, cylinders, and core-shell particles are currently rejected.

Lognormal sizes
---------------

A lognormal number distribution can be constructed with::

    distribution = ParticleSizeDistribution.lognormal(
        median_diameter=200 * ureg.nanometer,
        geometric_std=1.3,
        sampling=32,
    )

Build the ``SphereSet`` using these new ``distribution.diameters`` before
averaging. The median is the **number median diameter**, and ``geometric_std``
is dimensionless: the standard deviation of log diameter is its natural
logarithm. A geometric standard deviation of one produces a single diameter.

Gauss-Hermite quadrature samples the full lognormal distribution without an
explicit cutoff. Nodes and weights are integration points, not histogram bins.
Increase ``sampling`` and compare the optical results to check convergence;
broad distributions and sharp Mie resonances can require more points. Tail
nodes can extend far beyond an intended particle-size range and increase
computation cost; use a discrete distribution for explicitly bounded
populations.

More distribution families
--------------------------

All constructors return the same ``ParticleSizeDistribution`` type, so use
``distribution.diameters`` in the ``SphereSet`` and call
``experiment.average_size_distribution(distribution, "Csca")`` as above.
All weights represent particle number. The non-interacting approximation and
PackLab guidance at the top of this page apply to every family.

.. list-table:: Available constructors
   :header-rows: 1
   :widths: 25 45 30

   * - Constructor
     - Parameters
     - Interpretation
   * - ``monodisperse``
     - ``diameter``
     - All particles have one diameter.
   * - ``uniform``
     - ``minimum_diameter``, ``maximum_diameter``, ``sampling=32``
     - Constant number density per unit diameter within the bounds.
   * - ``lognormal``
     - ``median_diameter``, ``geometric_std``, ``sampling=32``
     - Normal distribution of log diameter.
   * - ``truncated_normal``
     - ``mean_diameter``, ``standard_deviation``, required positive bounds,
       ``sampling=64``
     - Normal density conditioned on the bounded interval.
   * - ``triangular``
     - ``minimum_diameter``, ``mode_diameter``, ``maximum_diameter``,
       ``sampling=16``
     - Linear density on either side of a peak.
   * - ``mixture``
     - ``distributions``, ``number_weights``
     - Combine populations by their particle counts or number fractions.

For example::

    single_size = ParticleSizeDistribution.monodisperse(200 * ureg.nanometer)
    uniform = ParticleSizeDistribution.uniform(
        100 * ureg.nanometer, 300 * ureg.nanometer,
    )
    normal = ParticleSizeDistribution.truncated_normal(
        mean_diameter=200 * ureg.nanometer,
        standard_deviation=40 * ureg.nanometer,
        minimum_diameter=100 * ureg.nanometer,
        maximum_diameter=300 * ureg.nanometer,
    )
    triangular = ParticleSizeDistribution.triangular(
        minimum_diameter=100 * ureg.nanometer,
        mode_diameter=180 * ureg.nanometer,
        maximum_diameter=300 * ureg.nanometer,
    )
    bimodal = ParticleSizeDistribution.mixture(
        distributions=[single_size, uniform],
        number_weights=[3, 1],
    )

The normal parameters describe the **untruncated** distribution. Conditioning
on the specified bounds generally changes its actual mean and standard
deviation. Bounds are explicit so negative particle sizes are never generated.
Use ``monodisperse`` when the width is zero.

Uniform and truncated-normal constructors use Gauss-Legendre quadrature.
Triangular distributions integrate each nonempty side of the peak separately:
``sampling`` is the point count **per side**, so an interior peak produces twice
that many nodes. An endpoint peak is supported. Always compare optical results
at increasing sampling, particularly for a normal distribution on an interval
much wider than its standard deviation or when Mie resonances are sharp.

Mixtures preserve each component's quadrature nodes. In the example, 75% of
particles belong to ``single_size`` and 25% to ``uniform``, regardless of how
many nodes either component has. Duplicate diameters are retained. These are
mixtures of sizes, not a joint distribution of size and composition.

What the average means
----------------------

Let :math:`w_i` be normalized number fractions, :math:`d_i` diameters, and
:math:`A_i=\pi d_i^2/4` projected areas. The supported reductions are:

.. math::

   \overline{C} = \sum_i w_i C_i,
   \qquad
   Q_{\mathrm{ensemble}} =
   \frac{\sum_i w_i C_i}{\sum_i w_i A_i},
   \qquad
   g_{\mathrm{ensemble}} =
   \frac{\sum_i w_i C_{\mathrm{sca},i} g_i}
        {\sum_i w_i C_{\mathrm{sca},i}}.

Cross sections ``Csca``, ``Cext``, ``Cabs``, ``Cback``, ``Cforward``, and
``Cpr`` are means **per particle**, in area units. Their corresponding
``Q`` measures divide by mean projected area; they are not arithmetic number
averages of individual efficiencies. Asymmetry ``g`` is weighted by scattered
power through ``Csca``; it is undefined when the mean scattering cross section
is zero, in which case the method raises ``ValueError``.

``coupling`` is the number-averaged **single-particle detector power**, in
watts, for identical illumination and detector geometry. Even for a coherent
mode detector, this averages individual powers, never complex amplitudes.
There is no integration over particle positions, illumination gradients, or
relative optical phases. It is not the total detected power of a suspension.
Complex multipole coefficients, ``Qratio``, and ``Cratio`` are not supported.

A sample's particle count, concentration, and path length are not supplied or
inferred. Multiplying a mean cross section by number density gives an
independent-scattering coefficient; predicting transmitted or detected light
through a sample requires a suitable propagation model and its own validity
assessment.

The result's ``attrs`` (or ``metadata`` for typed results) records the
approximation, assumptions, weighting rule, diameter nodes, and number
fractions. These records survive unit conversion and labeled-array/dataframe
export.

For the cross-section integration and scattering-weighted asymmetry
conventions, see `Geer et al. (2021), Section 2
<https://gmd.copernicus.org/articles/14/7497/2021/gmd-14-7497-2021.html#section2>`_.

API
---

.. autoclass:: PyMieSim.ParticleSizeDistribution
   :members:

.. automethod:: PyMieSim.Experiment.average_size_distribution
