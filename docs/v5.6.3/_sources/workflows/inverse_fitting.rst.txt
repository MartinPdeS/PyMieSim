.. _inverse_fitting:

Inverse parameter fitting
=========================

PyMieSim includes a small NumPy-only, derivative-free fitting layer. It adds no
required optimization dependency. The forward model remains a normal PyMieSim
``Simulation`` or ``Experiment``; the fitting layer only supplies parameters,
collects predictions, and compares them with observations.

``Parameter`` defines a named initial value and mandatory bounds. ``Observation``
holds measured values and optional positive standard uncertainties. The user
supplies a callable receiving a mapping of parameter names to values and
returning a prediction with the same shape as the observation.

The default coordinate search is deterministic and intended for a small number
of parameters. It is derivative-free and has no global-optimum guarantee; use
physical bounds and repeat fits from several initial guesses for multimodal
resonant problems.

Set ``show_progress=True`` to display the iteration, current objective, and
normalized search step while a fit is running. Progress output is disabled by
default, which keeps scripts and notebooks quiet.

``FitResult`` contains fitted parameters, predictions, residuals, objective
value, iteration and evaluation counts, and convergence status. Unit-bearing
parameters are passed to the model with their original units, and uncertainties
are used to weight residuals.

The fitted values are conditional on the forward model. A fit to a dense
suspension using isolated-particle Mie calculations can compensate for missing
multiple scattering or particle correlations. For correlation-aware scattering
models, refer to `PackLab's scattering workflow
<https://martinpdes.github.io/PackLab/docs/latest/scattering.html>`_.

Examples
--------

The examples below introduce the fitting interface from a small analytic model
to a unit-aware PyMieSim calculation. All use the same deterministic coordinate
search. It is best suited to a few bounded parameters; inspect the objective,
residuals, and fitted values before drawing physical conclusions.

1. :ref:`sphx_glr_gallery_inverse_fit_linear_scale.py` fits a scalar scale
   factor and shows the minimal API.
2. :ref:`sphx_glr_gallery_inverse_fit_unit_aware_diameter.py` fits a parameter
   expressed in nanometres while the optimizer works on its magnitude.
3. :ref:`sphx_glr_gallery_inverse_fit_sphere_scattering.py` recovers an
   isolated-sphere diameter from a scattering-efficiency observation.

The third example uses an isolated-particle forward model. Measurements from a
concentrated or correlated suspension require a forward model that represents
those interactions; refer to PackLab for those approximations.

.. autoclass:: PyMieSim.Parameter
   :members:
.. autoclass:: PyMieSim.Observation
   :members:
.. autoclass:: PyMieSim.FitResult
   :members:
.. autofunction:: PyMieSim.fit_parameters
