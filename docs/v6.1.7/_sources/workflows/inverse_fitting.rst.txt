.. _inverse_fitting:

Inverse parameter fitting
=========================

PyMieSim includes a small NumPy-only, derivative-free fitting layer. It adds no
required optimization dependency. The forward model remains a normal PyMieSim
``Simulation`` or ``Experiment``; the fitting layer only supplies parameters,
collects predictions, and compares them with observations.

``FitParameter`` defines a named initial value and mandatory bounds.
``MeasuredData`` holds measured values and optional positive standard
uncertainties. The user supplies a callable receiving a mapping of parameter
names to values and returning a prediction with the same shape as the observation.

The earlier names ``Parameter``, ``Observation``, and ``FitResult`` remain
available as aliases for ``FitParameter``, ``MeasuredData``, and
``ParameterFitResult``, respectively.

The default coordinate search is deterministic and intended for a small number
of parameters. It is derivative-free and has no global-optimum guarantee; use
physical bounds and repeat fits from several initial guesses for multimodal
resonant problems.

Set ``show_progress=True`` to display the iteration, current objective, and
normalized search step while a fit is running. Progress output is disabled by
default, which keeps scripts and notebooks quiet.

``ParameterFitResult`` contains fitted parameters, predictions, residuals,
objective value, iteration and evaluation counts, and convergence status.
Unit-bearing parameters are passed to the model with their original units,
and uncertainties are used to weight residuals.

Quantity bounds are converted to the initial parameter's units. Quantity
predictions and uncertainties are converted to the observation's units before
computing residuals; incompatible dimensions raise a Pint
``DimensionalityError``. Bare numbers are interpreted as magnitudes in the
corresponding reference units. When the reference is a bare number, quantities
must be dimensionless. Zero-valued initial quantities are supported.

``ParameterFitResult.prediction`` and ``ParameterFitResult.observed`` are numerical
arrays in the observation's units. ``ParameterFitResult.residuals`` contains
their numerical difference, divided by the uncertainty when one is supplied.

The fitted values are conditional on the forward model. A fit to a dense
suspension using isolated-particle Mie calculations can compensate for missing
multiple scattering or particle correlations. For correlation-aware scattering
models, refer to `PackLab's scattering workflow
<https://martinpdes.github.io/PackLab/docs/latest/scattering.html>`_.

Fit a noisy scattering spectrum
-------------------------------

The complete script below estimates diameter and refractive index from nine
synthetic wavelength channels. It starts with a known sphere, adds Gaussian
noise with a fixed seed, then fits a second sphere model. Wavelengths and
diameter carry length units; ``Qsca`` and its uncertainties are dimensionless.
No material catalog or external measurements are needed.

1. ``scattering_spectrum`` adapts single-particle simulations to an array of
   predictions with shape ``(9,)``.
2. ``MeasuredData`` pairs those observations with standard uncertainties of
   the same shape. One deliberately biased channel is assigned a larger
   uncertainty to demonstrate weighting.
3. ``FitParameter`` supplies physical bounds and initial guesses. The model
   receives diameter as a quantity and refractive index as a plain number.
4. The overlay compares predictions with measurements; the lower panel exposes
   systematic differences that a single objective value can hide.

.. literalinclude:: ../../examples/interpretation/inverse_spectrum_diagnostics.py
   :language: python
   :start-at: import matplotlib.pyplot

The minimized objective is the sum of squared standardized residuals:

.. math::

   r_i=\frac{f_i(\mathbf{p})-y_i}{\sigma_i},\qquad
   J(\mathbf{p})=\sum_i r_i^2.

Without uncertainties, residuals are numerical prediction-minus-observation
values in the observation's units. ``objective`` is the raw sum, not a mean,
reduced chi-square, or parameter uncertainty estimate. Each ``sigma`` must be
positive; correlations between measurement channels are not modeled.

Interpreting convergence and ambiguity
--------------------------------------

For this script, expect a fitted spectrum near the noisy observations and
parameters near the generating values (310 nm and 1.58), rather than exact
recovery. There are nine predictions and nine residuals. A positive residual
means the model lies above the observation; the deliberately biased channel
can retain a negative residual even after convergence.

``success=True`` means the coordinate search reached its step tolerance.
It does not establish a global optimum, a unique physical solution, or a good
model. Check ``message``, ``objective``, and ``evaluations`` as well as the plot.
If the fit stops at its iteration limit, inspect the trajectory and bounds
before increasing ``max_iterations``.

Repeat the fit from several diameter and refractive-index starting points
within the same bounds. Compare objectives and predictions: different parameter
pairs with similar objectives indicate ambiguity. A parameter close to a bound
can indicate an active physical constraint or an overly restrictive search
range. Plot an objective grid across the two parameters to inspect valleys and
competing minima; add informative wavelength channels when the data permit it.
The returned result does not provide confidence intervals or a covariance matrix.

A visually good fit remains conditional on the isolated-particle forward model.
Use residual structure to investigate calibration, material dispersion, or
missing physics before interpreting fitted values as material properties.
See the :ref:`rendered diagnostic example
<sphx_glr_gallery_interpretation_inverse_spectrum_diagnostics.py>` for its plots.

Smaller examples
----------------

* :ref:`sphx_glr_gallery_inverse_fit_linear_scale.py` introduces the minimal API.
* :ref:`sphx_glr_gallery_inverse_fit_unit_aware_diameter.py` demonstrates quantity parameters.
* :ref:`sphx_glr_gallery_inverse_fit_sphere_scattering.py` fits one sphere diameter.

API
---

.. autoclass:: PyMieSim.FitParameter
   :members:
.. autoclass:: PyMieSim.MeasuredData
   :members:
.. autoclass:: PyMieSim.ParameterFitResult
   :members:
.. autofunction:: PyMieSim.fit_parameters
