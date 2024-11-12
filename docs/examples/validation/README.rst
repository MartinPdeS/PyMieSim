.. _validation_index:

Validation Module
=================

The Validation module in PyMieSim is designed to ensure the accuracy and reliability of the simulation results.
This module contains a set of examples and tests that compare PyMieSim’s output with well-established theoretical results, such as those from *Bohren and Huffman* or *PyMieScatt*.
By providing these validation cases, the module allows users to verify the correctness of their simulations and offers confidence in the consistency of the results produced by the library.

.. _key_components_validation:

Key Components
--------------

1. **Bohren and Huffman Validation**:
   PyMieSim includes several examples that compare its results against the classical results presented in *Bohren and Huffman's* "Absorption and Scattering of Light by Small Particles".
   These examples are used to validate calculations of scattering cross-sections, extinction efficiencies, and more.

   - `figure_87`: Validates scattering efficiency against theoretical values.
   - `figure_88`: Compares extinction cross-section results.

2. **PyMieScatt Validation**:
   PyMieSim is validated against the *PyMieScatt* library, a Python package that implements Mie theory for single-sphere scattering problems.
   The module includes several test cases comparing the results of PyMieSim with those of PyMieScatt, ensuring consistent calculations.

   - `coreshell_0`: Compares scattering properties of core-shell particles.
   - `sphere_0`: Validates scattering results for a single homogeneous sphere.

3. **Internal Validation**:
   In addition to external references, PyMieSim includes internal validation cases that test the correctness of specific components of the library, such as the implementation of phase functions, detector configurations, and source models.

   - `phase_function_detector`: Tests the accuracy of phase function calculations for different detector types.
   - `energy_flow_vs_coupling`: Verifies the energy flow in coupling experiments.
