API Reference
=============

Start with :doc:`public_api` to choose an interface and understand its return
values. ``Simulation`` configures one particle and illumination; ``Experiment``
evaluates parameter combinations. An experiment is a parameter study, not a
many-particle electromagnetic solver.

.. toctree::
    :maxdepth: 1
    :caption: Code Structure

    public_api.rst
    single.rst
    experiment.rst
    material.rst
    polarization.rst

Classes and Components
----------------------

- **Source**: Defines the source of light for the simulation, with options such as ``PlaneWave`` and ``Gaussian``. The source class is responsible for setting the wavelength, polarization, and other beam parameters.

- **Scatterer**: The object that interacts with the light source, producing scattered light. PyMieSim supports various scatterer types, including spheres, cylinders, and core-shell particles.

- **Detector**: Captures the scattered light and measures properties like intensity, polarization, and phase. The detector can be placed at different angles or configured to capture light from specific regions around the scatterer.

- **Experiment**: Combines source, scatterer, and optional detector sets to evaluate a parameter grid.
