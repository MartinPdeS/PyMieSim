API Reference
=============

Start with :doc:`public_api` to choose an interface and understand its return
values. ``Simulation`` configures one particle and illumination; ``Experiment``
evaluates parameter combinations. An experiment is a parameter study, not a
many-particle electromagnetic solver.

Classes and Components
----------------------

- **Source**: Defines the source of light for the simulation, with options such as ``PlaneWave`` and ``Gaussian``. The source class is responsible for setting the wavelength, polarization, and other beam parameters.

- **Scatterer**: The object that interacts with the light source, producing scattered light. PyMieSim supports various scatterer types, including spheres, cylinders, and core-shell particles.

- **Detector**: Captures the scattered light and measures properties like intensity, polarization, and phase. The detector can be placed at different angles or configured to capture light from specific regions around the scatterer.

- **Experiment**: Combines source, scatterer, and optional detector sets to evaluate a parameter grid.

Explore the public interfaces
-----------------------------

.. grid:: 1 2 3 3
   :gutter: 2

   .. grid-item-card:: Stable Python API
      :link: public_api
      :link-type: doc

      Recommended entry point for the supported public interface.

   .. grid-item-card:: Single module
      :link: single
      :link-type: doc

      Configure one particle, source, and detector.

   .. grid-item-card:: Experiment module
      :link: experiment
      :link-type: doc

      Sweep combinations of sources, scatterers, and detectors.

   .. grid-item-card:: Materials
      :link: material
      :link-type: doc

      Define refractive-index models and material data.

   .. grid-item-card:: Polarization
      :link: polarization
      :link-type: doc

      Describe polarization states and polarization sets.

.. toctree::
   :hidden:
   :maxdepth: 1

   public_api.rst
   single.rst
   experiment.rst
   material.rst
   polarization.rst
