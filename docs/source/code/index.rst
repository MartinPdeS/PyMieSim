Source Code Structures
=======================

This section provides an overview of the classes and components available in the PyMieSim package. Each component listed below inherits relevant attributes and methods, allowing users to build flexible and customizable scattering simulations.
This overview serves as a guide for understanding how the package is organized and how its parts can be used to simulate different optical scattering scenarios.

The PyMieSim package is divided into two main modules:

1. **Single Module**: Focuses on simulating scattering by individual particles such as spheres, cylinders, and core-shell structures. It provides detailed control over single-particle interactions and their scattering behavior.

2. **Experiment Module**: Supports more complex experimental setups with multiple components, such as sources, detectors, and scatterers, enabling users to simulate sophisticated scattering experiments.

Each module contains a variety of classes and methods that can be customized to suit the user's needs, whether you're analyzing basic scattering properties or configuring advanced detector systems.

.. toctree::
    :maxdepth: 1
    :caption: Code Structure

    single/index.rst
    experiment/index.rst

Classes and Components
----------------------

- **Source**: Defines the source of light for the simulation, with options such as PlaneWave and GaussianBeam. The source class is responsible for setting the wavelength, polarization, and other beam parameters.

- **Scatterer**: The object that interacts with the light source, producing scattered light. PyMieSim supports various scatterer types, including spheres, cylinders, and core-shell particles.

- **Detector**: Captures the scattered light and measures properties like intensity, polarization, and phase. The detector can be placed at different angles or configured to capture light from specific regions around the scatterer.

- **Experiment**: A class that brings together the source, scatterer, and detector, creating a complete experimental setup for optical scattering simulations.

Inheritance Structure
---------------------

Each class in PyMieSim inherits relevant attributes and methods from parent classes, allowing for a flexible and modular design.
This structure ensures that users can easily extend the functionality of individual components without rewriting the core logic of the package.

