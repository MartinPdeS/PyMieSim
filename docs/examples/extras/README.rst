Extras Examples
===============

PyMieSim can be used in many ways. Below is a gallery of examples showing different ways to use the library. These examples demonstrate more advanced usage of PyMieSim, showcasing custom setups and configurations that go beyond standard scattering simulations. The goal is to provide users with a variety of real-world applications that can inspire and guide further experimentation.

Key Examples
------------

1. **Sphere Properties**:
   This example demonstrates how to simulate and analyze the optical properties of spherical scatterers, such as their scattering and absorption cross-sections, as well as their interaction with different types of light sources.

   - *SphereProperties.py*: A script that computes and plots the scattering efficiency for a sphere of varying sizes.

2. **Coupling Heatmap**:
   Explore the coupling efficiency between a scatterer and a detector using different configurations. This example generates a heatmap to visualize how efficiency varies as a function of system parameters.

   - *plot_coupling_heatmap.py*: Generates a heatmap for coupling efficiency across a range of numerical apertures and scatterer properties.

3. **System Plotting**:
   PyMieSim allows you to create plots that illustrate the behavior of complex optical systems. This example demonstrates how to plot the entire system setup, including the source, scatterer, and detectors, providing a clear visualization of the experimental setup.

   - *plot_system.py*: A script that visualizes the experimental configuration for a scattering simulation.

4. **Scattering as a Function of Permittivity**:
   Investigate how scattering efficiency varies with the permittivity and size parameter of a particle. This is useful for understanding material-specific scattering behavior in optical simulations.

   - *plot_Qsca_vs_permittivity_vs_size_parameter.py*: Plots scattering cross-sections as a function of the scatterer's permittivity and size parameter.
