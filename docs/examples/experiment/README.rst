.. _experiment_index:

Experiment Module
=================

The Experiment module in PyMieSim provides the foundational tools to simulate complex optical scattering experiments. It enables users to define custom setups by specifying sources, detectors, and scatterers. This modular approach allows for the simulation of diverse configurations ranging from single scatterers to complex, multi-source systems.

Key Components
--------------

1. **Sources**:
   The experiment requires a light source, and PyMieSim supports various source types including plane waves, Gaussian beams. Each source can be configured with parameters such as wavelength, polarization, amplitude, providing flexibility for different simulation needs.

   - `planewave`: Simulates a uniform plane wave, commonly used in scattering experiments.
   - `gaussian`: Simulates a Gaussian beam source, useful for focused beam studies.

2. **Scatterers**:
   This component defines the particle or object that interacts with the light. The module supports different shapes and material configurations such as spheres, core-shell particles, and cylinders, allowing for customizable scattering scenarios.

   - `sphere`: Defines a spherical scatterer, the most common shape in Mie scattering theory.
   - `cylinder`: Allows for the study of cylindrical scatterers, which are useful for simulating fibers or other elongated structures.
   - `core_shell`: Supports core-shell structures, enabling the study of particles with layered materials.

3. **Detectors**:
   PyMieSim allows the placement of detectors around the scatterer to capture scattered light at various angles. These detectors can measure properties like intensity, polarization, and coherence.

   - `photodiode`: Simulates a basic detector to measure intensity.
   - `coherent_mode`: Captures information related to the coherence properties of the scattered field.
   - `integrating_sphere`: Simulates an integrating sphere detector, commonly used to measure total scattered light in all directions.

4. **Experiment Setup**:
   The module allows you to combine the source, scatterer, and detector into a complete experimental setup. You can simulate the interaction of light with the scatterer and capture the resulting scattered light using the specified detector.

Example Usage
-------------

Below is an example of how to set up a basic experiment using the `experiment` module:

.. code-block:: python

    from PyMieSim.experiment.scatterer import Sphere
    from PyMieSim.experiment.source import Gaussian
    from PyMieSim.experiment import Setup
    from PyMieSim.units import nanometer, degree, watt, AU, RIU

    source = Gaussian(
        wavelength=[500., 1000., 1500.] * nanometer,
        polarization=30. * degree,
        optical_power=1e-3 * watt,
        NA=0.2 * AU
    )

    scatterer = Sphere(
        diameter=800. * nanometer,
        property=np.linspace(1.3, 1.9, 150) * RIU,
        medium_property=1. * RIU,
        source=source
    )

    experiment = Setup(scatterer=scatterer, source=source)

    dataframe = experiment.get('Qsca')

    dataframe.plot_data(x="scatterer:property")


The `Experiment` module simplifies the process of creating and running scattering simulations, providing a streamlined workflow for exploring complex optical interactions.
