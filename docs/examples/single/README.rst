.. _single_index:

Single Module
=============

The Single module in PyMieSim is designed for simulating light scattering by individual particles.
It provides tools for defining single scatterers, their sources, and detectors, offering a detailed analysis of optical properties at the level of a single interaction.
This module is ideal for researchers and engineers interested in exploring the scattering behavior of individual particles in isolation.

Key Components
--------------

1. **Sources**:
   The `single` module offers several types of sources to illuminate the scatterer.
   These sources can be configured to match specific experimental conditions, such as wavelength, polarization, and propagation direction.

   - `planewave`: Simulates a plane wave source, a uniform beam commonly used in scattering studies.
   - `gaussian`: Simulates a Gaussian beam, useful for focused light interactions with scatterers.

2. **Scatterers**:
   Defines the individual particles being studied.
   PyMieSim supports various geometries such as spheres, cylinders, and core-shell structures, allowing for the investigation of how shape, size, and material composition affect scattering properties.

   - `sphere`: A single spherical particle, modeled using Mie theory.
   - `cylinder`: A cylindrical scatterer, ideal for elongated particles like fibers.
   - `core_shell`: Defines a layered particle with different materials in the core and shell, enabling studies of multi-material systems.

3. **Detectors**:
   Detectors capture the scattered light from the single scatterer.
   Various detector types are available, each offering different capabilities for measuring the properties of scattered light.

   - `photodiode`: A simple detector for measuring intensity.
   - `coherent`: Captures the coherence properties of the scattered field, providing insights into phase and amplitude relationships.
   - `uncoherent`: Measures intensity without regard to phase information.

4. **Single Scatterer Setup**:
   The module allows for a flexible combination of source, scatterer, and detector, offering detailed simulations of how light interacts with individual particles.
   You can explore a wide range of properties, such as scattering cross-sections, phase functions, polarization effects, and more.

   Example configurations include:
   - A spherical scatterer illuminated by a Gaussian beam with intensity measurements at various angles.
   - Core-shell particles analyzed for extinction cross-sections as a function of shell thickness.
   - Cylindrical scatterers examined for scattering efficiencies across different wavelengths.

Example Usage
-------------

Below is an example of how to set up and run a simulation using the `single` module in PyMieSim:

.. code-block:: python

    from PyMieSim.single.scatterer import Sphere
    from PyMieSim.single.source import Gaussian
    from PyMieSim.single.detector import Photodiode
    from PyMieSim.units import nanometer, degree, watt, AU, RIU

    source = Gaussian(
        wavelength=450 * nanometer,
        polarization=0 * degree,
        optical_power=1 * watt,
        NA=0.3 * AU
    )

    scatterer = Sphere(
        diameter=6 * nanometer,  # 6 nm
        source=source,
        medium_property=1.0 * RIU,
        property=1.4 * RIU
    )

    detector = Photodiode(
        NA=0.1 * AU,
        phi_offset=0 * degree,
        gamma_offset=0 * degree,
        sampling=200 * AU,
        polarization_filter=None
    )

    coupling = detector.coupling(scatterer)

    print(coupling)


The `Single` module simplifies the process of analyzing light scattering at the individual particle level, providing an intuitive interface for exploring the interaction between light and single particles.
