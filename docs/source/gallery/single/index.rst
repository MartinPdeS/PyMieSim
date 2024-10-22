:orphan:

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



.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    </div>

Detectors
~~~~~~~~~


.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the initialization and visualization of an Integrating Sphere detecto...">

.. only:: html

  .. image:: /gallery/single/detector/images/thumb/sphx_glr_integrating_sphere_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_detector_integrating_sphere.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Integrating sphere</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the initialization and visualization of an LP11 Mode detector using P...">

.. only:: html

  .. image:: /gallery/single/detector/images/thumb/sphx_glr_LP11_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_detector_LP11.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">LP11 Mode Detector</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the initialization and visualization of an LP01 Mode detector using P...">

.. only:: html

  .. image:: /gallery/single/detector/images/thumb/sphx_glr_LP01_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_detector_LP01.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">LP01 Mode Detector</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the initialization and visualization of an LP02 Mode detector using P...">

.. only:: html

  .. image:: /gallery/single/detector/images/thumb/sphx_glr_LP02_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_detector_LP02.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">LP02 Mode Detector</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the initialization and visualization of HG01 Mode detector using PyMi...">

.. only:: html

  .. image:: /gallery/single/detector/images/thumb/sphx_glr_HG01_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_detector_HG01.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Hermite-Gauss 01 Mode Detector</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the initialization and visualization of HG01 Mode detector using PyMi...">

.. only:: html

  .. image:: /gallery/single/detector/images/thumb/sphx_glr_LG11_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_detector_LG11.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Laguerre-Gauss 2-3 Mode Detector</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the initialization and visualization of a Photodiode detector using P...">

.. only:: html

  .. image:: /gallery/single/detector/images/thumb/sphx_glr_photodiode_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_detector_photodiode.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Photodiode Detector</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the initialization and visualization of HG31 Mode detector using PyMi...">

.. only:: html

  .. image:: /gallery/single/detector/images/thumb/sphx_glr_HG11_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_detector_HG11.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Hermite-Gauss 31 Mode Detector</div>
    </div>


.. raw:: html

    </div>

Scatterers
~~~~~~~~~~



.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to visualize the properties of a light source using PyMieSim.">

.. only:: html

  .. image:: /gallery/single/scatterer/images/thumb/sphx_glr_source_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_scatterer_source.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Source Plottings</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the computation of scattering properties using PyMieSim.">

.. only:: html

  .. image:: /gallery/single/scatterer/images/thumb/sphx_glr_properties_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_scatterer_properties.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Print properties</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the computation and visualization of the Stokes parameters using PyMi...">

.. only:: html

  .. image:: /gallery/single/scatterer/images/thumb/sphx_glr_stokes_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_scatterer_stokes.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Stokes Parameters Computation</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the S1 and S2 scattering functions using...">

.. only:: html

  .. image:: /gallery/single/scatterer/images/thumb/sphx_glr_s1s2_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_scatterer_s1s2.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">S1 S2 Function Computation</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the process of computing and visualizing the far-fields of a scattere...">

.. only:: html

  .. image:: /gallery/single/scatterer/images/thumb/sphx_glr_farfield_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_scatterer_farfield.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Far-Fields Computation and Visualization</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates the computation and visualization of the Scattering Phase Function (S...">

.. only:: html

  .. image:: /gallery/single/scatterer/images/thumb/sphx_glr_spf_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_scatterer_spf.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">SPF Computation</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the footprint of a scatterer using PyMie...">

.. only:: html

  .. image:: /gallery/single/scatterer/images/thumb/sphx_glr_footprint_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_single_scatterer_footprint.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Scatterer Footprint Calculation and Visualization</div>
    </div>


.. raw:: html

    </div>


.. toctree::
   :hidden:
   :includehidden:


   /gallery/single/detector/index.rst
   /gallery/single/scatterer/index.rst



.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
