:orphan:

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



.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    </div>

Core-Shell
~~~~~~~~~~


.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the B1 scattering parameter as a functio...">

.. only:: html

  .. image:: /gallery/experiment/coreshell/images/thumb/sphx_glr_coreshell_b1_vs_corediameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_coreshell_coreshell_b1_vs_corediameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">CoreShell: B1 vs Core Diameter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the B1 scattering parameter as a functio...">

.. only:: html

  .. image:: /gallery/experiment/coreshell/images/thumb/sphx_glr_coreshell_a1_vs_corediameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_coreshell_coreshell_a1_vs_corediameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">CoreShell: A1 vs Core Diameter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the backscattering efficiency (Qback) as...">

.. only:: html

  .. image:: /gallery/experiment/coreshell/images/thumb/sphx_glr_coreshell_Qback_vs_corediameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_coreshell_coreshell_Qback_vs_corediameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">CoreShell: Qback vs Core Diameter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the coupling efficiency as a function of...">

.. only:: html

  .. image:: /gallery/experiment/coreshell/images/thumb/sphx_glr_coreshell_coupling_vs_corediameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_coreshell_coreshell_coupling_vs_corediameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">CoreShell: Coupling vs Diameter</div>
    </div>


.. raw:: html

    </div>

Cylinder
~~~~~~~~


.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Cylinder: Qsca vs wavelength std">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_Qsca_vs_wavelength_std_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_Qsca_vs_wavelength_std.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: Qsca vs wavelength std</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the A1 scattering coefficient as a funct...">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_a11_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_a11_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: A1 Scattering Coefficient</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the B1 scattering coefficient as a funct...">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_b11_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_b11_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: B1 Scattering Coefficient</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a fu...">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_Qsca_vs_index_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_Qsca_vs_index.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: Qsca vs Index</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a fu...">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_Qsca_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_Qsca_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: Qsca vs Diameter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a fu...">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_Qsca_vs_wavelength_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_Qsca_vs_wavelength.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: Qsca vs Wavelength</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the scattering efficiency (Qsca) as a fu...">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_Qabs_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_Qabs_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: Qabs vs Diameter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the coupling efficiency as a function of...">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_coupling_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_coupling_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: Coupling vs Diameter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to use a goniometer setup to measure and visualize the coupling e...">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_coupling_vs_phioffset_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_coupling_vs_phioffset.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: Goniometer</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example demonstrates how to compute and visualize the coupling efficiency as a function of...">

.. only:: html

  .. image:: /gallery/experiment/cylinder/images/thumb/sphx_glr_cylinder_coupling_vs_wavelength_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_cylinder_cylinder_coupling_vs_wavelength.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder: Coupling vs Wavelength</div>
    </div>


.. raw:: html

    </div>

Sphere
~~~~~~


.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: B1 scattering coefficient">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_b1_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_b1_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: B1 scattering coefficient</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Qsca vs index">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_Qsca_vs_index_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_Qsca_vs_index.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Qsca vs index</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Qsca vs wavelength">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_Qsca_vs_wavelength_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_Qsca_vs_wavelength.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Qsca vs wavelength</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: A1 scattering coefficient">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_a1_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_a1_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: A1 scattering coefficient</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Qabs vs diameter">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_Qabs_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_Qabs_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Qabs vs diameter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Qsca vs wavelength STD">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_Qsca_vs_wavelength_std_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_Qsca_vs_wavelength_std.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Qsca vs wavelength STD</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Qsca vs diameter">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_Qsca_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_Qsca_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Qsca vs diameter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Coupling vs sampling">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_coupling_vs_sampling_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_coupling_vs_sampling.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Coupling vs sampling</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Coherent mode field rotation">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_coherent_coupling_vs_rotation_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_coherent_coupling_vs_rotation.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Coherent mode field rotation</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Goniometer">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_coupling_vs_phioffset_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_coupling_vs_phioffset.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Goniometer</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: coherent coupling vs sampling">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_coherent_coupling_vs_sampling_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_coherent_coupling_vs_sampling.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: coherent coupling vs sampling</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Coupling vs polarization filter">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_coupling_vs_filter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_coupling_vs_filter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Coupling vs polarization filter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Coupling vs diameter">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_coupling_vs_diameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_coupling_vs_diameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Coupling vs diameter</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Coupling vs wavelength">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_coupling_vs_wavelength_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_coupling_vs_wavelength.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Coupling vs wavelength</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Coherent Goniometer">

.. only:: html

  .. image:: /gallery/experiment/sphere/images/thumb/sphx_glr_sphere_coherent_coupling_vs_phioffset_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_experiment_sphere_sphere_coherent_coupling_vs_phioffset.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Coherent Goniometer</div>
    </div>


.. raw:: html

    </div>


.. toctree::
   :hidden:
   :includehidden:


   /gallery/experiment/coreshell/index.rst
   /gallery/experiment/cylinder/index.rst
   /gallery/experiment/sphere/index.rst



.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
