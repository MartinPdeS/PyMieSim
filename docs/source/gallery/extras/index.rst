:orphan:

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



.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="PyMieSim makes it easy to create a source and a scatterer. With these objects defined, it is po...">

.. only:: html

  .. image:: /gallery/extras/images/thumb/sphx_glr_plot_Qsca_vs_permittivity_vs_size_parameter_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_extras_plot_Qsca_vs_permittivity_vs_size_parameter.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Scattering efficiency of a sphere</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Samples Properties">

.. only:: html

  .. image:: /gallery/extras/images/thumb/sphx_glr_SphereProperties_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_extras_SphereProperties.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Samples Properties</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere: Coupling vs numerical aperture">

.. only:: html

  .. image:: /gallery/extras/images/thumb/sphx_glr_plot_coupling_vs_NA_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_extras_plot_coupling_vs_NA.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere: Coupling vs numerical aperture</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Example Script: Using the plot_system Function">

.. only:: html

  .. image:: /gallery/extras/images/thumb/sphx_glr_plot_system_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_extras_plot_system.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Plot system</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="PyMieSim makes it easy to create a source and a scatterer. With these objects defined, it is po...">

.. only:: html

  .. image:: /gallery/extras/images/thumb/sphx_glr_plot_coupling_heatmap_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_extras_plot_coupling_heatmap.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Coupling heatmap of a sphere</div>
    </div>


.. raw:: html

    </div>


.. toctree::
   :hidden:

   /gallery/extras/plot_Qsca_vs_permittivity_vs_size_parameter
   /gallery/extras/SphereProperties
   /gallery/extras/plot_coupling_vs_NA
   /gallery/extras/plot_system
   /gallery/extras/plot_coupling_heatmap



.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
