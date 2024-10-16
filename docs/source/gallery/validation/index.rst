:orphan:

Validation Module
=================

The Validation module in PyMieSim is designed to ensure the accuracy and reliability of the simulation results.
This module contains a set of examples and tests that compare PyMieSim’s output with well-established theoretical results, such as those from *Bohren and Huffman* or *PyMieScatt*.
By providing these validation cases, the module allows users to verify the correctness of their simulations and offers confidence in the consistency of the results produced by the library.

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



.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    </div>



Validation with Bohren & Huffmann
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Cylinder Scatterer Bohren-Huffman figure 8.10">

.. only:: html

  .. image:: /gallery/validation/bohren_huffman/images/thumb/sphx_glr_figure_810_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_validation_bohren_huffman_figure_810.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder Scatterer Bohren-Huffman figure 8.10</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Cylinder Scatterer Bohren-Huffman figure 8.7">

.. only:: html

  .. image:: /gallery/validation/bohren_huffman/images/thumb/sphx_glr_figure_87_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_validation_bohren_huffman_figure_87.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder Scatterer Bohren-Huffman figure 8.7</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Cylinder Scatterer Bohren-Huffman figure 8.8">

.. only:: html

  .. image:: /gallery/validation/bohren_huffman/images/thumb/sphx_glr_figure_88_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_validation_bohren_huffman_figure_88.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Cylinder Scatterer Bohren-Huffman figure 8.8</div>
    </div>


.. raw:: html

    </div>



Internal validation
~~~~~~~~~~~~~~~~~~~


.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Goniometric Coupling vs S1 S2 Comparison">

.. only:: html

  .. image:: /gallery/validation/internal/images/thumb/sphx_glr_phase_function_detector_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_validation_internal_phase_function_detector.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Goniometric Coupling vs S1 S2 Comparison</div>
    </div>


.. raw:: html

    </div>



Validation with PyMieScatt
~~~~~~~~~~~~~~~~~~~~~~~~~~


.. raw:: html

    <div class="sphx-glr-thumbnails">


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere Particles: 0">

.. only:: html

  .. image:: /gallery/validation/pymiescatt/images/thumb/sphx_glr_sphere_0_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_validation_pymiescatt_sphere_0.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere Particles: 0</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Sphere Particles: 1">

.. only:: html

  .. image:: /gallery/validation/pymiescatt/images/thumb/sphx_glr_sphere_1_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_validation_pymiescatt_sphere_1.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Sphere Particles: 1</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Core-Shell Particles: 0">

.. only:: html

  .. image:: /gallery/validation/pymiescatt/images/thumb/sphx_glr_coreshell_0_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_validation_pymiescatt_coreshell_0.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Core-Shell Particles: 0</div>
    </div>


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Core-Shell Particles: 1">

.. only:: html

  .. image:: /gallery/validation/pymiescatt/images/thumb/sphx_glr_coreshell_1_thumb.png
    :alt:

  :ref:`sphx_glr_gallery_validation_pymiescatt_coreshell_1.py`

.. raw:: html

      <div class="sphx-glr-thumbnail-title">Core-Shell Particles: 1</div>
    </div>


.. raw:: html

    </div>


.. toctree::
   :hidden:
   :includehidden:


   /gallery/validation/bohren_huffman/index.rst
   /gallery/validation/internal/index.rst
   /gallery/validation/pymiescatt/index.rst



.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
