
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "gallery/experiment/cylinder/cylinder_coupling_vs_wavelength.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_gallery_experiment_cylinder_cylinder_coupling_vs_wavelength.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_gallery_experiment_cylinder_cylinder_coupling_vs_wavelength.py:


Cylinder: Coupling vs Wavelength
================================

This example demonstrates how to compute and visualize the coupling efficiency as a function of wavelength for cylindrical scatterers using PyMieSim.

.. GENERATED FROM PYTHON SOURCE LINES 9-10

Importing the package dependencies: numpy, PyMieSim

.. GENERATED FROM PYTHON SOURCE LINES 10-18

.. code-block:: python3

    import numpy as np
    from PyMieSim.experiment.detector import CoherentMode
    from PyMieSim.experiment.scatterer import Cylinder
    from PyMieSim.experiment.source import Gaussian
    from PyMieSim.experiment import Setup
    from PyOptik import Material
    from PyMieSim.units import nanometer, degree, watt, AU, RIU








.. GENERATED FROM PYTHON SOURCE LINES 19-20

Defining the source

.. GENERATED FROM PYTHON SOURCE LINES 20-27

.. code-block:: python3

    source = Gaussian(
        wavelength=np.linspace(950, 1050, 300) * nanometer,  # Wavelengths ranging from 950 nm to 1050 nm
        polarization=0 * degree,  # Linear polarization angle in radians
        optical_power=1e-3 * watt,  # 1 milliwatt
        NA=0.2 * AU  # Numerical Aperture
    )








.. GENERATED FROM PYTHON SOURCE LINES 28-30

Defining the scatterer distribution
Here we look at cylinders with a set diameter, refractive index, and medium.

.. GENERATED FROM PYTHON SOURCE LINES 30-37

.. code-block:: python3

    scatterer = Cylinder(
        diameter=np.linspace(100, 8000, 5) * nanometer,  # Diameters ranging from 100 nm to 8000 nm
        property=Material.BK7,  # Material of the cylinder
        medium_property=1 * RIU,  # Refractive index of the surrounding medium
        source=source
    )








.. GENERATED FROM PYTHON SOURCE LINES 38-39

Defining the detector

.. GENERATED FROM PYTHON SOURCE LINES 39-49

.. code-block:: python3

    detector = CoherentMode(
        mode_number="LP11",  # Specifying the LP11 mode
        NA=[0.05, 0.01] * AU,  # Array of Numerical Apertures for the detector
        phi_offset=-180 * degree,  # Phi offset in degrees
        gamma_offset=0 * degree,  # Gamma offset in degrees
        polarization_filter=None,  # No polarization filter
        sampling=300 * AU,  # Number of sampling points
        rotation=0 * degree,  # Rotation of the mode field
    )








.. GENERATED FROM PYTHON SOURCE LINES 50-51

Setting up the experiment

.. GENERATED FROM PYTHON SOURCE LINES 51-53

.. code-block:: python3

    experiment = Setup(scatterer=scatterer, source=source, detector=detector)








.. GENERATED FROM PYTHON SOURCE LINES 54-55

Measuring the coupling efficiency

.. GENERATED FROM PYTHON SOURCE LINES 55-57

.. code-block:: python3

    dataframe = experiment.get('coupling', scale_unit=True)





.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    dict_keys(['source:wavelength', 'source:polarization', 'source:NA', 'source:optical_power', 'scatterer:medium_property', 'scatterer:diameter', 'scatterer:property', 'detector:mode_number', 'detector:NA', 'detector:phi_offset', 'detector:gamma_offset', 'detector:sampling', 'detector:rotation', 'detector:polarization_filter'])




.. GENERATED FROM PYTHON SOURCE LINES 58-60

Plotting the results
Visualizing how the coupling efficiency varies with the wavelength.

.. GENERATED FROM PYTHON SOURCE LINES 60-60

.. code-block:: python3

    dataframe.plot_data(x="source:wavelength", std='scatterer:diameter')


.. image-sg:: /gallery/experiment/cylinder/images/sphx_glr_cylinder_coupling_vs_wavelength_001.png
   :alt: cylinder coupling vs wavelength
   :srcset: /gallery/experiment/cylinder/images/sphx_glr_cylinder_coupling_vs_wavelength_001.png
   :class: sphx-glr-single-img






.. rst-class:: sphx-glr-timing

   **Total running time of the script:** (0 minutes 1.278 seconds)


.. _sphx_glr_download_gallery_experiment_cylinder_cylinder_coupling_vs_wavelength.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example




    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: cylinder_coupling_vs_wavelength.py <cylinder_coupling_vs_wavelength.py>`

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: cylinder_coupling_vs_wavelength.ipynb <cylinder_coupling_vs_wavelength.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
