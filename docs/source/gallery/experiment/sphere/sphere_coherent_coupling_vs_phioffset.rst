
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "gallery/experiment/sphere/sphere_coherent_coupling_vs_phioffset.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_gallery_experiment_sphere_sphere_coherent_coupling_vs_phioffset.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_gallery_experiment_sphere_sphere_coherent_coupling_vs_phioffset.py:


Sphere: Coherent Goniometer
===========================

.. GENERATED FROM PYTHON SOURCE LINES 9-10

Importing the package dependencies: numpy, PyMieSim

.. GENERATED FROM PYTHON SOURCE LINES 10-19

.. code-block:: python3

    import numpy
    from PyMieSim.experiment.detector import CoherentMode
    from PyMieSim.experiment.scatterer import Sphere
    from PyMieSim.experiment.source import Gaussian
    from PyMieSim.units import degree
    from PyMieSim.experiment import Setup
    from PyOptik import Material
    from PyMieSim.units import nanometer, degree, watt, AU, RIU








.. GENERATED FROM PYTHON SOURCE LINES 20-21

Defining the source to be employed.

.. GENERATED FROM PYTHON SOURCE LINES 21-27

.. code-block:: python3

    source = Gaussian(
        wavelength=1200 * nanometer,
        polarization=90 * degree,
        optical_power=1e-3 * watt,
        NA=0.2 * AU
    )







.. GENERATED FROM PYTHON SOURCE LINES 28-29

Defining the ranging parameters for the scatterer distribution

.. GENERATED FROM PYTHON SOURCE LINES 29-36

.. code-block:: python3

    scatterer = Sphere(
        diameter=2000 * nanometer,
        property=Material.BK7,
        medium_property=1 * RIU,
        source=source
    )








.. GENERATED FROM PYTHON SOURCE LINES 37-38

Defining the detector to be employed.

.. GENERATED FROM PYTHON SOURCE LINES 38-48

.. code-block:: python3

    detector = CoherentMode(
        mode_number='LP11',
        NA=[0.5, 0.3, 0.1, 0.05] * AU,
        phi_offset=numpy.linspace(-180, 180, 300) * degree,
        gamma_offset=0 * degree,
        sampling=400 * AU,
        polarization_filter=10 * degree,
        rotation=0 * degree,  # Rotation of the mode field
    )








.. GENERATED FROM PYTHON SOURCE LINES 49-50

Defining the experiment setup

.. GENERATED FROM PYTHON SOURCE LINES 50-52

.. code-block:: python3

    experiment = Setup(scatterer=scatterer, source=source, detector=detector)








.. GENERATED FROM PYTHON SOURCE LINES 53-54

Measuring the properties

.. GENERATED FROM PYTHON SOURCE LINES 54-58

.. code-block:: python3

    dataframe = experiment.get('coupling')

    # # %%
    # # Plotting the results
    dataframe.plot_data(x="detector:phi_offset")


.. image-sg:: /gallery/experiment/sphere/images/sphx_glr_sphere_coherent_coupling_vs_phioffset_001.png
   :alt: sphere coherent coupling vs phioffset
   :srcset: /gallery/experiment/sphere/images/sphx_glr_sphere_coherent_coupling_vs_phioffset_001.png
   :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    dict_keys(['source:wavelength', 'source:polarization', 'source:NA', 'source:optical_power', 'scatterer:medium_property', 'scatterer:diameter', 'scatterer:property', 'detector:mode_number', 'detector:NA', 'detector:phi_offset', 'detector:gamma_offset', 'detector:sampling', 'detector:rotation', 'detector:polarization_filter'])





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** (0 minutes 0.649 seconds)


.. _sphx_glr_download_gallery_experiment_sphere_sphere_coherent_coupling_vs_phioffset.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example




    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: sphere_coherent_coupling_vs_phioffset.py <sphere_coherent_coupling_vs_phioffset.py>`

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: sphere_coherent_coupling_vs_phioffset.ipynb <sphere_coherent_coupling_vs_phioffset.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
