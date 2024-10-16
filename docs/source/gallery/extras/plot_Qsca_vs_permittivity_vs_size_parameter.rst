
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "gallery/extras/plot_Qsca_vs_permittivity_vs_size_parameter.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_gallery_extras_plot_Qsca_vs_permittivity_vs_size_parameter.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_gallery_extras_plot_Qsca_vs_permittivity_vs_size_parameter.py:


Scattering efficiency of a sphere
=================================

PyMieSim makes it easy to create a source and a scatterer. With these objects
defined, it is possible to use PyMieSim to find the scattering efficiency of the
scatterer. This feature can be used to plot a graph of the scattering efficiency
of a sphere as a function of the permittivity and the size parameter.

.. GENERATED FROM PYTHON SOURCE LINES 12-13

Importing the package: PyMieSim

.. GENERATED FROM PYTHON SOURCE LINES 13-63

.. code-block:: python3

    import numpy

    from PyMieSim.experiment.scatterer import Sphere
    from PyMieSim.experiment.source import Gaussian
    from PyMieSim.experiment import Setup
    from PyMieSim.units import degree, watt, AU, nanometer, RIU

    import matplotlib.pyplot as plt


    permitivity = numpy.linspace(-10, 50, 400)

    index = numpy.sqrt(permitivity.astype(complex)) * RIU

    diameter = numpy.linspace(1, 200, 400) * nanometer

    source = Gaussian(
        wavelength=400 * nanometer,
        polarization=90 * degree,
        optical_power=1e-3 * watt,
        NA=0.2 * AU
    )


    scatterer = Sphere(
        diameter=diameter,
        property=index,
        medium_property=1 * RIU,
        source=source
    )

    experiment = Setup(
        scatterer=scatterer,
        source=source
    )

    data = experiment.get('Qsca', add_units=False).squeeze().values.reshape([permitivity.size, diameter.size])

    figure, ax = plt.subplots(1, 1)
    ax.set(
        xlabel="Permittivity",
        ylabel=r'Diameter [$\mu$ m]',
        title="Scattering efficiency of a sphere"
    )

    image = ax.pcolormesh(permitivity, diameter, numpy.log(data))

    plt.colorbar(mappable=image)

    plt.show()



.. image-sg:: /gallery/extras/images/sphx_glr_plot_Qsca_vs_permittivity_vs_size_parameter_001.png
   :alt: Scattering efficiency of a sphere
   :srcset: /gallery/extras/images/sphx_glr_plot_Qsca_vs_permittivity_vs_size_parameter_001.png
   :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    dict_keys(['source:wavelength', 'source:polarization', 'source:NA', 'source:optical_power', 'scatterer:medium_property', 'scatterer:diameter', 'scatterer:property'])





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** (0 minutes 2.534 seconds)


.. _sphx_glr_download_gallery_extras_plot_Qsca_vs_permittivity_vs_size_parameter.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example




    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: plot_Qsca_vs_permittivity_vs_size_parameter.py <plot_Qsca_vs_permittivity_vs_size_parameter.py>`

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: plot_Qsca_vs_permittivity_vs_size_parameter.ipynb <plot_Qsca_vs_permittivity_vs_size_parameter.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
