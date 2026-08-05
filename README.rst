|logo|

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - Badge
     - Status
   * - Python versions
     - |python|
   * - Documentation
     - |docs|
   * - Scientific article
     - |article|
   * - Continuous integration
     - |ci/cd|
   * - Test coverage
     - |coverage|
   * - Google Colab
     - |colab|
   * - PyPI package
     - |PyPI|
   * - PyPI downloads
     - |PyPI_download|
   * - Anaconda package
     - |anaconda|
   * - Anaconda downloads
     - |anaconda_download|
   * - Latest Anaconda release
     - |anaconda_date|

PyMieSim
========

**PyMieSim** is an open-source Python package for fast and flexible Mie scattering simulations.
It supports spherical, cylindrical and core--shell particles and provides helper classes for custom sources and detectors.
The project targets both quick single-scatterer studies and large parametric experiments.

Try the live web GUI: `PyMieSim Parameter Sweep Lab <https://pymiesim.onrender.com/>`_.

Features
--------
- Solvers for spheres, cylinders and core--shell geometries.
- Built-in models for plane wave and Gaussian sources.
- Multiple detector types including photodiodes and coherent modes.
- Simple data analysis with pandas DataFrame outputs.

Installation
------------
PyMieSim is available on PyPI and Anaconda.  Install it with:

.. code-block:: bash

   pip install PyMieSim
   conda install PyMieSim  --channels MartinPdeS

See the `online documentation <https://martinpdes.github.io/PyMieSim/>`_ for detailed usage and additional examples.

Quick example
-------------
Below is a short example computing the scattering efficiency of a sphere.

.. code-block:: python

    from PyMieSim import Gaussian, PolarizationState, Simulation, Sphere, ureg

    source = Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1e-3 * ureg.watt,
        numerical_aperture=0.2,
    )

    scatterer = Sphere(
        diameter=200 * ureg.nanometer,
        material=4 + 1j,
        medium=1.0,
    )

    simulation = Simulation(scatterer=scatterer, source=source)
    qsca = simulation.run("Qsca")
    print(qsca)

For wavelength, size, or material sweeps, use the experiment API described in
the `experiment examples <https://martinpdes.github.io/PyMieSim/examples.html>`_.


.. image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/resonances.png
    :width: 1000
    :align: center
    :alt: Scattering efficiency of a 200 nm sphere with refractive index 4.0.


Code structure
---------------
Here is the architecture for a standard workflow using PyMieSim:

.. image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/code_structure.png
   :width: 1000
   :align: center
   :alt: Code structure of a standard workflow using PyMieSim.

Building from source
--------------------
For development or manual compilation, clone the repository and run:

.. code-block:: bash

   git submodule update --init
   mkdir build && cd build
   cmake ../ -G"Unix Makefiles"
   sudo make install
   cd ..
   python -m pip install .

Testing
-------
Run the unit tests with:

.. code-block:: bash

   pip install PyMieSim[testing]
   pytest

Citing PyMieSim
---------------
If you use PyMieSim in academic work, please cite:

.. code-block:: none

   @article{PoinsinetdeSivry-Houle:23,
       author = {Martin Poinsinet de Sivry-Houle and Nicolas Godbout and Caroline Boudoux},
       journal = {Opt. Continuum},
       title = {PyMieSim: an open-source library for fast and flexible far-field Mie scattering simulations},
       volume = {2},
       number = {3},
       pages = {520--534},
       year = {2023},
       doi = {10.1364/OPTCON.473102},
   }

Contact
-------
For questions or contributions, contact `martin.poinsinet.de.sivry@gmail.com <mailto:martin.poinsinet.de.sivry@gmail.com>`_.

.. |logo| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/logo.png
    :alt: PyOptik logo
.. |python| image:: https://img.shields.io/pypi/pyversions/pymiesim.svg
    :alt: Python
    :target: https://www.python.org/
.. |article| image:: https://img.shields.io/badge/Optics%20Continuum-PyMieSim-green.svg
    :alt: Scientific article
    :target: https://opg.optica.org/optcon/viewmedia.cfm?uri=optcon-2-3-520&html=true
.. |colab| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Google Colab
    :target: https://colab.research.google.com/github/MartinPdeS/PyMieSim/blob/master/notebook.ipynb
.. |docs| image:: https://github.com/martinpdes/pymiesim/actions/workflows/deploy_documentation.yml/badge.svg
    :target: https://martinpdes.github.io/PyMieSim/
    :alt: Documentation Status
.. |PyPI| image:: https://badge.fury.io/py/PyMieSim.svg
    :alt: PyPI version
    :target: https://badge.fury.io/py/PyMieSim
.. |PyPI_download| image:: https://api.pepy.tech/badge/PyMieSim/month
    :alt: PyPI downloads
    :target: https://pepy.tech/projects/pymiesim
.. |coverage| image:: https://raw.githubusercontent.com/MartinPdeS/PyMieSim/python-coverage-comment-action-data/badge.svg
    :alt: Unittest coverage
    :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html
.. |ci/cd| image:: https://github.com/martinpdes/pymiesim/actions/workflows/deploy_coverage.yml/badge.svg
    :alt: Unittest Status
.. |wikipedia_example| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/wikipedia_example.png
    :width: 800
    :alt: Example wikipedia
.. |example_plasmon| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/plasmonic_resonances.png
    :width: 800
    :alt: Plasmonic resonances
.. |example_qsca| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/Qsca_diameter.png
    :width: 800
    :alt: Qsca vs diameter
.. |anaconda| image:: https://anaconda.org/martinpdes/pymiesim/badges/version.svg
    :alt: Anaconda version
    :target: https://anaconda.org/martinpdes/pymiesim
.. |anaconda_download| image:: https://anaconda.org/martinpdes/pymiesim/badges/downloads.svg
    :alt: Anaconda downloads
    :target: https://anaconda.org/martinpdes/pymiesim
.. |anaconda_date| image:: https://anaconda.org/martinpdes/pymiesim/badges/latest_release_relative_date.svg
    :alt: Latest release date
    :target: https://anaconda.org/martinpdes/pymiesim
