|logo|

.. list-table::
   :widths: 10 25 25 25
   :header-rows: 0

   * - Meta
     - |python|
     - |docs|
     - |zenodo|
   * - Testing
     - |ci/cd|
     - |coverage|
     - |colab|
   * - PyPI
     - |PyPI|
     - |PyPI_download|
     -
   * - Anaconda
     - |anaconda|
     - |anaconda_download|
     - |anaconda_date|

PyMieSim
========

**PyMieSim** is an open-source Python package for fast and flexible Mie scattering simulations.
It supports spherical, cylindrical and core--shell particles and provides helper classes for custom sources and detectors.
The project targets both quick single-scatterer studies and large parametric experiments.

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

   import numpy as np
   from TypedUnit import ureg

   from PyMieSim.experiment.scatterer import Sphere
   from PyMieSim.experiment.source import Gaussian
   from PyMieSim.experiment import Setup

   source = Gaussian(
       wavelength=np.linspace(400, 1000, 500) * ureg.nanometer,
       polarization=0 * ureg.degree,
       optical_power=1e-3 * ureg.watt,
       NA=0.2 * ureg.AU,
   )

   scatterer = Sphere(
       diameter=[200] * ureg.nanometer,
       property=[4] * ureg.RIU,
       medium_property=1 * ureg.RIU,
       source=source,
   )

   experiment = Setup(scatterer=scatterer, source=source)
   df = experiment.get("Qsca")
   df.plot(x="source:wavelength")

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
.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5593704.svg
    :alt: Scientific article
    :target: https://doi.org/10.5281/zenodo.4556074
.. |colab| image:: https://colab.research.google.com/assets/colab-badge.svg
    :alt: Google Colab
    :target: https://colab.research.google.com/github/MartinPdeS/PyMieSim/blob/master/notebook.ipynb
.. |docs| image:: https://github.com/martinpdes/pymiesim/actions/workflows/deploy_documentation.yml/badge.svg
    :target: https://martinpdes.github.io/PyMieSim/
    :alt: Documentation Status
.. |PyPI| image:: https://badge.fury.io/py/PyMieSim.svg
    :alt: PyPI version
    :target: https://badge.fury.io/py/PyMieSim
.. |PyPI_download| image:: https://img.shields.io/pypi/dm/PyMieSim?style=plastic&label=PyPI%20downloads&labelColor=hex&color=hex
    :alt: PyPI downloads
    :target: https://pypistats.org/packages/pymiesim
.. |coverage| image:: https://raw.githubusercontent.com/MartinPdeS/PyMieSim/python-coverage-comment-action-data/badge.svg
    :alt: Unittest coverage
    :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html
.. |ci/cd| image:: https://github.com/martinpdes/pymiesim/actions/workflows/deploy_coverage.yml/badge.svg
    :alt: Unittest Status
.. |example_gui| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/example_gui.png
    :width: 800
    :alt: Structure of the library
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
