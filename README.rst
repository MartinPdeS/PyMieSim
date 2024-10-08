|logo|

|python|
|zenodo|
|colab|
|coverage|
|docs|
|PyPi|
|PyPi_download|


PyMieSim
========

**PyMieSim** is a Python library designed to provide a robust and flexible framework for performing Mie scattering simulations. The software is easy to install and operate, making it accessible to both new users and experienced researchers. PyMieSim enables users to explore the scattering properties of particles under various configurations, and is tailored for investigating single scattering events, as well as conducting large-scale parametric experiments.

At its core, PyMieSim includes three solvers optimized for different types of scatterers:

- **Spherical particles**
- **Infinite cylindrical particles**
- **Core-shell spherical particles**

The software also allows the user to customize the light source and detector attributes, depending on the specific simulation needs. The package is modular and provides an intuitive interface for users to model complex scattering scenarios with minimal effort.

|code_structure|

### Main Submodules

PyMieSim is organized into two primary submodules:

1. **single**: Focused on analyzing individual scattering events, such as:
   - Far-field distributions
   - Scattering phase functions
   - Stokes parameters

2. **experiment**: Designed for exploring how scattering parameters, such as `Qsca`, `Qext`, `g`, and `coupling (power)`, behave over large datasets, incorporating variations in sources, scatterers, and detectors.

Both submodules work seamlessly together, making PyMieSim adaptable for a wide range of scattering simulations.


----

## Getting Started

To use PyMieSim in Python, simply install the package and begin incorporating it into your scripts.

### Installation

PyMieSim supports Windows, Linux, macOS (including Apple M1/M2 chips), and ARM architectures. To install the package, use pip:

.. code-block:: bash

    pip install PyMieSim

For more details, visit the `documentation <https://pymiesim.readthedocs.io/en/latest/>`_ for a comprehensive guide on how to use the package.

----

## Example Code

Here is an example of how to use PyMieSim for a simple Mie scattering simulation. This example demonstrates how to configure a light source, scatterer, and detector, and retrieve the scattering data:

.. code-block:: python

    import numpy as np

    from PyMieSim.experiment.scatterer import Sphere
    from PyMieSim.experiment.source import Gaussian
    from PyMieSim.experiment import Setup
    from PyMieSim.units import nanometer, degree, watt, AU, RIU

    source = Gaussian(
        wavelength=np.linspace(400, 1000, 500) * nanometer,
        polarization=0 * degree,
        optical_power=1e-3 * watt,
        NA=0.2 * AU
    )

    scatterer = Sphere(
        diameter=[200] * nanometer,
        property=[4] * RIU,
        medium_property=1 * RIU,
        source=source
    )

    experiment = Setup(scatterer=scatterer, source=source)

    dataframe = experiment.get('Qsca')

    dataframe.plot_data(x="wavelength")


It produces the following figure which is equivalent to the one found on `wikipedia <https://en.wikipedia.org/wiki/Mie_scattering#/media/File:N4wiki.svg>`_.

|wikipedia_example|


This is just one example of PyMieSim in action. You can find more examples in the
`examples section <https://pymiesim.readthedocs.io/en/master/gallery/index.html>`_ of the documentation.

----

## Examples

Here are a few more examples showcasing the capabilities of PyMieSim:

### Example 1: Plasmonic Resonances for CoreShell Particles

|example_plasmon|

### Example 2: Scattering Efficiency vs Diameter for Spherical Particles

|example_qsca|

----

## Manual Building

If you prefer or need to build the project manually (e.g., for Apple silicon devices), ensure you have a C++ compiler (such as gcc) and Fortran installed, as well as Python 3.7+.

### Build Instructions

Linux/MacOS
~~~~~~~~~~~

.. code-block:: bash

    git clone https://github.com/MartinPdeS/PyMieSim.git
    cd PyMieSim
    git submodule init && git submodule update
    mkdir build
    cd build
    cmake ../ -G"Unix Makefiles"
    sudo make install
    cd ..
    python -m pip install .

For Windows, use `MinGW Makefiles` instead of `Unix Makefiles` when invoking CMake.

----

## Testing

You can test the local version of PyMieSim by running the following commands:

.. code-block:: bash

    git clone https://github.com/MartinPdeS/PyMieSim.git
    cd PyMieSim
    pip install PyMieSim[testing]
    pytest

This will run the suite of unit tests and provide coverage details.

----

## Google Colab

In 2024, running code on your local machine is optional! You can leverage the power of Google Colab to run PyMieSim remotely. Use the provided
`Colab notebook <https://colab.research.google.com/github/MartinPdeS/PyMieSim/blob/master/notebook.ipynb>`_ for an interactive experience.

|colab|

----

## Citing PyMieSim

If PyMieSim contributes to your research, we kindly ask that you cite the following paper:

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

You can access the full article `here <https://opg.optica.org/optcon/fulltext.cfm?uri=optcon-2-3-520&id=526697>`_

----

## Experimental Graphical User Interface (GUI)

Since version 1.7.0, PyMieSim offers an experimental GUI for users who prefer a graphical approach to simulations. While still under development, the GUI can be installed and accessed as follows:

.. code-block:: bash

    pip install PyMieSim
    python -m PyMieSim

The GUI is not yet as robust as the core Python API, but it provides a simplified interface for generating simulations.

|example_gui|

----

## Contact Information

PyMieSim is actively developed and maintained by Martin Poinsinet de Sivry-Houle. If you're interested in contributing or have questions, feel free to reach out.

Email: `martin.poinsinet.de.sivry@gmail.ca <mailto:martin.poinsinet.de.sivry@gmail.ca?subject=PyMieSim>`_

----


.. |python| image:: https://img.shields.io/pypi/pyversions/pymiesim.svg
    :target: https://www.python.org/

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5593704.svg
    :target: https://doi.org/10.5281/zenodo.4556074

.. |colab| image:: https://colab.research.google.com/assets/colab-badge.svg
    :target: https://colab.research.google.com/github/MartinPdeS/PyMieSim/blob/master/notebook.ipynb

.. |docs| image:: https://github.com/martinpdes/pymiesim/actions/workflows/deploy_documentation.yml/badge.svg
    :target: https://martinpdes.github.io/PyMieSim/
    :alt: Documentation Status

.. |PyPi| image:: https://badge.fury.io/py/PyMieSim.svg
    :target: https://badge.fury.io/py/PyMieSim

.. |logo| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/logo.png

.. |example_plasmon| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/plasmonic_resonances.png

.. |example_qsca| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/Qsca_diameter.png

.. |PyPi_download| image:: https://img.shields.io/pypi/dm/PyMieSim.svg
    :target: https://pypistats.org/packages/pymiesim

.. |code_structure| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/code_structure.png
    :width: 800
    :alt: Structure of the library

.. |example_gui| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/example_gui.png
    :width: 800
    :alt: Structure of the library

.. |coverage| image:: https://raw.githubusercontent.com/MartinPdeS/PyMieSim/python-coverage-comment-action-data/badge.svg
    :alt: Unittest coverage
    :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyMieSim/blob/python-coverage-comment-action-data/htmlcov/index.html

.. |wikipedia_example| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/wikipedia_example.png

