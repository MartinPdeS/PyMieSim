|Logo|

|python|
|zenodo|
|colab|
|unittest|
|docs|
|PyPi|
|PyPi_download|

PyMieSim
========

PyMieSim is a software designed for comprehensive Mie scattering analysis, featuring a user-friendly installation and operation process. The characterization of the scattering event within PyMieSim is determined by a set of specific properties, as illustrated in the subsequent figure.

Currently, PyMieSim integrates three distinct solvers tailored to three different types of scatterers: spherical particles, infinite cylindrical particles, and core-shell spherical particles. Additional parameters governing the scattering event are contingent upon the attributes of the light source and the detector (when applicable). The attributes pertinent to each of these components are delineated in the ensuing figure.


.. image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/code_structure.png
  :width: 800
  :alt: Structure of the library

The package also lets you construct an **Experiment** using **SphereSet**/**CoreShellSet**/**CylinderSet**, **SourceSet** and **DetectorSet**.
Those class define the type of scatterers, light sources and detectors you want to study.


----

Documentation
**************
All the latest available documentation is available `here <https://pymiesim.readthedocs.io/en/latest/>`_ or you can click the following badge:

|docs|

----

Google Colab
**************
It's 2023, you don't need to run all your code on you computer anymore. Google Colab is a platform which allows to write/use python scripts remotely.
You can open the PyMieSim.ipynb in the file to access it or click on the following "Open in Colab" badge:

|colab|

----


Installation
************

For common version of Windows, Linux and MacOS, (on x86_64 architecture), the package can readily be installed using pip;

.. code-block:: python

   >>> pip install PyMieSim

The ready to install wheel is not available for arm chip of the newer mac M1, M2 product. You can however install manually the package.


If, however, this fail you can build the package from scratch following the steps on the **Manual building** section.

**Note:** Wheel support now extended to `manylinux2014 <https://www.python.org/dev/peps/pep-0599/>`_.


----



Manual building
***************

To manually buld the project on your computer make sure that you do have gcc installed (c++ and fortran compiler), plus python version 3.7+.
For windows system I recommend install MingGW with g++ and fortran compiler.

This being done, the following commands should do the trick.

Linux / MacOs
~~~~~~~~~~~~~

.. code-block:: python

   >>> git clone https://github.com/MartinPdeS/PyMieSim.git
   >>> cd PyMieSim
   >>> git submodule init && git submodule update
   >>> mkdir build
   >>> cd build
   >>> cmake ../ -G"Unix MakeFiles" (macOS, Linux)
   >>> cmake ../ -G"MinGW MakeFiles" (Windows)
   >>> sudo make install
   >>> cd ..
   >>> python -m pip install .

----

Testing
*******

To test localy (with cloning the GitHub repository) you'll need to install the dependencies and run the coverage command as

.. code:: console

   pip install -r requirements/requirements.txt
   coverage run --source=<package> --module pytest --verbose <test-files-dirs> coverage report --show-missing

----


Coding examples
***************
Plenty of examples are available online, I invite you to check the `examples <https://pymiesim.readthedocs.io/en/latest/examples.html>`_
section of the documentation.


----

Scientific article
******************
The associated article is free of access on this link `article <https://opg.optica.org/optcon/fulltext.cfm?uri=optcon-2-3-520&id=526697>`_


Citing this work?
******************
I spent a full year to develop this tool for you to use so if it helped you in your research, I would greatly appreciate you citing the article associated to my work. Many thanks!

.. code-block:: none

   @article{PoinsinetdeSivry-Houle:23,
       author = {Martin Poinsinet de Sivry-Houle and Nicolas Godbout and Caroline Boudoux},
       journal = {Opt. Continuum},
       keywords = {Light scattering; Mie theory; Optical coherence tomography; Radiation pressure; Scattering theory; Surface plasmon resonance},
       number = {3},
       pages = {520--534},
       publisher = {Optica Publishing Group},
       title = {PyMieSim: an open-source library for fast and flexible far-field Mie scattering simulations},
       volume = {2},
       month = {Mar},
       year = {2023},
       url = {https://opg.optica.org/optcon/abstract.cfm?URI=optcon-2-3-520},
       doi = {10.1364/OPTCON.473102},
       abstract = {},
   }

----



Contact Information
************************
As of 2023, the project is still under development. If you want to collaborate, it would be a pleasure! I encourage you to contact me.

PyMieSim was written by `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_  .

Email:`martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=PyMieSim>`_ .


.. |travis| image:: https://img.shields.io/travis/com/MartinPdeS/PyMieSim/master?label=Travis%20CI
   :target: https://travis-ci.com/github/MartinPdeS/PyMieSim

.. |python| image:: https://img.shields.io/pypi/pyversions/pymiesim.svg
   :target: https://www.python.org/

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5593704.svg
   :target: https://doi.org/10.5281/zenodo.4556074

.. |colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/drive/1FUi_hRUXxCVvkHBY10YE1yR-nTATcDei?usp=sharing

.. |docs| image:: https://readthedocs.org/projects/pymiesim/badge/?version=latest
   :target: https://pymiesim.readthedocs.io/en/latest/

.. |PyPi| image:: https://badge.fury.io/py/PyMieSim.svg
    :target: https://badge.fury.io/py/PyMieSim

.. |Logo| image:: https://github.com/MartinPdeS/PyMieSim/raw/master/docs/images/logo.png

.. |PyPi_download| image:: https://img.shields.io/pypi/dm/PyMieSim.svg
   :target: https://pypi.org/project/PyMieSim/

.. |unittest| image:: https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/MartinPdeS/f0955be398d59efac69042c1b0fbece2/raw/b9e9d7f275011abe01054263702fc6f4a83273f3/PyMieSimcoverage_badge.json

