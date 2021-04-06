
.. image:: ./docs/images/Logo.png



|codecov|
|travis|
|python|
|zenodo|
|colab|
|docs|

PyMieSim
========



PyMieSim is a very easy to install/use tool for extensive Mie scattering analysis. It allows to study the light scattering
on different kind of object (scatterer), at the moment I only implemented spherical scatterers.
Using this package, one can easily set a **Source** a **Scatterer** and a **Detector** within a very wide range of parameters such as:

1. **Source** structure (e.g. plane wave or Gaussian focused)
2. **Source** wavelength
3. **Source** Polarization
4. **Scatterer** diameter
5. **Scatterer** refractive index
6. Medium refractive index
7. **Detector** type (photodiode or LPMode)
8. **Detector** numerical aperture
9. **Detector** angle offfset in polariation parallel axis (&phi;)
10. **Detector** angle offfset in polariation perpendicular axis (&theta;)
11. **Detector** coupling mode (Mean coupling or centered coupling)



The package also let you use a **ScattererSet** which define a range of scatterer diameter and a range of refractive index
in order to study how light scattered by such Set will be coupling in different situations.


----

Documentation
**************
All the latest available documentation is available `here <https://pymiesim.readthedocs.io/en/latest/>`_ or you can click the following badge:

|docs|

----

Google Colab
**************
It's 2021, you don't need to run all codes on you computer anymore. Google Colab is a platform which allows to write/use python script remotely.
You can open the PyMieSim.ipynb in the file (or click on the "Open in Colab" badge) to access it or click on the following badge:

|colab|

----

**Important** At the moment there is a problem with the colab compilation (it seems to be specific to boost compilation into ubuntu 18.04)


Dependencies
************
In order to install the package you first need to install some dependencies, which are the c++ `boost library <https://boost.org>`_ and some plotting library. To install one can use the command line such as:


Linux (Debian)
--------------

.. code-block:: python

   sudo apt-get install libboost-all-dev

MacOs
-----

.. code-block:: python

   brew install boost


Windows
-------
`Boost installation guide <https://www.boost.org/doc/libs/1_62_0/more/getting_started/windows.html>`_


----

Installation
************
It's pretty simple:

.. code-block:: python

   pip install PyMieSim

----

Running Unittest
*****************
To run the Unit-tests one need the coverage library.

.. code-block:: python

   python -m unittest tests/Unittest.py

----

Usage
******
Here is an example on how to use the library.

.. code-block:: python

  from PyMieSim.Source import PlaneWave
  from PyMieSim.Detector import LPmode
  from PyMieSim.Scatterer import Sphere

  LightSource = PlaneWave(Wavelength   = 450e-9,
                         Polarization = 0,
                         E0           = 1)

  Detector = LPmode(Mode         = (0, 1,'h'),
                   Sampling     = 201,
                   NA           = 0.2,
                   GammaOffset  = 0,
                   PhiOffset    = 0,
                   CouplingMode = 'Centered')


  Scat = Sphere(Diameter    = 300e-9,
               Source      = LightSource,
               Index       = 1.4)

  Coupling = Detector.Coupling(Scatterer = Scat)

  print(Coupling) # output: 1.66e+02 nWatt

For more examples I invite you to check the `examples <https://pymiesim.readthedocs.io/en/latest/Examples.html>`_
section of the documentations.


----

To-Do List
**********
- Adding T-matrix formalism
- Addind cylindrical scatterer
- Adding docstring
- Adding Stokes parameter representations
- Adding more unittests
- Adding monotonic metric to optimizer class
- Comments on c++ codes
- verify if changes of NA for <LPmode> class can be simplified (it takes way too much time)
- adding travis and codecov [DONE]



----

Citing this work?
******************
|zenodo|


----

Contact Information
************************
As of 2021 the project is still under development if you want to collaborate it would be a pleasure. I encourage you to contact me.

PyMieSim was written by `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_  .

Email:`martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=PyMieSim>`_ .



.. |codecov| image:: https://codecov.io/gh/MartinPdeS/PyMieSim/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/MartinPdeS/PyMieSim

.. |travis| image:: https://img.shields.io/travis/com/MartinPdeS/PyMieSim/master?label=Travis%20CI
   :target: https://travis-ci.com/github/numpy/numpy

.. |python| image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
   :target: https://www.python.org/

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4556074.svg
   :target: https://doi.org/10.5281/zenodo.4556074

.. |colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/drive/1FUi_hRUXxCVvkHBY10YE1yR-nTATcDei?usp=sharing

.. |docs| image:: https://readthedocs.org/projects/pymiesim/badge/?version=latest
   :target: https://pymiesim.readthedocs.io/en/latest/?badge=latest