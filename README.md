[![Travis CI](https://img.shields.io/travis/com/numpy/numpy/master?label=Travis%20CI)](https://travis-ci.com/github/numpy/numpy)



Welcome to PyMieSim's documentation!
====================================

PyMieSim is a tool for extensive Mie scattering analysis. It allows to study the light scattering
on different kind of object, at the moment spherical scatterer and continous sample under the Born approximation.
Using this package, one can easily set a **LightSource** a **Scatterer** and a **Detector** within a very wide range of parameters such as:
<ol>
<li>LightSource wavelength</li>
<li>LightSource Polarization</li>
<li>Scatterer diameter</li>
<li>Scatterer refractive index</li>
<li>Detector type (photodiode or LPMode)</li>
<li>Detector numerical aperture</li>
<li>Detector angle offfset in polariation parallel axis</li>
<li>Detector angle offfset in polariation perpendicular axis</li>
<li>Detector coupling mode (Mean coupling or centered coupling)</li>
</ol>


The package is also let you use a **ScattererSet** which is define by a range of scatterer diameter and a range of refractive index
in order to study how light scattered by such Set will be coupling in differents situations.


Documentation
=============
For the moment, the documentation for the package is in the Docs/build/html/index,html file.
I invite you to open it, in order to checkout some pre-defined examples.


Google Colab
============
It's 2021, you don't need to run all codes on you computer anymore. Google Colab is a platform which allows to write/use python script remotely.
You can open the PyMieSim.ipynb in the file to access it.


dependencies
============
In order to install the package you first need to install some dependencies, which are:
```console
sudo apt-get install libproj-dev proj-data proj-bin  
sudo apt-get install libgeos-dev
sudo apt-get install libboost-all-dev
```

Installation
============
```console
pip3 install -v git+https://github.com/MartinPdeS/PyMieSim.git
```

## Documentation

To-Do List
----------

- Adding docstring
- Adding monotonic metric to optimizer class
- Comments on c++ codes
- verify if changes of NA for <LPmode> class can be simplified (it takes way too much time)



Running Unittest
================

To run the Unittests one need the coverage library.

```console
   pip3 install coverage
   coverage run -m unittest Test/Unittest.py
```
