.. PyMieSim documentation master file, created by
   sphinx-quickstart on Wed Feb  3 17:01:20 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

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

dependencies
============
```console 
sudo apt-get install libproj-dev proj-data proj-bin  
sudo apt-get install libgeos-dev
```

Installation
============
```console
pip3 install -v git+https://gitlab.com/Martth/miecoupling.git
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
