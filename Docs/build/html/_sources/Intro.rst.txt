Introduction
============

This project (SuPyModes) aims to create a python library that gives the necesary tool to simulate propagation mode of light inside fibers coupler.


Getting started
---------------


Here is the list of packages you need to use this library:
    - Python (3+)
    - Numpy
    - Cython
    - Scipy
    - Pandas
    - Fibermodes: https://github.com/cbrunet/fibermodes
    - Mayavi


Those depedencies are included in the "requirements.txt" file


Installing package
------------------

Soon enough you will be able to use "pip" to install PyMieSim but for the moment one can install it manually
from the `github repository <https://github.com/MartinPdS/PyMieSim>`_


Once downloaded the class command "python setup.py install" should do the trick. If not i invite you o contact me, I answer quickly!


As of today (February 2021) the package was only tested on Ubuntu 20.14lts and as it
necessitate compilation of c++ core it might not be easily exportable to Windows or MacOS.
However I am currently working on the Google COLAB notebook to share with any interested
user so you wont need to install anything to use the package.



Run test simulation
-------------------

Here an example of command to run a simulation on linux (ubuntu 16.04).

.. code-block:: console
   :linenos:
