

Introduction
============

This project aims to be a flexible and usefull tool for Mie scattering analysis.
At the moment the Lorenz-Mie Theory (LMT) is available and the framework for
the Generalized Lorenz-Mie Theory (GLMT) is under development.
I invite you to go to the Examples section to see what the package is capable of.


Getting started
---------------


Here is the list of packages you need to use this library:
    - Python (3+)
    - Numpy
    - Cython
    - Scipy
    - Pandas
    - Fibermodes @ git+https://github.com/cbrunet/fibermodes#egg=fibermodes-0.2.0
    - Mayavi
    - Pybind11


Those depedencies are included in the "requirements.txt" file and can be installed using the command:

.. code-block:: console
   :linenos:

   pip3 install -r requirements.txt


Installing package
------------------

First of all, the package has some c++ dependencies that can be installed using the command:

.. code-block:: console
   :linenos:

   sudo apt-get install libboost-all-dev   -> Boost library


You can now install the package using the command

.. code-block:: console
   :linenos:

    pip3 install PyMieSim


If there are some probleme with the installation I invite you o contact me, I answer quickly!


As of today (February 2021) the package was only tested on Ubuntu 20.14lts and as it
necessitate compilation of c++ core it might not be easily exportable to Windows or MacOS.
However I am currently working on the Google COLAB notebook to share with any interested
user so you won't need to install anything to use the package.



Run tests
---------

To runs tests use the following command (ubuntu 16.04).

.. code-block:: console
   :linenos:


   python3 ./tests/Unittest.py


Citing this work?
-----------------

DOI: https://doi.org/10.5281/zenodo.4556074


Author contact information
--------------------------

As of 2021 the project is still under development if you want to collaborate it would be a pleasure. I encourage you to contact me.

PyMieSim was written by `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_.

Email: `martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=PyMieSim>`_
