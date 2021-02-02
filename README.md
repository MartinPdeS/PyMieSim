Online user's guide for the Python package for Mie Scattering analysis (PyMieAnalysis)
======================================================================

Documentation is always under development, but mostly complete.

**NOTE TO USERS:** When using PyMieAnalysis, pay close attention that all units are in S1, meaning that light wavelength and scatterer diameter is always in unit of meter.



Install PyMieAnalysis
------------------

NOTE: You must install `Shapely <https://shapely.readthedocs.io/>`_ first, preferably from GitHub. Users have reported difficulty installing it with pip. Conda works, too.

The current version is 1.8.0. You can install PyMieScatt from `The Python Package Index (PyPI) <https://pypi.python.org/pypi/PyMieScatt>`_ with ::

   $ pip install PyMieScatt


or from `GitHub <https://github.com/bsumlin/PyMieScatt>`_. Clone the repository and then run ::

   $ python setup.py install

Revision Notes - version 1.8.1 (1 February, 2021)
------------------------------------------------------------------------------

  -

Revision History
----------------


Revisions in Progress
---------------------

- Would like to re-write the inversion functions to be as general as possible, i.e., if I pass scattering, absorption, particle size, and refractive index, it would solve for the wavelength.
- Ablility to pass array objects directly to all functions (within reason).
- Auto-graphing capabilities for sacttering functions.

Documentation To-Do List
------------------------

- More example scripts, I guess?
- As a few function names and parameter names get updated, there may be some typos in old examples. I'll catch those as they crop up.

PyMieScatt To-Do List
---------------------

- Upload package to Anaconda cloud.

Publications Using PyMieAnalysis
-----------------------------

If you use PyMieScatt in your research, please let me know and I'll link the publications here.

- `Google scholar link with all citations. <https://scholar.google.com/scholar?cites=17069755164099851469&as_sdt=5,36&sciodt=0,36&hl=en>`_

- My own work using PyMieAnalysis:



Author Contact Information
--------------------------
PyMieAnalysis was written by `Martin Poinsinet de Sivry-Houle (École Polytechnique Montréal)`_.

Email: `martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=PyMieScatt>`_
