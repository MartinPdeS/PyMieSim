Reproducible studies
====================

PyMieSim's validation and benchmark examples are executable Python sources.
Sphinx-Gallery runs those sources during the documentation build and publishes
the resulting Python files and Jupyter notebooks.

Record the environment
----------------------

For a result that can be reproduced later, record:

* PyMieSim and Python versions;
* operating system, architecture, compiler, and OpenMP runtime;
* source, scatterer, detector, material, and unit-bearing parameters;
* angular/spatial sampling and numerical settings;
* requested measures and output mode.

The benchmark example prints Python and platform metadata together with the
grid shape and timing::

   python docs/examples/benchmarks/reproducible_parameter_sweep.py

Benchmarks
----------

The benchmark scripts deliberately warm up the solver, run repeated samples,
and report the median runtime. Runtime is hardware-dependent, so compare runs
only with matching environments and configurations.

Validation
----------

The validation gallery contains comparisons against Bohren--Huffman reference
data, PyMieScatt, and internal energy-flow checks. The compact energy
conservation example can be run directly::

   python docs/examples/validation/energy_conservation.py

It asserts the identity ``Qext = Qsca + Qabs`` across a wavelength sweep.

Notebook downloads
------------------

Each Sphinx-Gallery example has downloadable Python and notebook outputs. The
generated notebooks are intended for exploration; the Python source files are
the canonical, reviewable reproducibility records.
