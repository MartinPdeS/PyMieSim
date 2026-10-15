Compatibility matrix
====================

This page describes packaging targets for the documentation's source version.
The table summarizes the targets configured in ``pyproject.toml``. A configured
target is not a guarantee that an uploaded wheel exists for every release.

Python and wheel support
------------------------

.. list-table:: Configured CPython wheel targets
   :header-rows: 1
   :widths: 25 25 50

   * - Operating system
     - Architecture
     - Python versions
   * - Linux
     - x86_64
     - 3.11, 3.12, 3.13, 3.14
   * - macOS
     - arm64
     - 3.11, 3.12, 3.13, 3.14
   * - Windows
     - AMD64
     - 3.11, 3.12, 3.13, 3.14

Package metadata requires Python ``>=3.10``. Python 3.10 requires a source build
unless a matching wheel is available separately.

Linux wheels target ``manylinux_2_28``. macOS builds target Apple Silicon
with a minimum deployment version of macOS 15 and bundle their OpenMP runtime.
Windows builds target AMD64. Python versions and architectures above reflect
the configured build targets, including newly added interpreter versions.

To require a matching wheel instead of silently attempting a local compilation::

   python -m pip install --only-binary=:all: PyMieSim

If no matching distribution is found, check the interpreter version and CPU
architecture against the table before choosing a source build.

Source-build requirements
-------------------------

Building from source requires:

* Python development headers and a working C++20 compiler;
* a Fortran compiler for the Bessel subroutine;
* CMake and Ninja or Make;
* pybind11 and scikit-build-core;
* an OpenMP implementation;
* ``libomp`` on macOS, typically installed with Homebrew.

The native extensions are compiled into the ``PyMieSim`` package. A source
build is recommended for unsupported architectures or when developing the C++
backend.

Dependency policy
-----------------

PyMieSim currently pins several numerical dependencies. When diagnosing an
installation problem, inspect the resolved environment with::

   python -m pip show PyMieSim numpy matplotlib pandas
   python -m pip check

The compiled extensions must be built for the same Python implementation and
architecture as the interpreter importing them. Mixing extensions from a
different Python version, architecture, or build directory can produce import
errors such as missing symbols or incompatible binary formats.

Release verification
--------------------

Each wheel target should be verified in a clean environment by:

1. installing the wheel from outside the source checkout;
2. importing ``PyMieSim`` and every compiled extension;
3. running one minimal sphere simulation;
4. running ``python -m pip check``;
5. inspecting native dependencies with ``auditwheel``, ``delocate``, or
   ``dumpbin`` as appropriate for the platform.

See :doc:`troubleshooting` for common compiler, OpenMP, and import failures.
