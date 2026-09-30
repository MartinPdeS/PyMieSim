Compatibility matrix
====================

This page describes the platforms targeted by the current packaging
configuration. A platform marked **wheel** has a compiled binary distribution;
**source** means that PyMieSim must be built locally with its native toolchain.

Python and wheel support
------------------------

.. list-table:: Current release targets
   :header-rows: 1
   :widths: 18 16 16 18 32

   * - Operating system
     - Architecture
     - Python
     - Distribution
     - Notes
   * - Linux
     - x86_64
     - 3.11
     - manylinux wheel
     - Built against ``manylinux_2_28``.
   * - Linux
     - x86_64
     - 3.12
     - manylinux wheel
     - Built against ``manylinux_2_28``.
   * - Linux
     - x86_64
     - 3.13
     - manylinux wheel
     - Built against ``manylinux_2_28``.
   * - macOS
     - arm64
     - 3.11
     - macOS wheel
     - Requires the packaged OpenMP runtime; deployment target is macOS 15.
   * - macOS
     - arm64
     - 3.12
     - macOS wheel
     - Requires the packaged OpenMP runtime; deployment target is macOS 15.
   * - macOS
     - arm64
     - 3.13
     - macOS wheel
     - Requires the packaged OpenMP runtime; deployment target is macOS 15.
   * - Windows
     - AMD64
     - 3.11
     - Windows wheel
     - Built with the configured MinGW toolchain and OpenMP support.
   * - Windows
     - AMD64
     - 3.12
     - Windows wheel
     - Built with the configured MinGW toolchain and OpenMP support.
   * - Windows
     - AMD64
     - 3.13
     - Windows wheel
     - Built with the configured MinGW toolchain and OpenMP support.
   * - Any supported OS
     - Native architecture
     - 3.10
     - Source build
     - Declared by the package metadata, but not currently included in cibuildwheel targets.

The package metadata currently declares ``Python >= 3.10``. Python 3.10 users
should therefore expect to compile from source unless a separately published
wheel is available. The project should either add Python 3.10 wheel builds or
raise ``requires-python`` to ``>=3.11`` before removing this row.

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
