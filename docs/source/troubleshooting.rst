Troubleshooting
===============

Import and installation errors
------------------------------

Install a wheel matching your Python version and platform first::

   python -m pip install --upgrade PyMieSim
   python -c "import PyMieSim; print(PyMieSim.__version__)"

If importing a compiled extension fails, check that the interpreter running
the command is the same one where PyMieSim was installed::

   python -m pip show PyMieSim
   python -c "import sys, PyMieSim; print(sys.executable); print(PyMieSim.__file__)"

When building from source, a C++20 compiler, Fortran compiler, CMake, pybind11,
and OpenMP are required. Prefer a released wheel unless source changes require
a local build.

OpenMP and macOS
----------------

PyMieSim's native solver uses OpenMP. If CMake reports that OpenMP cannot be
found, install the platform toolchain below and rebuild with the same Python
interpreter used by your environment.

Linux
~~~~~

For Debian or Ubuntu, install the compiler, Fortran, CMake, and OpenMP
development packages before running ``make PYTHON=python editable``::

   sudo apt-get update
   sudo apt-get install build-essential gfortran cmake libomp-dev

GCC normally supplies the runtime as ``libgomp``. Clang-based builds may use
``libomp`` instead. If CMake finds the compiler but not the runtime, inspect
the configure output and ensure the selected compiler and OpenMP runtime come
from the same toolchain.

macOS
~~~~~

Install Homebrew's OpenMP runtime and a Fortran compiler::

   brew install libomp gcc cmake

On Apple Silicon, the project configuration searches Homebrew's standard
``/opt/homebrew/opt/libomp`` location. On Intel macOS, the corresponding
location is usually ``/usr/local/opt/libomp``. If you use a non-standard
installation, set ``LIBOMP_PREFIX`` to its prefix before rebuilding. Do not
reuse extensions built for another Python version or CPU architecture.

Windows
~~~~~~~

Use either Visual Studio with its OpenMP support or a MinGW-w64 toolchain with
GCC's ``libgomp``. Install CMake, a Fortran compiler such as MinGW-w64
gfortran, and the matching Python development environment. Run the build from
the corresponding developer shell so CMake sees the intended compiler. If
OpenMP is detected during configuration but import fails afterward, verify
that the OpenMP runtime DLL is on ``PATH``.

After correcting the toolchain, remove the build directory and rebuild::

   make clean
   make PYTHON=python editable

Avoid manually copying shared libraries into ``PyMieSim``; this can mix
artifacts from different builds and hide the original configuration error.

Units and input validation
--------------------------

Lengths, wavelengths, powers, and angles should carry explicit units::

   633 * ureg.nanometer
   1e-3 * ureg.watt
   0 * ureg.degree

Dimensionless refractive indices may be plain real or complex values. If a
constructor reports a unit error, inspect the quantity with ``.units`` and
convert it before passing it to PyMieSim.

Material wavelength ranges
--------------------------

Built-in and tabulated materials declare a wavelength range. The practical
material helpers default to strict range checking::

   from PyMieSim import load_material

   material = load_material("BK7")

Use ``extrapolation="linear"`` only when endpoint-slope extrapolation is
physically justified. ``validate_material`` and ``validate_wavelength`` can
be used to diagnose range and passivity problems.

Coupling errors
--------------

``coupling`` requires a detector. Add ``Photodiode``, ``CoherentMode``, or
``IntegratingSphere`` to a single simulation, or the corresponding detector
set to an experiment. Use ``available_measures`` to confirm that coupling is
enabled.

Unexpected memory use or slow sweeps
------------------------------------

Print ``array_shape`` and ``total_iterations`` before running an experiment.
Reduce the grid, request fewer measures, use ``.as_numpy()``, or process the
study in chunks. Near-field and far-field sampling can add substantial arrays;
start with a coarse sampling value.

Plots and headless environments
-------------------------------

For servers and CI, select a non-interactive Matplotlib backend before
importing plotting code::

   import matplotlib
   matplotlib.use("Agg")

Save figures explicitly instead of relying on ``show()``. Matplotlib font or
cache warnings are environment warnings and do not generally indicate a
scattering failure.

Reporting a bug
---------------

Include the PyMieSim version, Python version, operating system, installation
method, minimal inputs, requested measure, and complete traceback. For a
numerical discrepancy, include the wavelength range, material model, geometry,
and sampling settings.
