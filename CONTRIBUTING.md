# Contributing to PyMieSim

Thank you for improving PyMieSim. Focused bug fixes, validation cases,
documentation corrections, and feature contributions are welcome.

## Development setup

PyMieSim requires Python 3.10 or newer, CMake, a C++ compiler, pybind11, and
a Fortran compiler. On macOS, install the compiler and OpenMP runtime required
by the local CMake configuration before building.

```bash
git clone https://github.com/MartinPdeS/PyMieSim.git
cd PyMieSim
make editable
make test
```

After changing native sources, rebuild the editable installation before
testing. Do not hand-edit `PyMieSim/_version.py`; it is generated from Git
metadata during releases and builds.

## Repository layout

`PyMieSim/` contains the importable Python package. Native C++ and Fortran
sources, including pybind11 bindings, live in the top-level `cpp/` directory.
CMake installs compiled extensions into their corresponding package locations;
`PyMieSim._pint` is the private runtime unit-binding extension.

## Contribution expectations

- Add tests for every behaviour change. Validation changes should compare to a
  documented analytical, numerical, or experimental reference where possible.
- Exercise public Python APIs rather than only private extension modules.
- Keep units explicit at public boundaries and document physical assumptions,
  conventions, and numerical limits.
- Add NumPy-style docstrings for Python APIs and useful docstrings for bound
  native classes.
- Keep documentation examples small enough for the gallery build.
- Do not commit generated documentation, native build products, caches, or
  Python bytecode.

Run the relevant focused tests first, then the project checks before opening a
pull request:

```bash
make quality
make test
make release-check
```

Native regression tests cover NumPy copy ownership and error propagation from
parallel far-field sweeps. Enable them when configuring a development build:

```bash
make configure
cmake -S . -B build -DPYMIESIM_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Pull requests and releases

Keep pull requests focused. Explain the scientific or technical motivation,
the public behaviour affected, and the verification performed. Report platform
and compiler details for native changes.

Release tags use `vMAJOR.MINOR.PATCH`. `make tag VERSION=vX.Y.Z` regenerates
the source version file, creates a release commit, and makes an annotated tag;
it never pushes. Review it, then push the commit and tag explicitly.
