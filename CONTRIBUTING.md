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
testing. Do not hand-edit `PyMieSim/_version.py`; the release helper updates it alongside
`pyproject.toml`, the Conda recipe, and Zenodo metadata.

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

Release tags use `vMAJOR.MINOR.PATCH`. The release helper keeps the explicit
`pyproject.toml` version, Conda recipe, Zenodo version and publication date,
and tracked `PyMieSim/_version.py` aligned, following PackLab's release rules.
Do not edit `_version.py` by hand.

`make release patch`, `make release minor`, and `make release major` increment
the highest existing local semantic-version tag, create a release commit and
annotated tag, and push that commit and tag to `origin`. Fetch remote tags
before choosing a release version. The helper requires a clean working tree
and refuses an existing tag. A pushed `v*` tag triggers the shared MPSActions
`publish_release.yml@v5` workflow for PyPI and Conda, plus versioned documentation.

Run `make release-check` to check all four version declarations, or
`make release-check VERSION=vX.Y.Z` to check an expected version.

For a local release without pushing, `make tag VERSION=vX.Y.Z` (or
`make tag vX.Y.Z`) updates the release metadata and updates the source version file, creates a release commit, and makes an annotated tag;
it never pushes. Review it, then push the commit and tag explicitly.

## Shared release tooling

Install the pinned release tool in the interpreter used by Make, without
rebuilding this package:

```sh
python -m pip install "martinpdes-maintenance==0.2.2"
make release-check PYTHON=python
```

The same dependency is available through the project's `release` extra.
`tools/check_release.py`, `tools/next_release_version.py`, and
`tools/release_tag.py` are compatibility wrappers around `mps-maintenance`.
Package-specific metadata paths live in `[tool.mps-maintenance]` in
`pyproject.toml`. Release tooling is not a runtime dependency.

Preview the next release's metadata changes from a clean, full clone:

```sh
python tools/release_tag.py vX.Y.Z --dry-run
```

Fetch tags before selecting a version. Existing `make tag` and `make release`
commands keep their local-tag and push behavior respectively. The shared tool
requires an unused stable tag and a clean working tree. It stages only configured
release metadata and never pushes itself. If committing or tagging fails,
inspect `git status` and `git log -1` before recovering or retrying.

Package versions are explicit in `[project].version`. Builds do not derive
versions from Git tags or require setuptools-scm. Release commands update the
project version, package version file, Conda recipe, and Zenodo metadata together.
