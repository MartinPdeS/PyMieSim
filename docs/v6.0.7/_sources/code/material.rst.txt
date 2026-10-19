Materials
=========

PyOptik 3.5.1 owns catalog discovery, file loading, metadata and physical
validation. Use ``PyMieSim.from_pyoptik`` to copy a PyOptik model into a native
solver material. Wavelength sweeps then run entirely in C++.

Catalog materials
-----------------

Initialize the RefractiveIndex.INFO catalog with ``python -m PyOptik setup``.
Use explicit PyOptik catalog IDs to select a reproducible dataset::

    import PyOptik
    from PyMieSim import from_pyoptik

    glass = from_pyoptik(
        PyOptik.load_material(name="specs/SCHOTT-optical/P-BK7")
    )
    silver = from_pyoptik(
        PyOptik.load_material(name="main/Ag/Johnson")
    )
    water = from_pyoptik(
        PyOptik.load_material(name="main/H2O/Daimon-19.0C"), medium=True
    )

These are canonical catalog IDs, including the source page, rather than
PyMieSim aliases. For example, SCHOTT P-BK7 is
``specs/SCHOTT-optical/P-BK7`` and SCHOTT N-BK7 is
``specs/SCHOTT-optical/N-BK7``; they are distinct datasets.
Use ``PyOptik.find_materials`` for discovery and inspect a loaded model's
``catalog_id`` and ``provenance`` for identity and source metadata.
The native material name retains the catalog ID.

Additional datasets used in the examples are ``main/Au/Olmon-ev`` for gold,
``main/Al/Rakic`` for aluminium and ``main/SiO2/Malitson`` for fused silica.
To use an existing snapshot elsewhere, pass ``data_root`` or a
``PyOptik.MaterialCatalog`` to ``PyOptik.load_material``.

Measured materials
------------------

Create models through PyOptik's public constructors::

    from TypedUnit import ureg

    source = PyOptik.TabulatedMaterial.from_arrays(
        name="measured", wavelength=[500, 600] * ureg.nanometer,
        n=[1.5, 1.6], k=[0.01, 0.02],
    )
    measured = from_pyoptik(source)

    source = PyOptik.TabulatedMaterial.from_csv(
        path="measurement.csv", wavelength_unit=ureg.nanometer,
    )
    measured = from_pyoptik(source)

CSV columns and missing-value policies are configured directly in PyOptik.
For example, ``missing_k=0.0`` explicitly interprets blank extinction cells as
zero. PyMieSim supplies no CSV/JSON parsers or column aliases.

Native solver constraints
-------------------------

Conversion requests strict PyOptik index validation by default. Use
``validate=False`` to skip index sign and passivity checks; finite data and
positive wavelengths remain required by PyOptik. Surrounding media must always
be lossless because their native representation stores real indices.

The solver supports linear tabulated interpolation and formula types 1, 2, 5
and 6. Separate n/k grids are exported by PyOptik on their common grid, retaining
both linear responses. Cubic interpolation and unsupported formula types are
rejected during conversion.

Native range checking defaults to strict. ``extrapolation="linear"`` permits
table endpoint-slope extrapolation or analytical formula evaluation outside
its source range. Formula index validation samples eight points within declared
bounds and does not prove validity at every wavelength.

API
---

.. autofunction:: PyMieSim.from_pyoptik

.. automodule:: PyMieSim.material
    :members:
    :undoc-members:
    :show-inheritance:
