Materials
=========

.. automodule:: PyMieSim.material
    :exclude-members: BaseMedium, BaseMaterial
    :member-order: bysource

Practical material helpers
--------------------------

PyMieSim also exposes a curated, wavelength-aware registry through
``PyMieSim.materials``::

    from PyMieSim import available_materials, load_material, load_tabulated

    print(available_materials("sellmeier"))
    glass = load_material("BK7")
    silver = load_material("silver", extrapolation="linear")

Tabulated data can be loaded from CSV or JSON. CSV files should contain
``wavelength,n,k`` columns; suffixes such as ``wavelength_nm`` and
``wavelength_um`` select the wavelength unit automatically. JSON accepts
either sample records or array fields. The default extrapolation policy is
``"error"`` and the supported interpolation policy is explicit linear
interpolation.

Use ``validate_material`` or ``validate_refractive_indices`` to check finite,
positive, passive data. PyMieSim uses the convention ``imag(n) >= 0`` for
passive materials.

