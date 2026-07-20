#!/usr/bin/env python
"""Regression tests for dataframe unit handling."""

import warnings

import numpy as np
import pandas as pd
import pytest

from PyMieSim.experiment.dataframe_subclass import PyMieSimDataFrame
from PyMieSim.units import ureg


def test_scale_column_unit_converts_values_and_metadata():
    dataframe = PyMieSimDataFrame({"distance": [1.0, 2.0]})
    dataframe.attrs["units"] = {"distance": ureg.meter}

    converted = dataframe.scale_column_unit("distance", ureg.centimeter)

    assert converted["distance"].tolist() == [100.0, 200.0]
    assert converted.attrs["units"]["distance"] == ureg.centimeter
    assert dataframe["distance"].tolist() == [1.0, 2.0]


def test_to_compact_ignores_all_nan_columns_without_warning():
    dataframe = PyMieSimDataFrame({"distance": [np.nan, np.nan], "zero": [0.0, 0.0]})
    dataframe.attrs["units"] = {"distance": ureg.meter, "zero": ureg.meter}

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        compact = dataframe.to_compact()

    assert not caught
    assert compact["distance"].isna().all()
    assert compact["zero"].tolist() == [0.0, 0.0]


def test_scale_column_unit_rejects_unregistered_column():
    dataframe = PyMieSimDataFrame(pd.DataFrame({"distance": [1.0]}))

    with pytest.raises(ValueError, match="No unit registered"):
        dataframe.scale_column_unit("distance", ureg.centimeter)
