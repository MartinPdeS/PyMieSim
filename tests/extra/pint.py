"""The registry cannot be replaced after native bindings capture it."""
import importlib

import pytest

from PyMieSim import _pint, units, ureg


def test_registry_initialization_is_idempotent():
    _pint.set_ureg(ureg)
    importlib.reload(units)
    assert _pint.get_ureg() is ureg is units.ureg


@pytest.mark.parametrize("replacement", [None, object()])
def test_registry_rejects_replacement_without_changing_state(replacement):
    with pytest.raises(RuntimeError, match="None|cannot replace"):
        _pint.set_ureg(replacement)
    assert _pint.get_ureg() is ureg
