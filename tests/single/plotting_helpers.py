"""Protect angular ordering and normalization conventions during plot refactors."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pytest

from PyMieSim import ureg
from PyMieSim.single.representations import FarFields, SPF, Stokes


@pytest.mark.parametrize('cls', [FarFields, SPF, Stokes])
def test_shared_plot_array_order_units_and_axis_limits(cls):
    representation = cls.__new__(cls)
    representation.sampling = 2
    values = np.arange(4) * ureg.meter
    array = representation._quantity_to_magnitude_array(values, 'centimeter')
    np.testing.assert_array_equal(representation._as_square_array(array), [[0, 200], [100, 300]])
    with pytest.raises(ValueError, match='Expected 4 values'):
        representation._as_square_array(np.arange(3))
    figure = plt.figure()
    axis = figure.add_subplot(projection='3d')
    representation._format_3d_axis(axis, 'white', True, 20, 30)
    assert axis.get_xlabel() == 'x'
    np.testing.assert_allclose(axis.get_xlim(), [-1.15, 1.15])
    assert representation._resolve_colormap('viridis', 'intensity') if cls is FarFields else representation._resolve_colormap('viridis')
    plt.close(figure)


def test_representation_orientation_is_preserved():
    far = FarFields.__new__(FarFields); far.sampling = 2
    stokes = Stokes.__new__(Stokes); stokes.sampling = 2
    spf = SPF.__new__(SPF); spf.sampling = 2; spf.SPF = np.arange(4)
    np.testing.assert_array_equal(far._complex_field_array(np.arange(4)), [[0, 1], [2, 3]])
    np.testing.assert_array_equal(stokes._stokes_array(np.arange(4)), [[0, 1], [2, 3]])
    np.testing.assert_array_equal(spf._spf_array(), [[0, 2], [1, 3]])


def test_normalization_preserves_signed_phase_and_spf_zero_clipping():
    far = FarFields.__new__(FarFields)
    stokes = Stokes.__new__(Stokes)
    spf = SPF.__new__(SPF)
    values = np.array([0., 0., 2., 4., np.nan])
    assert far._get_normalization(values, 'intensity', 'linear', 50).vmax == 3
    assert spf._get_normalization(values, 'linear', 50).vmax == 1
    log = far._get_normalization(values, 'intensity', 'log', None)
    assert isinstance(log, LogNorm)
    assert (log.vmin, log.vmax) == (2, 4)
    signed = stokes._get_normalization(np.array([-4., 2.]), None)
    assert (signed.vmin, signed.vmax) == (-4, 4)
    phase = far._get_normalization(values, 'phase', 'linear', None)
    assert (phase.vmin, phase.vmax) == (-np.pi, np.pi)
    for values in (np.zeros(3), np.full(3, np.nan)):
        signed = stokes._get_normalization(values, None)
        assert (signed.vmin, signed.vmax) == (-1, 1)
        intensity = spf._get_normalization(values, 'log', None)
        assert (intensity.vmin, intensity.vmax) == (0, 1)
