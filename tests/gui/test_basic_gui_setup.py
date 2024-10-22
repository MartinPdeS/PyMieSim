from pytest import raises
import pytest
from unittest.mock import patch
import tkinter
from tkinter.ttk import Notebook, Frame
from PyMieSim.gui.main_window import PyMieSimGUI
from PyMieSim.gui.singleton import datashelf
import PyMieSim


def set_up_gui(foo):
    """
    This is a decorator that will set up the gui
    """
    def set_up():
        root = tkinter.Tk()
        root.geometry("750x600")
        gui = PyMieSimGUI(root)

        foo(gui=gui, root=root)

    return set_up


@set_up_gui
def test_on_close(**kwargs):
    """
    This test ensures that the .on_close method closes the tkinter window by calling a button once it has been closed and asserts that an error is raised.
    """
    gui = kwargs['gui']
    datashelf.control_tab.calculate_button.invoke()
    gui.on_close()
    with raises(tkinter.TclError):  # the gui is closed, so invoking the button should cause an error
        datashelf.control_tab.calculate_button.invoke()


@set_up_gui
def test_setup_notebook(**kwargs):
    """
    This test ensures that the setup_notebook method of PyMieSimGUI has created an instance of all the following classes.
    """
    gui = kwargs['gui']
    assert gui.component_notebook.__class__ == Notebook
    assert gui.axis_notebook.__class__ == Notebook
    assert gui.control_frame.__class__ == Frame
    assert datashelf.source_tab.__class__ == PyMieSim.gui.source_tab.SourceTab
    assert datashelf.detector_tab.__class__ == PyMieSim.gui.detector_tab.DetectorTab
    assert datashelf.scatterer_tab.__class__ == PyMieSim.gui.scatterer_tab.ScattererTab
    assert datashelf.axis_tab.__class__ == PyMieSim.gui.axis_tab.AxisTab
    assert datashelf.control_tab.__class__ == PyMieSim.gui.control_tab.ControlTab
    kwargs['root'].destroy()


@set_up_gui
@patch('matplotlib.backends.backend_tkagg.FigureCanvasTkAgg.draw')
def test_generate_figure(mock_draw, **kwargs):
    """
This test ensures that the `generate_figure` method, when called by the calculate button, does indeed invoke the canvas draw function.
    """
    datashelf.control_tab.calculate_button.invoke()
    assert mock_draw.call_count == 1, "mock_draw is not called by generate_figure"
    kwargs['root'].destroy()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
