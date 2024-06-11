from pytest import raises
import tkinter
from tkinter.ttk import Notebook
from PyMieSim.gui.main_window import PyMieSimGUI
import PyMieSim
from unittest.mock import patch, MagicMock

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
    gui.calculate_button.invoke()
    gui.on_close()
    with raises(tkinter.TclError): # the gui is closed, so invoking the button should cause an error
        gui.calculate_button.invoke()

@set_up_gui
def test_setup_notebook(**kwargs):
    """
    This test ensures that the setup_notebook method of PyMieSimGUI has created an instance of all the following classes.
    """
    gui = kwargs['gui']
    assert gui.notebook.__class__ == Notebook
    assert gui.notebook_2.__class__ == Notebook
    assert gui.source_tab.__class__ == PyMieSim.gui.source_tab.SourceTab
    assert gui.detector_tab.__class__ == PyMieSim.gui.detector_tab.DetectorTab
    assert gui.scatterer_tab.__class__ == PyMieSim.gui.scatterer_tab.ScattererTab
    assert gui.axis_tab.__class__ == PyMieSim.gui.axis_tab.AxisTab
    kwargs['root'].destroy()
    
