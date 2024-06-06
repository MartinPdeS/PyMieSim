from pytest import raises
import tkinter
import sys
sys.path.append(r'/Users/lodi/Desktop/git_project/PyMieSim')
from PyMieSim.gui.main_window import PyMieSimGUI
from unittest.mock import patch, MagicMock




@patch('numpy.savetxt')
@patch('tkinter.filedialog.asksaveasfilename')
def test_save_as_csv_button(mock_filepath, mock_save):
    """
    This test takes in a battery of widgets, and checks if the calculate button calls upon the draw method (i.e. creates a graph) for all the selections of x-axis buttons
    > At the moment, 2 widgets do not work properly
    > The testing time of this function is very long (50-60 seconds)
    """

    # setting up the environment
    root = tkinter.Tk()
    root.geometry("750x600")
    gui = PyMieSimGUI(root)

    # mocking the necessary variables
    gui.data = MagicMock()
    gui.filepath = MagicMock()

    # invoking the button
    gui.save_button.invoke()

    # the assertion
    assert mock_save.call_count == 1
