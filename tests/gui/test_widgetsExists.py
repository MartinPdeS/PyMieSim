from PyMieSim.gui.main_window import PyMieSimGUI
import tkinter as tk


root = tk.Tk()
root.geometry("750x600")
gui = PyMieSimGUI(root)

def test_inputwidgetsExist():
    """
    This function is checking if the number of widgets in each tab is the expected number of widgets
    """

    assert len(gui.source_tab.widget_collection.widgets) == 4

    for tab, widget_count in zip(gui.scatterer_tab.type_widget['values'],[3,3,5]):
        gui.scatterer_tab.type_widget.set(tab)
        gui.scatterer_tab.on_type_change()
        assert len(gui.scatterer_tab.widget_collection.widgets) == widget_count

    for tab, widget_count in zip(gui.detector_tab.type_widget['values'],[5,8]):
        gui.detector_tab.type_widget.set(tab)
        gui.detector_tab.on_type_change()
        assert len(gui.detector_tab.widget_collection.widgets) == widget_count

test_inputwidgetsExist()