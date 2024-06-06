from pytest import raises
import tkinter
from PyMieSim.gui.main_window import PyMieSimGUI
from unittest.mock import patch, MagicMock

def reset_std_button(possible_widgets, gui):

    for widget in possible_widgets:

        widget = gui.source_tab.widget_collection.widgets[0]
        widget.tk_radio_button_2.invoke()
        gui.reset_std_button.invoke()

        assert gui.STD_axis_label_widget.get() == 'None', f"reset button did not worl for {widget.tk_radio_button_2['value']} widget"

@patch('matplotlib.backends.backend_tkagg.FigureCanvasTkAgg.draw')
def x_axis_calculate_button(mock_draw, possible_widgets, gui):
    # Run the tests
    error_button_values = ['polarization_value', 'polarization_label'] #this is a list of the values of all the buttons that for which the calculate button does not work

    for x_widget in possible_widgets:
        x_widget.tk_radio_button_1.invoke()
        mock_draw.reset_mock()
        gui.calculate_button.invoke()
        assert mock_draw.call_count == 1, f"{x_widget.tk_radio_button_1['value']} button did not call the mock"

def test_in_all_combination_of_widgets():
    root = tkinter.Tk()
    root.geometry("750x600")
    gui = PyMieSimGUI(root)

    """
    This function is ment to cycle trough all combinations of tabs and generate the corresponding battery of widgets with which the tests should be ran
    It will then pass the tests for all tab combination
    """

    # This combination of for loops will create all the possible widget combinations (this is slow. Check if it is best to separate the function in multiple fonctions without the loops)
    source_widgets = gui.source_tab.widget_collection.widgets[0:2]

    for scatter_tab in gui.scatterer_tab.type_widget['values']:
        gui.scatterer_tab.type_widget.set(scatter_tab)
        gui.scatterer_tab.on_type_change()
        scatter_widgets = gui.scatterer_tab.widget_collection.widgets

        for detector_tab in gui.detector_tab.type_widget['values'][:-1]:
            gui.detector_tab.type_widget.set(detector_tab)
            gui.detector_tab.on_type_change()
            all_detector_widgets = gui.detector_tab.widget_collection.widgets
            detector_widgets =[widget for widget in all_detector_widgets if widget.component_label != "mean_coupling"]

            possible_widgets = [*source_widgets, *scatter_widgets, *detector_widgets] # The widgets it is possible to choose from for the std and x axis
            
            reset_std_button(possible_widgets=possible_widgets, gui=gui)
            x_axis_calculate_button(possible_widgets=possible_widgets, gui=gui)

    root.destroy()


@patch('tkinter.filedialog.asksaveasfilename')
@patch('tkinter.messagebox.showinfo')
def test_export_plot_button(mock_messagebox, mock_filepath):
    """
    This function is testing if the export_plot_button function is working as intended
    """
    
    # setting up the environment
    root = tkinter.Tk()
    root.geometry("750x600")
    gui = PyMieSimGUI(root)

    # mocking the necessary variables
    gui.figure = MagicMock()
    gui.filepath = MagicMock()

    # invoking the button
    gui.export_button.invoke()

    # the assertion
    assert gui.figure.savefig.call_count == 1

    root.destroy()


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
    
    root.destroy()




# -------------------------- LEGACY -------------------------
'''
@patch('PyMieSim.gui.main_window.PyMieSimGUI.save_data_as_csv')
def x_axis_save_as_csv_button(mock_save, widgets : list):
    """
    This test takes in a battery of widgets and checks if the save as CSV button calls upon the save_data_as_csv method for all selections of x-axis
    > This test is not ideal. Ideally, it should mock the np.savetxt method, but that requires to also mock two conditions to be true, and I don't know how to do that
    """

    # Initialize the root and the gui
    root1 = tkinter.Tk()
    root1.geometry("750x600")
    gui1 = PyMieSimGUI(root)

    # The Test
    for widget in widgets:
                mock_save.reset_mock()
                widget.tk_radio_button_1.invoke()
                gui1.save_button.invoke()
                assert mock_save.call_count == 1
'''