from PyMieSim.experiment import scatterer


class DataShelf:
    """
    This is a singleton that stores variables later accessed in another class.
    It is used to reduce coupling.

    Attributes it will hold during its lifetime:
    > x_axis_label_widget (tk.StringVar)
    > STD_axis_label_widget (tk.StringVar)
    > scatterer_tab_name (tk.StringVar)
    > measure_map (dict)
    > source_component (PyMieSim.experiment.source.Gaussian)
    > data (DataVisual.multi_array.Array)
    > figure (matplotlib.figure.Figure)
    > axis_tab (AxisTab), detector_tab (DetectorTab), source_tab (SourceTab), scatterer_tab (ScattererTab)
    """

    def __init__(self) -> None:
        # Default options for y_axis combobox
        self.measure_map = scatterer.Sphere.available_measure_list


# Instantiates the singleton. This is the only instance used throughout the rest of the code.
datashelf = DataShelf()

# -
