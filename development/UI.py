import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.setup import Setup
from PyMieSim.units import nanometer, degree, watt, RIU


def get_dataframe(diameter: float = 1):
    source = PlaneWave(
        wavelength=np.linspace(1000, 3000, 200) * nanometer,
        polarization=0 * degree,
        amplitude=0 * watt,
    )

    scatterer = Sphere(
        diameter=diameter * nanometer,
        property=1.5 * RIU,
        medium_property=1.0 * RIU,
        source=source
    )

    experiment = Setup(scatterer=scatterer, source=source)

    dataframe = experiment.get('Qsca')

    return dataframe


# # Streamlit app title
st.title("Basic Matplotlib UI with Streamlit")

# # Sidebar inputs for customizing the plot
st.sidebar.header("Plot Customization")
num_points = st.sidebar.slider("Number of Points", min_value=10, max_value=1000, value=100)
x_label = st.sidebar.text_input("X-axis Label", "X")
y_label = st.sidebar.text_input("Y-axis Label", "Y")
plot_color = st.sidebar.color_picker("Plot Color", "#1f77b4")
plot_title = st.sidebar.text_input("Plot Title", "Sample Plot")


diameter = st.sidebar.number_input("Scatterer diameter [nanometer]", 0)


dataframe = get_dataframe(diameter=diameter)

ax = dataframe.plot_data(x='source:wavelength', show=False)

figure = ax.get_figure()


# Display the plot in the Streamlit app
st.pyplot(figure)

# # Option to display raw data
# if st.checkbox("Show Raw Data"):
#     st.write("### Raw Data")
#     st.write({"X": x_data, "Y": y_data})

# if __name__ == "__main__":
    # from streamlit.web import cli

    # cli.main_run([f"{__file__}"])
