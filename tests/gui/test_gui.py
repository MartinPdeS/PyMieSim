import pytest
from dash.testing.application_runners import import_app

@pytest.fixture
def app():
    """
    Fixture to set up the Dash application for testing.

    Returns
    -------
    callable
        The callable that creates the Dash app instance.
    """
    from optical_setup_gui_docstrings import OpticalSetupGUI
    return OpticalSetupGUI().app


def test_layout(dash_duo, app):
    """
    Test if the layout of the app is correctly loaded.

    Parameters
    ----------
    dash_duo : DashDuo
        Dash testing utility fixture.
    app : Dash
        The Dash application instance.
    """
    dash_duo.start_server(app)

    # Check if the title exists
    assert dash_duo.find_element("h1").text == "Optical Simulation Interface"

    # Check if the source section exists
    source_title = dash_duo.find_element("h2", text="Source")
    assert source_title is not None

    # Check if the measure section exists
    measure_title = dash_duo.find_element("h2", text="Measure")
    assert measure_title is not None


def test_callbacks(dash_duo, app):
    """
    Test if callbacks are functioning as expected.

    Parameters
    ----------
    dash_duo : DashDuo
        Dash testing utility fixture.
    app : Dash
        The Dash application instance.
    """
    dash_duo.start_server(app)

    # Simulate user interactions for generating a plot
    dropdown = dash_duo.find_element("#measure-input")
    dropdown.select_by_value("Qsca")

    generate_plot_button = dash_duo.find_element("#generate-plot")
    generate_plot_button.click()

    # Verify that a plot is generated
    plot_image = dash_duo.find_element("#plot-image")
    assert plot_image.get_attribute("src") is not None


def test_save_data(dash_duo, app):
    """
    Test if the save data functionality works correctly.

    Parameters
    ----------
    dash_duo : DashDuo
        Dash testing utility fixture.
    app : Dash
        The Dash application instance.
    """
    dash_duo.start_server(app)

    # Simulate user interactions for saving data
    filename_input = dash_duo.find_element("#filename-input")
    filename_input.send_keys("test_output.csv")

    save_button = dash_duo.find_element("#save-data")
    save_button.click()

    # Verify that the file was saved (mock test, no actual file interaction)
    assert "Data saved" in dash_duo.get_logs()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
