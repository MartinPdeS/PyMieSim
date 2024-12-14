import pytest
import chromedriver_autoinstaller
from PyMieSim.simple_app import create_app

@pytest.fixture
def app():
    """
    Fixture to set up the Dash application for testing.

    Returns
    -------
    callable
        The Dash app instance.
    """
    return create_app()

def test_callbacks(dash_duo, app):
    """
    Test if the callback updates the output correctly.

    Parameters
    ----------
    dash_duo : DashDuo
        Dash testing utility fixture.
    app : Dash
        The Dash app instance.
    """
    # Ensure ChromeDriver is installed and available
    chromedriver_autoinstaller.install()

    dash_duo.start_server(app)

    # Simulate entering text
    input_box = dash_duo.find_element("#input-text")
    input_box.send_keys("Hello, Dash!")

    # Simulate clicking the submit button
    dash_duo.find_element("#submit-button").click()

    # Verify the output text
    dash_duo.wait_for_text_to_equal("#output-div", "You entered: Hello, Dash!", timeout=4)
