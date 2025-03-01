
from PyMieSim.gui.interface import OpticalSetupGUI


if __name__ == '__main__':
    _gui = OpticalSetupGUI()
    _gui.run(host='127.0.0.1', port='8050', open_browser=True)
