
import tkinter

from PyMieSim.gui.main_window import PyMieSimGUI


def main():
    root = tkinter.Tk()
    root.geometry("800x850")
    PyMieSimGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
