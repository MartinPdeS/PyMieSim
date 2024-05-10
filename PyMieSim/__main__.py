
import tkinter

from PyMieSim.gui.main_window import PyMieSimGUI


def main():
    root = tkinter.Tk()
    root.geometry("750x600")
    PyMieSimGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
