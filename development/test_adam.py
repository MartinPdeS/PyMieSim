import tkinter as tk
from tk import *

win = tk.Tk()

win.geometry('500x500')

button_var = "1"

class Input:
    
    def __init__(self, win) -> None:
        self.win = win

    def set_up(self, row):
        self.label = tk.Label(self.win, text="input")
        self.label.grid(row = row, column = 0)
        self.inp = tk.Entry(self.win)
        self.inp.grid(row = row, column = 1)
        self.button1 = tk.Radiobutton(self.win, variable= button_var, value=row, command= self.get_entry())
        self.button1.grid(row= row, column= 2)
        
    def get_entry(self):
        print(self.inp.get())
        
widgets = []

for i in range(2):
    A = Input(win)
    A.set_up(i)
    widgets.append(A)


win.mainloop()