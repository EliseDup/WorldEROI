# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import plot_methods
from world_grid import world_grid
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    df = world_grid()
    print_hi('PyCharm')
    plot_methods.plot_e_out_eroi_wind(df)
    plot_methods.plot_e_out_eroi_solar(df)


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
