# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import plot_methods
import wind_turbine_optimisation
from world_grid import world_grid, world_rooftop_pv
from model_methods import pv_efficiency
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    wind_turbine_optimisation.capacity_density_optimisation()
    #plot_methods.plot_e_out_eroi_pv(rooftop)


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
