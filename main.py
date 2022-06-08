from pandas import read_csv
import numpy as np
import model_methods
import model_params
import plot_methods
from world_grid import world_grid, world_rooftop_pv, world_grid_eroi, country

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    grid = world_grid_eroi()

    plot_methods.geo_plot(grid, "wind_eroi", "")


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
