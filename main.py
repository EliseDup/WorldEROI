import plot_methods
import wind_turbine_optimisation
from world_grid import world_grid, world_rooftop_pv, world_grid_eroi

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    grid = world_grid_eroi()
    plot_methods.plot_e_out_eroi_wind(grid)


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
