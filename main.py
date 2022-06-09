import model_methods
from world_grid import world_grid, world_rooftop_pv, world_grid_eroi, country

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    grid = world_grid_eroi()
    model_methods.print_potential(grid, 0)
    model_methods.print_potential(grid, 1)
    model_methods.print_potential(grid, 10)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
