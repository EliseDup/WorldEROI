import model_methods
import model_params
import tif_processing
from world_grid import world_grid, world_rooftop_pv, world_grid_eroi, country

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    world = world_grid_eroi()
    model_methods.print_results_country(["Australia", "Chile", "Oman", "Spain", "Morocco"], world)
    print("-----")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
