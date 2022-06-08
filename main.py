from pandas import read_csv
import numpy as np
import model_methods
import model_params
import plot_methods
from world_grid import world_grid, world_rooftop_pv, world_grid_eroi, country

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # grid = world_grid_eroi()
    for dni in range(100, 1000, 10):
        sm = model_methods.optimal_sm_csp(dni*8.76)
        area_sm = model_methods.reflective_area_csp(1E9, sm)
        e_in = (model_params.csp_life_time_inputs + area_sm / model_params.csp_default_aperture_area * model_params.csp_variable_inputs)
        print(dni, "\t", sm, "\t", model_methods.eroi_csp(dni*8.76 , sm), "\t",
              model_methods.life_time_efficiency(model_methods.efficiency_csp(dni*8.76, sm), 1.0, model_params.csp_degradation_rate,
                                   model_params.csp_life_time), "\t", model_methods.reflective_area_csp(1E9,sm), "\t", e_in

              )
    # plot_methods.geo_plot(grid, "csp_sm", "")
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
