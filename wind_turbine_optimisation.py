import math
from scipy.optimize import minimize
import model_methods
import model_params
from world_grid import world_grid_eroi


# ---- Capacity density optimisation based on specified EROI_min ---- #
#
# Method to optimise the wind farm design in each cell, i.e. the rated speed of the wind turbines (v_r [m/s])
# and the spacing between successive rows of wind turbines in the wind farm (n [-], the distance begin n*rotor
# diameter). The results are such as the resulting wind farm as an EROI >= EROI_min (1 as default value)
#
# ---- Optimisation Inputs ---- #
#
# The optimisation inpus are either in world grid data frame, or fixed parameters of
# the model defined in class model_params
#

def capacity_density_optimisation(res_file="results/res_opti", eroi_min=1):
    grid = world_grid_eroi()
    lats = grid.Lat.array
    lon = grid.Lon.array
    c = grid.c.array
    k = grid.k.array
    # total_area = grid.Area.array
    suitable_area = grid.wind_area_onshore.array + grid.wind_area_offshore.array
    air_density = grid.air_density.array
    dissipation = grid.Dissip.array
    # Energy inputs in MWh / MW Installed (i.e. Wh/W installed)
    # ! In the dataframe energy_inputs are in J / GW installed
    energy_inputs = (grid['inputs_gw_onshore']+grid['inputs_gw_offshore']).array / model_params.watth_to_joules / 1E6 / 1000

    operation_e = ((grid['Elev'] >= 0)*model_params.oe_wind_onshore + (grid['Elev'] < 0)*model_params.oe_wind_onshore).array
    avail_factor = ((grid['Elev'] >= 0) * model_params.availFactor_onshore + (
                grid['Elev'] < 0) * model_params.availFactor_offshore).array

    n = len(lats)

    output = open(res_file, 'w')
    for i in range(0, n):
        if i % 1000 == 0:
            print("Progress ", round(((i / n) % n) * 100, 0), "%")
        output.write(str(lats[i]) + '\t' + str(lon[i]))
        if suitable_area[i] > 0:
            res = maximizeNetEnergy(eroi_min, c[i], k[i], suitable_area[i], energy_inputs[i], operation_e[i], avail_factor[i],
                                    dissipation[i], air_density[i])
            if res[2]:
                output.write('\t' + str(res[0]) + '\t' + str(res[1]))
        output.write('\n')

    output.close()


def maximizeNetEnergy(eroi_min, c, k, area, energy_inputs, operation_e, avail_factor, dissipation, air_density):
    # x = (rated wind speed, turbine spacing n)
    def net_energy(x):
        e = eroi(c, k, x[0], x[1], energy_inputs, operation_e, area, avail_factor, dissipation)
        pd = productionDensity(c, k, x[0], x[1], avail_factor, air_density)
        if (e >= eroi_min) and (pd <= 1.5 * dissipation):
            return -(area * installedCapacityDensity(x[0], x[1], air_density) * (
                    energyPerYear1MW(c, k, x[0], x[1], avail_factor) * 25 * (1 - operation_e) - energy_inputs))
        else:
            return 1000

    res = minimize(net_energy, x0=(10, 20), bounds=[(10.0, 16.0), (1, 20)])  # options={'disp': True})

    return (res.x[0], res.x[1],
            eroi(c, k, res.x[0], res.x[1], energy_inputs, operation_e, area, avail_factor, air_density) >= eroi_min)


def eroi(c, k, vr, n, energy_inputs, operation_e, area, avail_factor, air_density):
    mw = installedCapacityDensity(vr, n, air_density) * area
    out = mw * energyPerYear1MW(c, k, vr, n, avail_factor) * 25
    return out / (mw * energy_inputs + out * operation_e)


def installedCapacityDensity(vr, n, air_density, cp=0.5):
    return (0.5 * cp * air_density * math.pi / 4.0 * vr ** 3) / (n ** 2)


def energyPerYear1MW(c, k, vr, n, avail_factor):
    return model_methods.capacity_factor(vr, c, k) * model_methods.array_efficiency(n) * avail_factor * 365 * 24


def productionDensity(c, k, vr, n, avail_factor, air_density):
    return installedCapacityDensity(vr, n, air_density) * model_methods.capacity_factor(vr, c,
                                                                                        k) * model_methods.array_efficiency(
        n) * avail_factor
