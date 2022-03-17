import math
from scipy.optimize import minimize
import model_methods
import model_params
from world_grid import world_grid_eroi

# ---- Optimisation Inputs ---- #
# Lat, lon, c, k, total area, suitable area wind, energy inputs pr MW [MWh], operationnal energy [% output], avail factor


def capacity_density_optimisation(eroi_min = 1):
    grid = world_grid_eroi()
    lats = grid.Lat.array
    lon = grid.Lon.array
    c = grid.c.array
    k = grid.k.array
    total_area = grid.Area.array
    suitable_area = grid.wind_area_onshore.array
    energy_inputs = grid['inputs_gw_onshore'].array / model_params.watth_to_joules / 1E6 / 1000
    # Fraction of the energy outputs
    operation_e = model_params.oe_wind_onshore
    avail_factor = model_params.availFactor_onshore
    air_density = grid.rho.array
    dissipation = grid.Dissip.array

    n = len(lats)

    output = open('results/res_opti', 'w')
    for i in range(0, n):
        if i % 1000 == 0:
            print(i)
        output.write(str(lats[i]) + '\t' + str(lon[i]))
        if suitable_area[i] > 0:
            res = maximizeNetEnergy(eroi_min, c[i], k[i], suitable_area[i], energy_inputs[i], operation_e, avail_factor, dissipation[i], air_density[i])
            if res[2]:
                output.write('\t' + str(res[0]) + '\t' + str(res[1]))
        output.write('\n')

    output.close()


def maximizeNetEnergy(eroi_min, c, k, area, energy_inputs, operation_e, avail_factor, dissipation, air_density):
    # x = (rated wind speed, turbine spacing n)
    def net_energy(x):
        if ((eroi(c, k, x[0], x[1], energy_inputs, operation_e, area, avail_factor, dissipation) >= eroi_min) & (
                productionDensity(c, k, x[0], x[1], avail_factor, air_density) <= 1.0)):  # <= 1.5*dissipation)): #1.0)):
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
    return out / (mw*energy_inputs+out*operation_e)


def installedCapacityDensity(vr, n, air_density, cp=0.5):
    return (0.5*cp*air_density*math.pi/4.0*vr**3) / (n**2)


def energyPerYear1MW(c, k, vr, n, avail_factor):
    return model_methods.capacity_factor(vr, c, k)*model_methods.array_efficiency(n)*avail_factor*365*24


def productionDensity(c, k, vr, n, avail_factor, air_density):
    return installedCapacityDensity(vr, n, air_density)*model_methods.capacity_factor(vr, c, k)*model_methods.array_efficiency(n)*avail_factor
