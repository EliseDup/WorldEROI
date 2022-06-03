import math

import model_params
import numpy as np
from math import sqrt, sin, cos, pi, exp
from scipy.special import gamma, gammainc
from pandas import DataFrame, read_csv

# Methods used in the WorldEROI model

# The area of the Earth between a line of latitude and the north pole (the area of a spherical cap):
# A = 2 PI R h       with  h = R * (1-sin(lat))

# So the area between two line of latitude is:
# A = 2 PI R^2 (1 - sin(lat1)) - 2 PI R^2 (1 - sin(lat2)) = 2 PI R^2 |sin(lat1)-sin(lat2)|

# The area of the lat long rectangle is proportionnal to the difference between the 2 longitudes
# => AreaRect = 2 PI R^2 *|sin(lat1)-sin(lat2)| * |lon1 - lon2| / 360
import world_grid


def area(latitude):
    return 2 * pi * model_params.earth_radius * model_params.earth_radius * abs(
        (latitude + model_params.resol / 2).apply(lambda x: sin(x * pi / 180)) - (
                latitude - model_params.resol / 2).apply(
            lambda x: sin(x * pi / 180))) * model_params.resol / 360  # [m²]


# Wind power calculations
# Capacity factor calculation depending on wind_onshore speed distribution and wind_onshore turbine specification
def capacity_factor(v_r, c, k):
    return -np.exp(-pow(model_params.v_f / c, k)) + 3 * pow(c, 3) * gamma(3 / k) / (
            k * (pow(v_r, 3) - pow(model_params.v_c, 3))) * (
                   gammainc(3 / k, pow(v_r / c, k)) - gammainc(3 / k, pow(model_params.v_c / c, k)))


# Wind farm array effect = -a exp(-b * lambda) avec lambda = pi / 4*n^2
#     a5, b5 = 0.9943, 5.2661
#     a10, b10 = 0.9871, 11.7542
#     a50, b50 = 0.9838, 42.5681
#     aInf, bInt = 0.9619, 88.9204
a50, b50 = 0.9838, 42.5681


def array_efficiency(n):
    return a50 * np.exp(-b50 * pi / (4 * n * n))  # We assume that array size = 50x50


# Calculation of the installed Rated Power on a given Area [W], given the optimal rated wind_onshore speed and spacing parameter n
# Relationship between rated power, rotor diameter and rated wind_onshore speed
# Power_rated = 1/2 * Cp_max * rho * PI / 4 * D^2 * v_rated^3
#   => v_rated = (Power_rated / (1/2 * Cp_max * rho * PI / 4 * D^2) )^(1/3)
#   => D = (Power_rated / (1/2 * Cp_max * rho * PI / 4 * v^3) )^(1/2)
# n = turbine spacing in # rotor diameter
# Capacity density [W/m2] = Pr / (nD)^2 =  1/2 * Cp_max * rho(elevation) * PI / (4 * n^2) * v_rated^3
def rated_power(v_r, n, rho, area):
    return 1 / 2 * model_params.C_pmax * rho * pi / (4 * n * n) * pow(v_r, 3) * area


# Energy produced on a given area over wind_onshore turbine life time [J/year]
# Installed Power [W] * Cf [-] * array effect [-] * availability factor [-]
def E_out_wind(v_r, n, c, k, rho, a, avail_factor):
    return rated_power(v_r, n, rho, a) * capacity_factor(v_r, c, k) * array_efficiency(
        n) * avail_factor * model_params.hours_in_year * model_params.watth_to_joules  # [J]


def E_out_onshore(v_r, n, c, k, rho, a):
    return E_out_wind(v_r, n, c, k, rho, a, model_params.availFactor_onshore)


def E_out_offshore(v_r, n, c, k, rho, a):
    return E_out_wind(v_r, n, c, k, rho, a, model_params.availFactor_offshore)


# Energy invested for a given available area a, based on the energy need for 1 GW wind_onshore farm [in J / GW]
def E_in_wind(v_r, n, rho, a, inputs_gw):
    return rated_power(v_r, n, rho, a) * 1e-9 * inputs_gw  # [J]


# From a given irradiation in kWh/m^2/day and an suitable area in m^2, compute the annual output [J/year]
def E_out_solar(irradiation, area):
    return irradiation * 365 * 1000 * pv_efficiency() * area * model_params.watth_to_joules


def pv_efficiency():
    return life_time_efficiency(model_params.pv_design_efficiency, model_params.pv_performance_ratio,
                                model_params.pv_degradation_rate, model_params.pv_life_time)


# Return the mean efficiency based on design (=lab) efficiency, performance ratio pr, annual degradation rate [%], and total life time
def life_time_efficiency(eta, pr, degradation_rate, life_time):
    return eta * pr * (1.0 - (1.0 - degradation_rate) ** life_time) / degradation_rate / life_time


# ---- CSP Computation
# Efficiency with DNI was extrapolated using SAM simulations. For the Solar Multiple usually used, and for the Solar
# Multiple that maximizes the EROI efficiency = (a*sm + b) ln DNI
# DNI in kWh/m2/year
def efficiency_csp(dni, sm):
    return (model_params.a_csp(sm) * math.log(dni) + model_params.b_csp(sm)) / 100.0


# Rated power of the power block of a CSP plant depending on a given solar multiple
# The power block is designed to provide its nominal power at a design irradiance (in general 950W/m^2)
# Then the collector (reflective) area is over / under scaled based on local conditions (via the solar multiple)
# Rated power [W] = Reflective Area [m^2] * Design Irradiance [W/m^2] * Design Efficiency [%] / SM
def rated_power_csp(collector_area, sm):
    return collector_area * (model_params.csp_design_irradiance * model_params.csp_design_efficiency) / sm


# Inverse relationship
# Reflective Area [m^2] = Rated power [W] / (Design Irradiance [W/m^2] * Design Efficiency [%] / SM)
def reflective_area_csp(rated_power, sm):
    return rated_power / (model_params.csp_design_irradiance * model_params.csp_design_efficiency / sm)


# Prevent the model to produce at CF > 100% (it can happen due to the extrapolation model with extreme values of DNI and SM)
def e_out_csp(area, dni, sm):
    e_out_max = rated_power_csp(area, sm) * 365 * 24 * model_params.watth_to_joules
    e_out = np.minimum(e_out_max, area * dni * life_time_efficiency(efficiency_csp(dni, sm), 1.0, model_params.csp_degradation_rate,
                                                  model_params.csp_life_time) * model_params.watth_to_joules * 1000)
    return e_out


# Return the EROI of a CSP plant (technology defined in model_params) (DNI in kWH/m2/year)
def eroi_csp(dni, sm):
    area_sm = reflective_area_csp(1E9, sm)
    area_ratio = area_sm/model_params.csp_default_aperture_area
    e_out = e_out_csp(area_sm, dni, sm)
    e_in = e_out * model_params.oe_csp + (model_params.csp_life_time_inputs + area_ratio * model_params.csp_variable_inputs)/model_params.csp_life_time
    return e_out/e_in



# Return the solar multiple that maximises the EROI for a given DNI in kWH/m2/year
def optimal_sm_csp(dni):
    if dni == 0:
        return model_params.csp_default_sm
    else:
        return model_params.sm_range[np.where(np.amax(eroi_csp(dni, model_params.sm_range)))[0]]


# Build cumulated E out [EJ/year] and EROI dataframe, based on the world grid dataframe df,
# And the name of the corresponding eout et eroi column in the world grid df
def df_cum_eout_eroi(df, eout, eroi):
    df_cum = DataFrame(index=df.index)
    df_cum[['e', 'eroi']] = df.loc[:, [eout, eroi]]
    df_cum = df_cum.sort_values(by=['eroi'], ascending=False)
    df_cum['e_cum'] = df_cum['e'].cumsum()
    return df_cum


# Data processing
# From a table of (land cover name, suitability factor), build a corresponding suitability factor per cell
# Inputs:   sf_table = file location of the table
#           name = of the new column created in the dataframe
#
def compute_sf(df, sf_table, name):
    sf = read_csv(sf_table, sep=',', header=None).transpose()
    df[name] = 0
    total = 0
    for i in range(len(sf.iloc[0, :])):
        df[name] += df[sf.iloc[0, i]] * sf.iloc[1, i]
        total += df[sf.iloc[0, i]]
    df[name] = df[name] / total
    # Correct the suitability factor to account for the proportion of protected areas in each cell
    df[name] *= (100 - df['protected']) / 100
    return df


def print_results_country(countries, grid):
    unit = 1  # model_params.ej_to_twh
    print("Country", "Area[1000km2]", "WindOnshore[EJ/y]", "InstalledCapacity[GW]", "CF[%]", "SuitableArea[km2]",
          "WindOffshore[EJ/y]", "InstalledCapacity[GW]", "CF[%]",
          "SuitableArea[1000km2]", "PV[EJ/y]", "InstalledCapacity[GW]", "CF[%]", "SuitableArea[1000km2]")
    for c in countries:
        cells = world_grid.country(c, grid)
        if len(cells) > 0:
            print(c, cells["Area"].sum() / 1E9, cells["wind_onshore_e"].sum() * unit, cells["wind_onshore_gw"].sum(),
              cells["wind_onshore_e"].sum() * model_params.ej_to_twh / (
                          cells["wind_onshore_gw"].sum() * 365 * 24 / 1000),
              cells["wind_area_onshore"].sum() / 1E9,
              cells["wind_offshore_e"].sum() * unit, cells["wind_offshore_gw"].sum(),
              cells["wind_offshore_e"].sum() * model_params.ej_to_twh / (
                          cells["wind_offshore_gw"].sum() * 365 * 24 / 1000),
              cells["wind_area_offshore"].sum() / 1E9,
              cells["pv_e"].sum() * unit, cells["pv_gw"].sum(),
              cells["pv_e"].sum() * model_params.ej_to_twh / (cells["pv_gw"].sum() * 365 * 24 / 1000),
              cells["pv_area"].sum() / 1E9)
