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


def area(latitude):
    return 2 * pi * model_params.earth_radius * model_params.earth_radius * abs(
        (latitude + model_params.resol / 2).apply(lambda x: sin(x * pi / 180)) - (
                    latitude - model_params.resol / 2).apply(
            lambda x: sin(x * pi / 180))) * model_params.resol / 360  # [m²]


# Wind power calculations
# Capacity factor calculation depending on wind speed distribution and wind turbine specification
def C_f(v_r, c, k):
    return -np.exp(-pow(model_params.v_f / c, k)) + 3 * pow(c, 3) * gamma(3 / k) / (
                k * (pow(v_r, 3) - pow(model_params.v_c, 3))) * (
                       gammainc(3 / k, pow(v_r / c, k)) - gammainc(3 / k, pow(model_params.v_c / c, k)))


# Wind farm array effect = -a exp(-b * lambda) avec lambda = pi / 4*n^2
#     a5, b5 = 0.9943, 5.2661
#     a10, b10 = 0.9871, 11.7542
#     a50, b50 = 0.9838, 42.5681
#     aInf, bInt = 0.9619, 88.9204
a50, b50 = 0.9838, 42.5681


def array_effect(n):
    return a50 * np.exp(-b50 * pi / (4 * n * n))  # We assume that array size = 50x50


# Calculation of the installed Rated Power on a given Area [W], given the optimal rated wind speed and spacing parameter n
# Relationship between rated power, rotor diameter and rated wind speed
# Power_rated = 1/2 * Cp_max * rho * PI / 4 * D^2 * v_rated^3
#   => v_rated = (Power_rated / (1/2 * Cp_max * rho * PI / 4 * D^2) )^(1/3)
#   => D = (Power_rated / (1/2 * Cp_max * rho * PI / 4 * v^3) )^(1/2)
# n = turbine spacing in # rotor diameter
# Capacity density [W/m2] = Pr / (nD)^2 =  1/2 * Cp_max * rho(elevation) * PI / (4 * n^2) * v_rated^3
def rated_power(v_r, n, rho, area):
    return 1 / 2 * model_params.C_pmax * rho * pi / (4 * n * n) * pow(v_r, 3) * area


# Energy produced on a given area over wind turbine life time [J]
# Installed Power [W] * Cf [-] * array effect [-] * availability factor [-] * 25 years
def E_out_wind(v_r, n, c, k, rho, a, avail_factor):
    return rated_power(v_r, n, rho, a) * C_f(v_r, c, k) * array_effect(
        n) * avail_factor * model_params.hours_in_year * model_params.life_time_wind * model_params.watth_to_joules  # [J]


def E_out_onshore(v_r, n, c, k, rho, a):
    return E_out_wind(v_r, n, c, k, rho, a, model_params.availFactor_onshore)


def E_out_offshore(v_r, n, c, k, rho, a):
    return E_out_wind(v_r, n, c, k, rho, a, model_params.availFactor_offshore)


# Energy invested for a given available area a, based on the energy need for 1 GW wind farm [in J / GW]
def E_in_wind(v_r, n, rho, a, inputsGW):
    return rated_power(v_r, n, rho, a) * 1e-9 * inputsGW  # [J]


def E_out_solar(solar):
    return solar * model_params.eta_monoSilicon * 365 * 24 * 3.6 * model_params.life_time_solar


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
    df[name] = df[name]/total
    return df
