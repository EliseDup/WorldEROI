import model_params
import numpy as np
from math import sqrt, sin, cos, pi, exp
from scipy.special import gamma, gammainc

# The area of the Earth between a line of latitude and the north pole (the area of a spherical cap):
# A = 2 PI R h       with  h = R * (1-sin(lat))

# So the area between two line of latitude is:
# A = 2 PI R^2 (1 - sin(lat1)) - 2 PI R^2 (1 - sin(lat2)) = 2 PI R^2 |sin(lat1)-sin(lat2)|

# The area of the lat long rectangle is proportionnal to the difference between the 2 longitudes
# => AreaRect = 2 PI R^2 *|sin(lat1)-sin(lat2)| * |lon1 - lon2| / 360


def area(latitude):
    return 2 * pi * model_params.earth_radius * model_params.earth_radius * abs((latitude + model_params.resol/2).apply(lambda x: sin(x * pi/180)) - (latitude - model_params.resol/2).apply(lambda x: sin(x * pi/180))) * model_params.resol / 360  # [m²]


# Wind power calculations
# Capacity factor calculation depending on wind speed distribution and wind turbine specification
def C_f(v_r, c, k):
    return -np.exp(-pow(model_params.v_f / c, k)) + 3 * pow(c, 3) * gamma(3 / k) / (k * (pow(v_r, 3) - pow(model_params.v_c, 3))) * (gammainc(3 / k, pow(v_r / c, k)) - gammainc(3 / k, pow(model_params.v_c / c, k)))


# Wind farm array effect = -a exp(-b * lambda) avec lambda = pi / 4*n^2
#     a5, b5 = 0.9943, 5.2661
#     a10, b10 = 0.9871, 11.7542
#     a50, b50 = 0.9838, 42.5681
#     aInf, bInt = 0.9619, 88.9204
def array_effect(n):
    return 0.983825 * np.exp(-42.568*pi/(4*n*n)) # We assume that array size = 50x50


# Wind turbine rated power [W]
def rated_power(v_r, n, rho, a):
    return model_params.C_pmax * pi * rho * pow(v_r, 3) / (8 * n * n) * a


# Energy porduced over life time [J]
def E_out_wind(v_r, n, c, k, rho, a, avail_factor):
    return C_f(v_r, c, k) * array_effect(n) * rated_power(v_r, n, rho, a) * avail_factor * 3600 * 24 * 365 * model_params.life_time_wind  # [J]


def E_out_onshore(v_r, n, c, k, rho, a):
    return E_out_wind(v_r, n, c, k, rho, a, model_params.availFactor_onshore)


def E_out_offshore(v_r, n, c, k, rho, a):
    return E_out_wind(v_r, n, c, k, rho, a, model_params.availFactor_offshore)


def E_in_wind(v_r, n, rho, a, inputs):
    return rated_power(v_r, n, rho, a) * 1e-9 * inputs  # [J]


def E_out_solar(solar):
    return solar * model_params.eta_monoSilicon * 365 * 24 * 3.6 * model_params.life_time_solar