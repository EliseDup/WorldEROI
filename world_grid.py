import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pandas import read_csv, DataFrame
import model_params
import model_methods
from math import gamma

# Build the world grid based on input files
def world_grid():
    col_names = read_csv('data/Col_names', sep='; ', header=None)
    df = pd.read_table('data/wind_solar_0_75', header=None)
    df = df.iloc[:, 0:46]
    df = df.rename(columns=col_names.iloc[0])

    # -------- Discard cells where no RES can presumably be installed --------#
    df = df.loc[(df['Country'].isnull() == False) & (df['Country'] != 'Antarctica') & (df['Country'] != 'Greenland') & (
            df['Country'] != 'French Southern & Antarctic Lands')]
    df = df.loc[~df['Country'].apply(lambda x: 'Island' in str(x))]
    df = df.loc[~df['Country'].apply(lambda x: 'Is.' in str(x))]
    # df = df.loc[df['Elev'] >= model_params.maxWaterDepth_wind] # Why remove these ?!
    df = df.loc[df['v_r_opti'].isnull() == False]  # For some of these cells, the suppression is very debatable (e.g. in the Caspian Sea)

    df['Area'] = model_methods.area(df['Lat'])

    # -------- Compute the suitable area in each cell for RES installation --------#

    wind_sf = read_csv('data/Wind_sf', sep=',', header=None).transpose()
    df['wind_sf_onshore'] = 0
    for i in range(len(wind_sf)):
        df['wind_sf_onshore'] += df[wind_sf.iloc[0, i]] * wind_sf.iloc[1, i]
    df['wind_sf_onshore'] = df['wind_sf_onshore'] / 5625
    # LC21 = water bodies
    df['wind_sf_offshore'] = df['LC21'] / 5625

    df.loc[df['Elev'] < -model_params.maxWaterDepth_wind, 'wind_sf_offshore'] = 0
    # EU Report 4 % of 0 - 10 km, 10 % of 10 - 50 km, 25 % of > 50 km
    # NREL Report 10 % of 0 - 5 Nm, 33 % of 5 - 20 Nm, 67 % of > 20 Nm
    df.loc[df['DistCoast'] >= 37.04, 'wind_sf_offshore'] *= 0.67
    df.loc[(df['DistCoast'] < 37.04) & (df['DistCoast'] >= 9.26), 'wind_sf_offshore'] *= 0.33
    df.loc[df['DistCoast'] < 9.26, 'wind_sf_offshore'] *= 0.1

    solar_sf = read_csv('data/PV_suitability_factor', sep=',', header=None).transpose()
    df['solar_sf'] = 0
    for i in range(len(solar_sf)):
        df['solar_sf'] += df[solar_sf.iloc[0, i]] * solar_sf.iloc[1, i]
    df['solar_sf'] = df['solar_sf'] / 5625

    slope_sf = read_csv('data/Slope_suitability_factor', sep=',', header=None).transpose()
    df['slope_total'] = 0
    df['slope_sf'] = 0
    for i in range(len(slope_sf)):
        df['slope_total'] += df[slope_sf.iloc[0, i]]
        df['slope_sf'] += df[slope_sf.iloc[0, i]] * slope_sf.iloc[1, i]
    df['slope_sf'] = df['slope_sf'] / df['slope_total']

    df['wind_sf_onshore'] *= (100 - df['protected']) / 100
    df['wind_sf_offshore'] *= (100 - df['protected']) / 100
    df['solar_sf'] *= (100 - df['protected']) / 100

    df['windArea_onshore'] = df['Area'] * df['wind_sf_onshore']  # [m²]
    df['windArea_offshore'] = df['Area'] * df['wind_sf_offshore']  # [m²]
    df['solarArea'] = df['Area'] * df['solar_sf']  # [m²]

    # -------- Compute wind energy inputs per GW installed --------#
    # Fixed value for onshore / offshore bottom fixed / offshore floating
    # fixedOnshore = Gigajoules(13744075)
    # fixedOffshoreFixed = Gigajoules(18185974)
    # fixedOffshoreFloating = Gigajoules(26670974)
    # Scaling factor fixed foundations : depth < 15 = 1; depth <= 20 = 1.08; depth <= 25 = 1.34; depth <= 30 = 1.57; depth <= 35 = 1.95; depth > 35 = 2.07
    # def offshoreFixedFoundations(depth: Length) = scaling factor fixed * (Gigajoules(16173 + 361962 + 10326 + 3477293))
    df['inputsGWOnshore'] = model_params.fixedOnshore + abs(df['DistCoast']) * (model_params.onshoreOMKm + model_params.onshoreInstallationKm)

    scaling_factor_fixed_foundations = 2.19 * ((df['Elev'] > -40) & (df['Elev'] <= -35)) + 1.95 * (
                (df['Elev'] > -35) & (df['Elev'] <= -30)) + 1.57 * ((df['Elev'] > -30) & (df['Elev'] <= -25)) + 1.34 * (
                                                       (df['Elev'] > -25) & (df['Elev'] <= -20)) + 1.08 * (
                                                       (df['Elev'] > -20) & (df['Elev'] <= -15)) + 1 * (
                                                       df['Elev'] > -15)

    df['inputsGWOffshore'] = (df['Elev'] <= -40) * model_params.fixedOffshoreFixed + (df['Elev'] > -40) * model_params.fixedOffshoreFloating
    df['inputsGWOffshore'] += scaling_factor_fixed_foundations * model_params.offshoreFixedFoundations
    df['inputsGWOffshore'] += abs(df['DistCoast']) * (model_params.offshoreOMKm + model_params.offshoreInstallationKm + model_params.offshoreCableKm)

    # -------- Build functions for wind energy outputs and wind energy inputs in each cell --------#

    df['WindMean100'] = (df['WindMean71'] + df['WindMean125']) / 2
    df['WindStd100'] = (df['WindStd71'] + df['WindStd125']) / 2
    df['k'] = pow(df['WindStd100'] / df['WindMean100'], -1.086)
    df['c'] = df['WindMean100'] / (1 + 1 / df['k']).apply(lambda x: gamma(x))
    # P = 1atm * (1 - 0.0065 z / T)^5.255
    # with z = elev + wind turbine height
    df['z'] = df['Elev'] * (df['Elev'] > 0) + 100
    df['P'] = 101325 * pow(1 - 0.0065 * df['z'] / 288.15, 5.255)  # [Pa]
    # rho := air density = P/RT
    df['rho'] = df['P'] / (287.05 * 288.15)  # [kg/m³]

    df['Wind_onshore E_out [J]'] = model_methods.E_out_wind(df.v_r_opti, df.n_opti, df.c, df.k, df.rho, df.windArea_onshore, model_params.availFactor_onshore)
    df['Wind_onshore E_out [J]'] *= (1 - model_params.operEnInputsOnshoreWind)  # Subtracting operational energy inputs
    df['Wind_onshore E_out [EJ/year]'] = df['Wind_onshore E_out [J]'] * 1e-18 / model_params.life_time_wind

    df['Wind_offshore E_out [J]'] = model_methods.E_out_wind(df.v_r_opti, df.n_opti, df.c, df.k, df.rho, df.windArea_offshore, model_params.availFactor_offshore)
    df['Wind_offshore E_out [J]'] *= (1 - model_params.operEnInputsOffshoreWind)  # Subtracting operational energy inputs
    df['Wind_offshore E_out [EJ/year]'] = df['Wind_offshore E_out [J]'] * 1e-18 / model_params.life_time_wind

    df['Wind E_out [J]'] = df['Wind_onshore E_out [J]'] + df['Wind_offshore E_out [J]']
    df['Wind E_out [EJ/year]'] = df['Wind_onshore E_out [EJ/year]'] + df['Wind_offshore E_out [EJ/year]']

    df['Wind_onshore E_in'] = model_methods.E_in_wind(df.v_r_opti, df.n_opti, df.rho, df.windArea_onshore, df.inputsGWOnshore)
    df['Wind_offshore E_in'] = model_methods.E_in_wind(df.v_r_opti, df.n_opti, df.rho, df.windArea_offshore, df.inputsGWOffshore)
    df['Wind E_in'] = df['Wind_onshore E_in'] + df['Wind_offshore E_in']

    df['Wind_onshore EROI'] = df['Wind_onshore E_out [J]'] / df['Wind_onshore E_in']
    df['Wind_offshore EROI'] = df['Wind_offshore E_out [J]'] / df['Wind_offshore E_in']
    df['Wind EROI'] = df['Wind E_out [J]'] / df['Wind E_in']

    plt.plot(df['Wind EROI'], df['Wind E_out [EJ/year]'])
    plt.show()
    return df


def world_rooftop():
    col_names = read_csv('Col_names_solarRooftop', sep=', ', header=None)
    df = pd.read_table('rooftop_area', header=None)
    df = df.rename(columns=col_names.iloc[0])
    df.set_index('Country', inplace=True)
    df['Residential E_out [J]'] = df['Area PV Residential'] * model_params.sf_residential * 1e6 * df[
        'GHI'] * model_params.eta_monoSilicon * 365 * 25 * 3.6e6
    df['Residential E_out [J]'] *= (1 - model_params.operEnInputsSolar)  # Subtracting operational energy inputs
    df['Residential E_out [EJ/year]'] = df['Residential E_out [J]'] * 1e-18 / 25
    df['Residential E_in'] = df['Area PV Residential'] * model_params.sf_residential * 1e6 * 1.620032 * 1e7 * 240  # [J]

    df['Commercial E_out [J]'] = df['Area PV Commercial'] * model_params.sf_commercial * 1e6 * df[
        'GHI'] * model_params.eta_monoSilicon * 365 * 25 * 3.6e6
    df['Commercial E_out [J]'] *= (1 - model_params.operEnInputsSolar)  # Subtracting operational energy inputs
    df['Commercial E_out [EJ/year]'] = df['Commercial E_out [J]'] * 1e-18 / 25

    df['Total E_out [EJ/year]'] = df['Residential E_out [EJ/year]'] + df['Commercial E_out [EJ/year]']
    df['EROI'] = df['Residential E_out [J]'] / df['Residential E_in'].apply(
        lambda x: max(x, 1))  # EROI_residential = EROI_commercial = EROI_total

    return df
