import pandas as pd
from pandas import read_csv
from math import gamma

import model_params
import model_methods

# Configuration
sf_files = "data/suitability_factors/"
data_files = "data/"

# Build the world grid based on input files
def world_grid():
    col_names = read_csv(data_files+'Col_names', sep='; ', header=None)
    df = pd.read_table(data_files+'wind_solar_0_75', header=None)
    df = df.iloc[:, 0:46]
    df = df.rename(columns=col_names.iloc[0])

    # -------- Discard cells where no RES can presumably be installed --------#
    df = df.loc[(df['Country'].isnull() == False) & (df['Country'] != 'Antarctica') & (df['Country'] != 'Greenland') & (
            df['Country'] != 'French Southern & Antarctic Lands')]
    df = df.loc[~df['Country'].apply(lambda x: 'Island' in str(x))]
    df = df.loc[~df['Country'].apply(lambda x: 'Is.' in str(x))]
    # df = df.loc[df['Elev'] >= model_params.maxWaterDepth_wind] # Why remove these ?!
    df = df.loc[df['v_r_opti'].isnull() == False]  # For some of these cells, the suppression is very debatable (e.g. in the Caspian Sea)

    # -------- Compute the total area of each cell [m²] --------#
    df['Area'] = model_methods.area(df['Lat'])

    # -------- Compute the suitable area in each cell for RES installation --------#
    df = model_methods.compute_sf(df, sf_files+'wind_onshore', 'wind_sf_onshore')
    # LC21 = water bodies
    df = model_methods.compute_sf(df, sf_files + 'wind_offshore', 'wind_sf_offshore')
    df.loc[df['Elev'] < -model_params.maxWaterDepth_wind, 'wind_sf_offshore'] = 0
    # EU Report 4 % of 0 - 10 km, 10 % of 10 - 50 km, 25 % of > 50 km
    # NREL Report 10 % of 0 - 5 Nm, 33 % of 5 - 20 Nm, 67 % of > 20 Nm
    df.loc[df['DistCoast'] >= 37.04, 'wind_sf_offshore'] *= 0.67
    df.loc[(df['DistCoast'] < 37.04) & (df['DistCoast'] >= 9.26), 'wind_sf_offshore'] *= 0.33
    df.loc[df['DistCoast'] < 9.26, 'wind_sf_offshore'] *= 0.1

    df = model_methods.compute_sf(df, sf_files+'pv', 'pv_sf')
    df = model_methods.compute_sf(df, sf_files+'slope_pv', 'slope_pv_sf')
    df['pv_sf'] *= df['slope_pv_sf']

    df = model_methods.compute_sf(df, sf_files + 'csp', 'csp_sf')
    df = model_methods.compute_sf(df, sf_files + 'slope_csp', 'slope_csp_sf')
    df['csp_sf'] *= df['slope_csp_sf']

    df['wind_area_onshore'] = df['Area'] * df['wind_sf_onshore']  # [m²]
    df['wind_area_offshore'] = df['Area'] * df['wind_sf_offshore']  # [m²]
    df['pv_area'] = df['Area'] * df['pv_sf']  # [m²]
    df['csp_area'] = df['Area'] * df['csp_sf']  # [m²]

    # -------- Build parameters for wind energy in each cell --------#
    df['WindMean100'] = (df['WindMean71'] + df['WindMean125']) / 2
    df['WindStd100'] = (df['WindStd71'] + df['WindStd125']) / 2
    df['k'] = pow(df['WindStd100'] / df['WindMean100'], -1.086)
    df['c'] = df['WindMean100'] / (1 + 1 / df['k']).apply(lambda x: gamma(x))
    # P = 1atm * (1 - 0.0065 z / T)^5.255
    # with z = elev + wind_onshore turbine height
    df['z'] = df['Elev'] * (df['Elev'] > 0) + 100
    df['P'] = 101325 * pow(1 - 0.0065 * df['z'] / 288.15, 5.255)  # [Pa]
    # rho := air density = P/RT
    df['rho'] = df['P'] / (287.05 * 288.15)  # [kg/m³]

    return df


def world_grid_eroi():
    df = world_grid()
    # -------- Compute wind energy inputs per GW installed --------#
    # Fixed value for onshore / offshore bottom fixed / offshore floating
    # fixedOnshore = Gigajoules(13744075)
    # fixedOffshoreFixed = Gigajoules(18185974)
    # fixedOffshoreFloating = Gigajoules(26670974)
    # Scaling factor fixed foundations : depth < 15 = 1; depth <= 20 = 1.08; depth <= 25 = 1.34; depth <= 30 = 1.57; depth <= 35 = 1.95; depth > 35 = 2.07
    # def offshoreFixedFoundations(depth: Length) = scaling factor fixed * (Gigajoules(16173 + 361962 + 10326 + 3477293))

    df['inputs_gw_onshore'] = model_params.fixedOnshore + abs(df['DistCoast']) * (model_params.onshoreOMKm + model_params.onshoreInstallationKm)

    scaling_factor_fixed_foundations = 2.19 * ((df['Elev'] > -40) & (df['Elev'] <= -35)) + 1.95 * (
                (df['Elev'] > -35) & (df['Elev'] <= -30)) + 1.57 * ((df['Elev'] > -30) & (df['Elev'] <= -25)) + 1.34 * (
                                                       (df['Elev'] > -25) & (df['Elev'] <= -20)) + 1.08 * (
                                                       (df['Elev'] > -20) & (df['Elev'] <= -15)) + 1 * (
                                                       df['Elev'] > -15)

    df['inputs_gw_offshore'] = (df['Elev'] <= -40) * model_params.fixedOffshoreFixed + (df['Elev'] > -40) * model_params.fixedOffshoreFloating
    df['inputs_gw_offshore'] += scaling_factor_fixed_foundations * model_params.offshoreFixedFoundations
    df['inputs_gw_offshore'] += abs(df['DistCoast']) * (model_params.offshoreOMKm + model_params.offshoreInstallationKm + model_params.offshoreCableKm)

    # Energy produced [EJ/year]
    df['wind_onshore_e'] = model_methods.E_out_wind(df.v_r_opti, df.n_opti, df.c, df.k, df.rho, df.wind_area_onshore, model_params.availFactor_onshore) * 1e-18
    df['wind_offshore_e'] = model_methods.E_out_wind(df.v_r_opti, df.n_opti, df.c, df.k, df.rho, df.wind_area_offshore, model_params.availFactor_offshore) * 1e-18
    if model_params.remove_operational_e :
        df['wind_onshore_e'] *= (1 - model_params.oe_wind_onshore)
        df['wind_offshore_e'] *= (1 - model_params.oe_wind_offshore)
    df['wind_e'] = df['wind_onshore_e'] + df['wind_offshore_e']
    # Energy invested [EJ/year]
    df['wind_onshore_e_in'] = model_methods.E_in_wind(df.v_r_opti, df.n_opti, df.rho, df['wind_area_onshore'],
                                                      df['inputs_gw_onshore']) * 1e-18 / model_params.wind_life_time
    df['wind_offshore_e_in'] = model_methods.E_in_wind(df.v_r_opti, df.n_opti, df.rho, df['wind_area_offshore'],
                                                       df['inputs_gw_offshore']) * 1e-18 / model_params.wind_life_time
    df['wind_e_in'] = df['wind_onshore_e_in'] + df['wind_offshore_e_in']

    df['wind_onshore_eroi'] = df['wind_onshore_e'] / df['wind_onshore_e_in']
    df['wind_offshore_eroi'] = df['wind_offshore_e'] / df['wind_offshore_e_in']
    df['wind_eroi'] = df['wind_e'] / df['wind_e_in']

    # -------- Compute the solar pv energy outputs, energy inputs and EROI --------#

    df['pv_e'] = model_methods.E_out_solar(df['GHI'], df['pv_area']* model_params.pv_gcr) * 1e-18
    if model_params.remove_operational_e:
        df['pv_e'] *= (1 - model_params.oe_pv)
    gw_installed = model_params.wc_pv_panel * df['pv_area'] * model_params.pv_gcr / 1E9
    df['pv_e_in'] = (model_params.pv_life_time_inputs / model_params.pv_life_time) * gw_installed * 1e-18
    df['pv_eroi'] = df['pv_e'] / df['pv_e_in']

    # -------- Compute the solar csp energy outputs, energy inputs and EROI --------#
    return df

def world_rooftop_pv():
    col_names = read_csv(data_files+'Col_names_solarRooftop', sep=', ', header=None)
    df = pd.read_table(data_files+'rooftop_area', header=None)
    df = df.rename(columns=col_names.iloc[0])
    df.set_index('Country', inplace=True)
    mean_ghi = (world_grid().groupby(['Country']).mean())['GHI']
    df = pd.concat([df, mean_ghi.reindex(df.index)], axis=1)
    df.loc['Singapore', 'GHI'] = df.loc['Malaysia', 'GHI']
    df.loc['Bahrain', 'GHI'] = df.loc['Qatar', 'GHI']
    df.loc['Chinese Taipei', 'GHI'] = df.loc['China', 'GHI']
    df.loc['Hong Kong, China', 'GHI'] = df.loc['China', 'GHI']
    df.loc['Kosovo', 'GHI'] = df.loc['Montenegro', 'GHI']
    df.loc[['Netherlands Antilles', 'Gibraltar'], 'GHI'] = 0

    df['residential_e'] = model_methods.E_out_solar(df['GHI'], df['Area PV Residential'] * 1E6 * model_params.sf_residential) * 1e-18
    df['commercial_e'] = model_methods.E_out_solar(df['GHI'], df['Area PV Commercial'] * 1E6 * model_params.sf_commercial) * 1e-18
    if model_params.remove_operational_e:
        df['residential_e'] *= (1 - model_params.oe_pv)
        df['commercial_e'] *= (1 - model_params.oe_pv)

    gw_installed_res = model_params.wc_pv_panel * df['Area PV Residential'] * 1E6 * model_params.sf_residential / 1E9
    df['residential_e_in'] = (model_params.pv_life_time_inputs / model_params.pv_life_time) * gw_installed_res * 1e-18
    gw_installed_com = model_params.wc_pv_panel * df['Area PV Commercial'] * 1E6 * model_params.sf_commercial / 1E9
    df['commercial_e_in'] = (model_params.pv_life_time_inputs / model_params.pv_life_time) * gw_installed_com * 1e-18

    df['pv_e'] = df['residential_e'] + df['commercial_e']
    df['pv_e_in'] = df['residential_e_in'] + df['commercial_e_in']
    df['pv_eroi'] = df['pv_e'] / df['pv_e_in'] #.apply(lambda x: max(x, 1))  # EROI_residential = EROI_commercial = EROI_total

    return df
