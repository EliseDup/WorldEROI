# -*- coding: utf-8 -*-
"""
Computation of global wind and solar EROI curves

Pierre JACQUES
April 2021
"""

# This python script recomputes from scratch the global wind and solar EROI curves, initally built by Elise Dupont
# Only the optimal values in each cell of v_r and n for EROI>=1 are given as input instead of being recomputed via an optimization problem, like Elise did

import numpy as np
from scipy.optimize import minimize, NonlinearConstraint
import matplotlib.pyplot as plt
import pandas as pd
from math import sqrt,sin,cos,gamma,pi, exp
from scipy.special import gammainc
from numpy import concatenate
from numpy import arange
from numpy import dot
from numpy import mean
from pandas import read_csv
from pandas import DataFrame
from pandas import concat


# Parameters values
maxWaterDepth_wind = -1000
R = 6371000 # Earth radius [m]
resol = 0.75
v_c = 3 # [m/s]
v_f = 25 # [m/s]
C_pmax = 0.5
availFactor_onshore = 0.97
availFactor_offshore = 0.95
operEnInputsOnshoreWind = 0.035 # Fraction of outputs [/]
operEnInputsOffshoreWind = 0.007 # Fraction of outputs [/]
operEnInputsSolar = 0.0097 # Fraction of outputs [/]
GCR_monoSilicon = 0.2
eta_monoSilicon = 0.186


#-------- Read data file --------#

col_names = read_csv('data/Col_names', sep='; ', header=None)
df = pd.read_table('data/wind_solar_0_75', header=None)
df = df.iloc[:,0:46]
df = df.rename(columns=col_names.iloc[0])
#df.drop(columns=['DNI'], inplace=True)


#-------- Discard cells where no RES can presumably be installed --------#

df = df.loc[(df['Country'].isnull() == False) & (df['Country'] != 'Antarctica') & (df['Country'] != 'Greenland') & (df['Country'] != 'French Southern & Antarctic Lands')]
df = df.loc[~df['Country'].apply(lambda x: 'Island' in str(x))]
df = df.loc[~df['Country'].apply(lambda x: 'Is.' in str(x))]
df = df.loc[df['Elev']>=maxWaterDepth_wind]
df = df.loc[df['v_r_opti'].isnull()==False] # For some of these cells, the suppression is very debatable (e.g. in the Caspian Sea)


#-------- Compute the total area of each cell --------#

# The area of the Earth between a line of latitude and the north pole (the area of a spherical cap):
# A = 2 PI R h       with  h = R * (1-sin(lat))

# So the area between two line of latitude is:
# A = 2 PI R^2 (1 - sin(lat1)) - 2 PI R^2 (1 - sin(lat2)) = 2 PI R^2 |sin(lat1)-sin(lat2)|

# The area of the lat long rectangle is proportionnal to the difference between the 2 longitudes
# => AreaRect = 2 PI R^2 *|sin(lat1)-sin(lat2)| * |lon1 - lon2| / 360

df['Area'] = pi * R*R * abs( (df['Lat']+resol/2).apply(lambda x: sin(x*pi/180)) - (df['Lat']-resol/2).apply(lambda x: sin(x*pi/180)) ) * resol  / 180  # [m²]


#-------- Compute the suitable area in each cell for RES installation --------#

df_windSuita = read_csv('data/Wind_suitability_factor', sep=',', header=None).transpose()
df_windSuita = df_windSuita.rename(columns=df_windSuita.iloc[0])
df_windSuita.drop([0], axis=0, inplace=True)
df['windSuitaFactor_onshore'] = df.iloc[:,12:35].dot(df_windSuita.transpose()) / 5625
df['windSuitaFactor_offshore'] = df.iloc[:,32] / 5625
df.loc[df['DistCoast']>=37.04, 'windSuitaFactor_offshore'] *= 0.67
df.loc[(df['DistCoast']<37.04) & (df['DistCoast']>=9.26), 'windSuitaFactor_offshore'] *=  0.33
df.loc[df['DistCoast']<9.26, 'windSuitaFactor_offshore'] *=  0.1

df_solarSuita = read_csv('data/PV_suitability_factor', sep=',', header=None).transpose()
df_solarSuita = df_solarSuita.rename(columns=df_solarSuita.iloc[0])
df_solarSuita.drop([0], axis=0, inplace=True)
df['solarSuitaFactor'] = df.iloc[:,12:35].dot(df_solarSuita.transpose()) / 5625

df_slopeSuita = read_csv('Slope_suitability_factor', sep=',', header=None).transpose()
df_slopeSuita = df_slopeSuita.rename(columns=df_slopeSuita.iloc[0])
df_slopeSuita.drop([0], axis=0, inplace=True)
df['slopeSuitaFactor'] = 1
df_slopeIntermediary = df.iloc[:,3:11].dot(df_slopeSuita.transpose())
df['slopeTotal'] = df.iloc[:,3:11].sum(axis=1)
df['Intermed'] = df_slopeIntermediary
df.loc[df['slopeTotal']>0, 'slopeSuitaFactor'] = df_slopeIntermediary.loc[df['slopeTotal']>0, 1] / df.loc[df['slopeTotal']>0, 'slopeTotal']
df['solarSuitaFactor'] *= df['slopeSuitaFactor']

df['windSuitaFactor_onshore'] *= (100-df['protected'])/100
df['windSuitaFactor_offshore'] *= (100-df['protected'])/100
df['solarSuitaFactor'] *= (100-df['protected'])/100

df['windArea_onshore'] = df['Area'] * df['windSuitaFactor_onshore']  # [m²]
df['windArea_offshore'] = df['Area'] * df['windSuitaFactor_offshore']  # [m²]
df['solarArea'] = df['Area'] * df['solarSuitaFactor']  # [m²]


#-------- Compute wind energy inputs per GW installed --------#

df['inputsOnshore'] = 1.2850855e16 + (6.05739375e11 + 2.13e10) * abs(df['DistCoast']) # J/GW
df['inputsOffshore'] = 1.7502782e16 + 8.485e15 * (df['Elev']<=-40) # J/GW
df['inputsOffshore'] += 3.656429e15 * ( 2.19 * ((df['Elev']>-40)&(df['Elev']<=-35)) + 1.95 * ((df['Elev']>-35)&(df['Elev']<=-30)) + 1.57 * ((df['Elev']>-30)&(df['Elev']<=-25)) + 1.34 * ((df['Elev']>-25)&(df['Elev']<=-20)) + 1.08 * ((df['Elev']>-20)&(df['Elev']<=-15)) + 1 * (df['Elev']>-15) )
df['inputsOffshore'] += 2.652527e13 * abs(df['DistCoast'])


#-------- Build functions for wind energy outputs and wind energy inputs in each cell --------#

df['WindMean100'] = (df['WindMean71'] + df['WindMean125']) / 2
df['WindStd100'] = (df['WindStd71'] + df['WindStd125']) / 2
df['k'] = pow( df['WindStd100']/df['WindMean100'], -1.086 )
df['c'] = df['WindMean100'] / (1+1/df['k']).apply(lambda x: gamma(x))
# P = 1atm * (1 - 0.0065 z / T)^5.255
# with z = elev + wind turbine height
df['z'] = df['Elev']*(df['Elev']>0) + 100
df['P'] = 101325 * pow(1-0.0065*df['z']/288.15, 5.255) # [Pa]
# rho := air density = P/RT
df['rho'] = df['P'] / (287.05 * 288.15) # [kg/m³]

def C_f(v_r, c, k):
    return -exp(-pow(v_f/c, k)) + 3*pow(c,3)*gamma(3/k) / (k*(pow(v_r,3)-pow(v_c, 3))) * ( gammainc(3/k, pow(v_r/c,k)) - gammainc(3/k, pow(v_c/c,k)) )

def eta(n):
    return 0.983825 * np.exp(-42.568*pi/(4*n*n)) # We assume that array size = 50x50

def MW(v_r, n, rho, A):
    return C_pmax*pi*rho*pow(v_r,3) / (8*n*n) * A

def E_out_onshore(v_r, n, c, k, rho, A_Onshore, A_Offshore):
    return C_f(v_r, c, k) * eta(n) * MW(v_r, n, rho, A_Onshore) * availFactor_onshore * 3600*24*365*25 # [J]

def E_out_offshore(v_r, n, c, k, rho, A_Onshore, A_Offshore):
    return C_f(v_r, c, k) * eta(n) * MW(v_r, n, rho, A_Offshore) * availFactor_offshore * 3600*24*365*25 # [J]

def E_in_onshore(v_r, n, rho, A_Onshore, A_Offshore, inputOnshore, inputOffshore):
    return MW(v_r, n, rho, A_Onshore) * 1e-9 * inputOnshore # [J]

def E_in_offshore(v_r, n, rho, A_Onshore, A_Offshore, inputOnshore, inputOffshore):
    return MW(v_r, n, rho, A_Offshore) * 1e-9 * inputOffshore # [J]


#-------- Example optimization problem with constraints --------#

# def fun(x,y,z):
#     return (x[0]-y)*(x[0]-y) + (x[1]-z)*(x[1]-z)

# def cons_f(x):
#     return x[0]*x[0] + x[1] - 1

# cnstr = NonlinearConstraint(cons_f, 0, 0)

# def opti_fun(y):
#     return minimize(fun, x0=[1000,1000], args=(y[0],y[1]), constraints=cnstr).x

# d = {'col1': [-1,2,3], 'col2': [3,4,7]}
# df_opti = pd.DataFrame(data=d)
# df_opti['sol']= df_opti.apply(opti_fun, axis=1)


#-------- Optimization problem to recompute v_r_opti and n_opti for wind --------#

# def Net_E_out(var, c, k, rho, A_Onshore, A_Offshore, inputOnshore, inputOffshore):
#     return - ( E_out(var[0], var[1], c, k, rho, A_Onshore, A_Offshore) - E_in(var[0], var[1], rho, A_Onshore, A_Offshore, inputOnshore, inputOffshore) )
    
# def opti_fun(x):
#     return minimize(Net_E_out, x0=[11.4,13], args=(x[0],x[1],x[2],x[3],x[4],x[5],x[6]))
# # Still need to add four constraints:
# # 1) 10 <= v_r <= 16
# # 2) 1 <= n <= 20
# # 3) EROI >= 1
# # 4) OUT <= Dissipation [W/m²]

# df_opti = DataFrame(index=df.index)
# df_opti[['c', 'k', 'rho', 'windArea_onshore', 'windArea_offshore', 'inputsOnshore', 'inputsOffshore']] = df.loc[:,['c', 'k', 'rho', 'windArea_onshore', 'windArea_offshore', 'inputsOnshore', 'inputsOffshore']]
# df_opti['sol'] = df_opti.apply(opti_fun, axis=1)
# df['v_r_opti_PJ'] = df_opti['sol'].apply(lambda y: y.x[0])
# df['n_opti_PJ'] = df_opti['sol'].apply(lambda y: y.x[1])


#-------- Given v_r_opti and n_opti, compute wind energy outputs, energy inputs and EROI --------#

df_compute_E_out = DataFrame(index=df.index)
df_compute_E_out[['v_r', 'n', 'c', 'k', 'rho', 'windArea_onshore', 'windArea_offshore']] = df.loc[:,['v_r_opti', 'n_opti', 'c', 'k', 'rho', 'windArea_onshore', 'windArea_offshore']]
def apply_E_out_onshore(x):
    return E_out_onshore(x[0], x[1], x[2], x[3], x[4], x[5], x[6])       

def apply_E_out_offshore(x):
    return E_out_offshore(x[0], x[1], x[2], x[3], x[4], x[5], x[6])

df['Wind_onshore E_out [J]'] = df_compute_E_out.apply(apply_E_out_onshore, axis=1)
df['Wind_onshore E_out [J]'] *= (1-operEnInputsOnshoreWind) # Subtracting operational energy inputs
df['Wind_onshore E_out [EJ/year]'] = df['Wind_onshore E_out [J]'] * 1e-18 / 25

df['Wind_offshore E_out [J]'] = df_compute_E_out.apply(apply_E_out_offshore, axis=1)
df['Wind_offshore E_out [J]'] *= (1-operEnInputsOffshoreWind)
df['Wind_offshore E_out [EJ/year]'] = df['Wind_offshore E_out [J]'] * 1e-18 / 25

df['Wind E_out [J]'] = df['Wind_onshore E_out [J]'] + df['Wind_offshore E_out [J]']
df['Wind E_out [EJ/year]'] = df['Wind_onshore E_out [EJ/year]'] + df['Wind_offshore E_out [EJ/year]']

df_compute_E_in = DataFrame(index=df.index)
df_compute_E_in[['v_r', 'n', 'rho', 'windArea_onshore', 'windArea_offshore', 'inputsOnshore', 'inputsOffshore']] = df.loc[:,['v_r_opti', 'n_opti', 'rho', 'windArea_onshore', 'windArea_offshore', 'inputsOnshore', 'inputsOffshore']]
def apply_E_in_onshore(x):
    return E_in_onshore(x[0], x[1], x[2], x[3], x[4], x[5], x[6])

def apply_E_in_offshore(x):
    return E_in_offshore(x[0], x[1], x[2], x[3], x[4], x[5], x[6])

df['Wind_onshore E_in'] = df_compute_E_in.apply(apply_E_in_onshore, axis=1)
df['Wind_offshore E_in'] = df_compute_E_in.apply(apply_E_in_offshore, axis=1)
df['Wind E_in'] = df['Wind_onshore E_in'] + df['Wind_offshore E_in']

df['Wind_onshore EROI'] = df['Wind_onshore E_out [J]'] / df['Wind_onshore E_in']
df['Wind_offshore EROI'] = df['Wind_offshore E_out [J]'] / df['Wind_offshore E_in']
df['Wind EROI'] = df['Wind E_out [J]'] / df['Wind E_in']


#-------- Compute the solar energy outputs, energy inputs and EROI --------#

df['Solar E_out [J]'] = df['solarArea'] * df['GHI'] * GCR_monoSilicon * eta_monoSilicon * 365 * 25 * 3.6e6
# same for Leccisi except that eta_monoSilicon = 0.205 and LT = 30 instead of 25
df['Solar E_out [J]'] *= (1 - operEnInputsSolar) # Subtracting operational energy inputs
df['Solar E_out [EJ/year]'] = df['Solar E_out [J]'] * 1e-18 / 25

df['Solar E_in'] = 1.620032 * 1e7 * 240 * df['solarArea'] * GCR_monoSilicon # [J]
# same for Leccisi except that Inputs = 1.3615 instead of 1.620032 and Wp = 205 instead of 240

df['Solar EROI'] = df['Solar E_out [J]'] / df['Solar E_in'].apply(lambda x: max(x,1))


#-------- Compute the rooftop solar energy output per country --------#

col_names_solarRooftop = read_csv('data/Col_names_solarRooftop', sep=', ', header=None)
rooftop_PV = pd.read_table('data/rooftop_area', header=None)
rooftop_PV = rooftop_PV.rename(columns=col_names_solarRooftop.iloc[0])
rooftop_PV.set_index('Country', inplace=True)
mean_GHI = (df.groupby(['Country']).mean())['GHI']
rooftop_PV = pd.concat([rooftop_PV, mean_GHI.reindex(rooftop_PV.index)], axis=1)
rooftop_PV.loc['Singapore','GHI'] = rooftop_PV.loc['Malaysia','GHI']
rooftop_PV.loc['Bahrain','GHI'] = rooftop_PV.loc['Qatar','GHI']
rooftop_PV.loc['Chinese Taipei','GHI'] = rooftop_PV.loc['China','GHI']
rooftop_PV.loc['Hong Kong, China','GHI'] = rooftop_PV.loc['China','GHI']
rooftop_PV.loc['Kosovo','GHI'] = rooftop_PV.loc['Montenegro','GHI']
rooftop_PV.loc[['Netherlands Antilles', 'Gibraltar'],'GHI'] = 0

rooftop_PV['Residential E_out [J]'] = rooftop_PV['Area PV Residential'] * 0.25 * 1e6 * rooftop_PV['GHI'] * eta_monoSilicon * 365 * 25 * 3.6e6
rooftop_PV['Residential E_out [J]'] *= (1 - operEnInputsSolar) # Subtracting operational energy inputs
rooftop_PV['Residential E_out [EJ/year]'] = rooftop_PV['Residential E_out [J]'] * 1e-18 / 25
rooftop_PV['Residential E_in'] = rooftop_PV['Area PV Residential'] * 0.25 * 1e6 * 1.620032 * 1e7 * 240 # [J]

rooftop_PV['Commercial E_out [J]'] = rooftop_PV['Area PV Commercial'] * 0.65 * 1e6 * rooftop_PV['GHI'] * eta_monoSilicon * 365 * 25 * 3.6e6
rooftop_PV['Commercial E_out [J]'] *= (1 - operEnInputsSolar) # Subtracting operational energy inputs
rooftop_PV['Commercial E_out [EJ/year]'] = rooftop_PV['Commercial E_out [J]'] * 1e-18 / 25

rooftop_PV['Total E_out [EJ/year]'] = rooftop_PV['Residential E_out [EJ/year]'] + rooftop_PV['Commercial E_out [EJ/year]']
rooftop_PV['EROI'] = rooftop_PV['Residential E_out [J]'] / rooftop_PV['Residential E_in'].apply(lambda x: max(x,1)) # EROI_residential = EROI_commercial = EROI_total



#-------- Plot Wind EROI curve --------#

df_plot_wind_onshore = DataFrame(index=df.index)
df_plot_wind_onshore[['E_out', 'EROI']] = df.loc[:,['Wind_onshore E_out [EJ/year]', 'Wind_onshore EROI']]
df_plot_wind_onshore = df_plot_wind_onshore.sort_values(by=['EROI'], ascending=False)
df_plot_wind_onshore['cum_E_out'] = df_plot_wind_onshore['E_out'].cumsum()

df_plot_wind_offshore = DataFrame(index=df.index)
df_plot_wind_offshore[['E_out', 'EROI']] = df.loc[:,['Wind_offshore E_out [EJ/year]', 'Wind_offshore EROI']]
df_plot_wind_offshore = df_plot_wind_offshore.sort_values(by=['EROI'], ascending=False)
df_plot_wind_offshore['cum_E_out'] = df_plot_wind_offshore['E_out'].cumsum()

df_plot_wind = DataFrame(index=df.index)
df_plot_wind[['E_out', 'EROI']] = df.loc[:,['Wind E_out [EJ/year]', 'Wind EROI']]
df_plot_wind = df_plot_wind.sort_values(by=['EROI'], ascending=False)
df_plot_wind['cum_E_out'] = df_plot_wind['E_out'].cumsum()

plt.figure()
plt.plot(df_plot_wind['cum_E_out'], df_plot_wind['EROI'], label='Total')
plt.plot(df_plot_wind_onshore['cum_E_out'], df_plot_wind_onshore['EROI'], label='Onshore')
plt.plot(df_plot_wind_offshore['cum_E_out'], df_plot_wind_offshore['EROI'], label='Offshore')
plt.grid(True, color="#93a1a1", alpha=0.3)
plt.legend(loc='upper right', fancybox=True, shadow=True)
plt.xlabel(r'Cumulated Production [EJ/year]', labelpad=8, fontsize=12)
plt.ylabel(r'EROI [/]', labelpad=6, fontsize=12)


#-------- Plot Solar EROI curve --------#

df_plot_solar = DataFrame(index=df.index)
df_plot_solar[['E_out', 'EROI']] = df.loc[:,['Solar E_out [EJ/year]', 'Solar EROI']]
df_plot_solar_rooftop = DataFrame(index=rooftop_PV.index)
df_plot_solar_rooftop[['E_out', 'EROI']] = rooftop_PV.loc[:,['Total E_out [EJ/year]', 'EROI']]
df_plot_solar = concat([df_plot_solar, df_plot_solar_rooftop])
df_plot_solar.drop(df_plot_solar[df_plot_solar['EROI']<=0].index, inplace=True)
df_plot_solar = df_plot_solar.sort_values(by=['EROI'], ascending=False)
df_plot_solar['cum_E_out'] = df_plot_solar['E_out'].cumsum()

plt.figure()
plt.plot(df_plot_solar['cum_E_out'], df_plot_solar['EROI'])
plt.grid(True, color="#93a1a1", alpha=0.3)
plt.xlabel(r'Cumulated Production [EJ/year]', labelpad=8, fontsize=12)
plt.ylabel(r'EROI [/]', labelpad=6, fontsize=12)


#-------- Plot Wind and Solar EROI curves --------#

plt.figure()
plt.plot(df_plot_wind['cum_E_out'], df_plot_wind['EROI'], label='Wind')
plt.plot(df_plot_solar['cum_E_out'], df_plot_solar['EROI'], label='Solar')
plt.grid(True, color="#93a1a1", alpha=0.3)
plt.legend(loc='upper right', fancybox=True, shadow=True)
plt.xlabel(r'Cumulated Production  [EJ/year]', labelpad=8, fontsize=12)
plt.ylabel(r'EROI [/]', labelpad=6, fontsize=12)


#-------- Plot curve of RES capital intensity (aka gamma_e) --------#

df_plot_total1 = DataFrame(index=df.index)
df_plot_total1[['E_out', 'EROI']] = df.loc[:,['Wind_onshore E_out [EJ/year]', 'Wind_onshore EROI']]
df_plot_total2 = DataFrame(index=df.index)
df_plot_total2[['E_out', 'EROI']] = df.loc[:,['Wind_offshore E_out [EJ/year]', 'Wind_offshore EROI']]
df_plot_total3 = DataFrame(index=df.index)
df_plot_total3[['E_out', 'EROI']] = df.loc[:,['Solar E_out [EJ/year]', 'Solar EROI']]
df_plot_total4 = DataFrame(index=rooftop_PV.index)
df_plot_total4[['E_out', 'EROI']] = rooftop_PV.loc[:,['Total E_out [EJ/year]', 'EROI']]
df_plot_total = concat([df_plot_total1, df_plot_total2, df_plot_total3, df_plot_total4])
df_plot_total = df_plot_total.sort_values(by=['EROI'], ascending=False)
df_plot_total['cum_E_out'] = df_plot_total['E_out'].cumsum()
df_plot_total.reset_index(drop=True, inplace=True)
df_plot_total.drop(df_plot_total.iloc[77518:].index, inplace=True)

# Transform marginal EROI curve into global EROI curve
df_plot_total['E_in'] = df_plot_total['E_out'] / df_plot_total['EROI']
df_plot_total['Global EROI'] = df_plot_total['E_out'].cumsum() / df_plot_total['E_in'].cumsum()

df_plot_total['Cap Intensity'] = 1/df_plot_total['Global EROI']
df_plot_total['Cap Intensity'] *= (25/0.005068)

plt.figure()
plt.plot(df_plot_total['cum_E_out'], df_plot_total['Global EROI'])
plt.grid(True, color="#93a1a1", alpha=0.3)
plt.xlabel(r'Cumulated Production [EJ/year]', labelpad=8, fontsize=12)
plt.ylabel(r'global EROI  [/]', labelpad=6, fontsize=12)


plt.figure()
plt.plot(df_plot_total['cum_E_out'], df_plot_total['Cap Intensity'])
# plt.xlim([0,2000])
# plt.ylim([0,1800])
# plt.xticks(np.arange(0, max(df_plot_total['cum_E_out'])+1, 400))
# plt.yticks(np.arange(0, max(df_plot_total['Cap Intensity'])+1, 400))
plt.grid(True, color="#93a1a1", alpha=0.3)
plt.xlabel(r'Cumulated Production [EJ/year]', labelpad=8, fontsize=12)
plt.ylabel(r'gamma_e   [US$2017/(GJ/year)]', labelpad=6, fontsize=12)
# plt.savefig('/home/pjacques/Desktop/EROI_curve_Wind_Solar', format='pdf')


#-------- Approximate the curve of gamma_e by a 5th order polynomial --------#

pol5 = np.polyfit(df_plot_total['cum_E_out'].astype(float).values, df_plot_total['Cap Intensity'].astype(float).values, 5)
pol5_vals = np.polyval(pol5, df_plot_total['cum_E_out'].astype(float).values)

plt.figure()
plt.plot(df_plot_total['cum_E_out'], df_plot_total['Cap Intensity'], label='gamma_e')
plt.plot(df_plot_total['cum_E_out'], pol5_vals, label='5th order approx.')
plt.legend(loc='lower right', fancybox=True, shadow=True)
plt.grid(True, color="#93a1a1", alpha=0.3)
plt.xlabel(r'Cumulated Production [EJ/year]', labelpad=8, fontsize=12)
plt.ylabel(r'Cap intensity   [US$2018/(GJ/year)]', labelpad=6, fontsize=12)
# plt.savefig('/home/pjacques/Desktop/PhD/Useful images/Latex_Global energy-economy model/gamma_re.pdf', format='pdf')

# #-------- Generate text file for colour plot on World map --------#

df_colourPlot = DataFrame(index=df.index)
df_colourPlot[[0,1]] = df.loc[:,['Lat','Lon']]
df_colourPlot['Solar EROI'] = df.loc[:,'Solar EROI']

data_5 = np.genfromtxt("technology_5", delimiter='\t', dtype=None)
data_5 = pd.DataFrame(data_5)

df_colourPlot = pd.merge(data_5, df_colourPlot, how="left", on=[0,1])
df_colourPlot.drop(columns=2,inplace=True)
df_colourPlot['Solar EROI'] = df_colourPlot['Solar EROI'].fillna(0)

df_colourPlot.to_csv('Solar_EROI', sep='\t', index=False, header=False)


#-------- Plot for the article --------#

# Transform marginal EROI curve into global EROI curve
df_plot_wind_onshore['E_in'] = df_plot_wind_onshore['E_out'] / df_plot_wind_onshore['EROI']
df_plot_wind_offshore['E_in'] = df_plot_wind_offshore['E_out'] / df_plot_wind_offshore['EROI']
df_plot_solar['E_in'] = df_plot_solar['E_out'] / df_plot_solar['EROI']
df_plot_wind_onshore['Global EROI'] = df_plot_wind_onshore['E_out'].cumsum() / df_plot_wind_onshore['E_in'].cumsum()
df_plot_wind_offshore['Global EROI'] = df_plot_wind_offshore['E_out'].cumsum() / df_plot_wind_offshore['E_in'].cumsum()
df_plot_solar['Global EROI'] = df_plot_solar['E_out'].cumsum() / df_plot_solar['E_in'].cumsum()

plt.figure()
plt.plot(df_plot_total['cum_E_out'], df_plot_total['Global EROI'], label='Total', color='tab:blue')
plt.plot(df_plot_wind_onshore['cum_E_out'], df_plot_wind_onshore['Global EROI'], label='Onshore wind', color='tab:green')
plt.plot(df_plot_wind_offshore['cum_E_out'], df_plot_wind_offshore['Global EROI'], label='Offshore wind', color='tab:purple')
plt.plot(df_plot_solar['cum_E_out'], df_plot_solar['Global EROI'], label='Solar', color='tab:orange')
plt.grid(True, color="#93a1a1", alpha=0.3)
plt.legend(loc='upper right', fancybox=True, shadow=True)
plt.xlabel(r'Final energy production  [EJ/year]', labelpad=8, fontsize=12)
plt.ylabel(r'EROI [/]', labelpad=6, fontsize=12)
plt.savefig('/home/pjacques/Desktop/PhD/Useful images/Latex_Global energy-economy model/EROI_curve_Wind_Solar_forPaper.pdf', format='pdf')


# Fossil EROI curve
A = 6.5
B = 0.00002
cum_y_nre = arange(1,140000)
EROI_fossil = 1 + A*np.exp(-B*cum_y_nre)

plt.figure()
plt.plot(cum_y_nre, EROI_fossil, color='tab:blue')
plt.grid(True, color="#93a1a1", alpha=0.3)
plt.xlabel(r'Cumulated final energy production [EJ]', labelpad=8, fontsize=12)
plt.ylabel(r'EROI [/]', labelpad=6, fontsize=12)
plt.savefig('/home/pjacques/Desktop/PhD/Useful images/Latex_Global energy-economy model/EROI_curve_Fossil_forPaper.pdf', format='pdf')

