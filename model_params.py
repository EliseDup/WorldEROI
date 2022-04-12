# Constants
earth_radius = 6371000  # Earth radius [m]
watth_to_joules = 3600  # Conversion from 1Wh to 1J
ej_to_twh = 277.778 # Conversion from 1EJ to 1TWh
hours_in_year = 365*24

resol = 0.75  # Resolution of the grid, in °
remove_operational_e = True

# RES life times
wind_life_time = 25
pv_life_time = 25
csp_life_time = 30

# Operational energy cost as a fraction of outputs [/]
oe_wind_onshore = 0.035
oe_wind_offshore = 0.007
oe_pv = 0.0097
oe_csp = 0.05 + 0.023

# Wind turbine specifications
maxWaterDepth_wind = -1000  # [m]
v_c = 3  # Cut-in speed wind_onshore turbine [m/s]
v_f = 25  # Cut-out speed wind_onshore turbine [m/s]
C_pmax = 0.5  # Power coefficient wind_onshore turbine (theoretical maximum value is 59%)
availFactor_onshore = 0.97  # % Time wind_onshore turbine is not shut down due to maintenance, etc
availFactor_offshore = 0.95

# Solar power plants specifications
wc_pv_panel = 240 # Wc / m2
pv_gcr = 1/5.0 # Ground cover ratio [%]
# eta_mono_silicon = 0.186  # Mean solar to electricity conversion efficiency over monosilicon pv lifetiem [%]
sf_commercial = 0.65
sf_residential = 0.25
eta_mono_si = 0.24
eta_poly_si = 0.17
pv_degradation_rate = 0.5/100
pv_performance_ratio = 0.81

csp_gcr = 1/7.5
csp_degradation_rate = 0.2/100
csp_performance_ratio = 1.0

# ----- Energy inputs data ----- #
# Wind
fixedOnshore = 13744075 * 1E9  # J / GW PJ : why 1.2850855e16 ?
fixedOffshoreFixed = 18185974 * 1E9  # J / GW PJ : why 1.7502782e16 ?
fixedOffshoreFloating = 26670974 * 1E9  # J / GW
offshoreFixedFoundations = (16173 + 361962 + 10326 + 3477293) * 1E9  # J / GW, cost of offshore fixed foundations, vayring with water depth
# Inputs per km (cables for offshore, distance to shore for onshore)
offshoreInstallationKm = 16904 * 1E9  # J / GW / km
offshoreCableKm = (4681 + 105) * 1E9  # J / GW / km
offshoreOMKm = 6615 * 1E9  # J / GW / km
onshoreInstallationKm = 605.74 * 1E9  # J / GW / km
onshoreOMKm = 21.3 * 1E9  # J / GW / km

# Solar PV
pv_life_time_inputs = (8880746 + 5596459 + 61279 + 61279 + 469257 + 45252) * 1E9 # J / GW



