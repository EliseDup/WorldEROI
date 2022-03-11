# Constants
earth_radius = 6371000  # Earth radius [m]
watth_to_joules = 3600 # Conversion from 1Wh to 1J
hours_in_year = 365*24

resol = 0.75  # Resolution of the grid, in °
# Wind turbine specifications
maxWaterDepth_wind = 1000  # [m]
v_c = 3  # Cut-in speed wind turbine [m/s]
v_f = 25  # Cut-out speed wind turbine [m/s]
C_pmax = 0.5  # Power coefficient wind turbine (theoretical maximum value is 59%)
availFactor_onshore = 0.97  # % Time wind turbine is not shut down due to maintenance, etc
availFactor_offshore = 0.95
operEnInputsOnshoreWind = 0.035  # Operational energy cost as a fraction of outputs [/]
operEnInputsOffshoreWind = 0.007  # Operational energy cost as a fraction of outputs [/]
life_time_wind = 25

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

# Solar power plants specifications
operEnInputsSolar = 0.0097  # Operational energy cost as a fraction of outputs [/]
GCR_monoSilicon = 0.2  # Ground cover ratio [%]
eta_monoSilicon = 0.186  # Solar to electricity conversion efficiency [%]
sf_commercial = 0.65
sf_residential = 0.25
life_time_solar = 25