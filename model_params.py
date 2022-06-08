import numpy as np

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
sf_commercial = 0.65
sf_residential = 0.25
pv_performance_ratio = 0.81

csp_gcr = 1/7.5
csp_degradation_rate = 0.2/100
csp_performance_ratio = 1.0
csp_design_irradiance = 950  # W/m^2

# ----- Wind energy inputs data ----- #
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

# ---- Solar PV : 2 technologies are available : mono or polycrystalline.
# Comment / Uncomment to select the technology you want to use in the model
# 1. Monocrystalline PV

pv_life_time_inputs = (8880746 + 5596459 + 61279 + 61279 + 469257 + 45252*pv_life_time) * 1E9  # J / GW
pv_design_efficiency = 0.24
pv_degradation_rate = 0.36/100

# 2. Polycrystalline PV
"""
pv_life_time_inputs = (12211495 + 4394479 + 71652 + 71652 + 652934 + 52911.81*pv_life_time) * 1E9  # J / GW
pv_design_efficiency = 0.17
pv_degradation_rate = 0.5 / 100
"""

# ---- Solar CSP : 3 technologies are available: here you can select the one you want by adding / removing block comments
# By default the "best" one is selected, i.e. solar tower with 12h of thermal energy storage

# Embodied energy specific calculations: transport and construction energy inputs depending on the actual size of the
# solar field

# 1. CSP Parabolic with no storage
"""
def a_csp(sm):
   return -3.38 * sm + 11.55
def b_csp(sm):
  return 23.85 * sm - 72.26
csp_life_time_inputs = (4742245 + 3178 + 114400 + 106001 + 732751 + 89118 * csp_life_time) * 1E9# J / GW
csp_design_efficiency = 0.22
csp_default_sm = 1.3
csp_transport_variable = 479183 * 1E9  # J
csp_construction_variable = (6033371 + 6732247) * 1E9 # J
sm_range = np.arange(0.5, 2.6, 0.1)
"""
# 2. CSP Parabolic with 12h of thermal storage
"""
def a_csp(sm):
    return -1.578 * sm + 11.17
def b_csp(sm):
    return 10.65 * sm - 66.33
csp_life_time_inputs = (5415779 + 13500041 + 220157 + 237600 + (756412 + 1080164) + 183720 * csp_life_time) * 1E9 # J / GW
csp_design_efficiency = 0.22
csp_default_sm = 2.7
csp_transport_variable = 930204 * 1E9  # J
csp_construction_variable = (11434969 + 4412439) * 1E9 # J
sm_range = np.arange(0.5, 4.1, 0.1)
"""

# 3. CSP Tower with 12h of thermal storage (by default)

def a_csp(sm):
    return -1.62 * sm + 8.742
def b_csp(sm):
    return 11.01 * sm - 46.86

csp_life_time_inputs = (8053825 + 9086794 + 220157 + 237600 + 1196178 + 727130 + 183720 * csp_life_time) * 1E9  # J / GW
csp_design_efficiency = 0.21
csp_default_sm = 2.7
csp_transport_variable = 457262 * 1E9  # J
csp_construction_variable = (7751182 + 3642818) * 1E9  # J
sm_range = np.arange(0.5, 4.1, 0.1)

# Common to all csp technologies
# The default aperture area (with default SM) corresponds to the value of transport and construction variable.
# These inputs are then scaled up / down based on the actual solar mutliple
csp_default_aperture_area = 1E9 * csp_default_sm / (csp_design_irradiance * csp_design_efficiency)
csp_variable_inputs = csp_transport_variable + csp_construction_variable

