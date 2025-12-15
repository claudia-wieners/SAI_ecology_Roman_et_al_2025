import numpy as np
import xarray as xr
import cesm_scenarios


""" LOAD DATA """

cases_cnt = { # monthly and daily data
            'b.e21.BSSP585cmip6.f09_g17.control.01': {
    'directory': '/gpfs/work4/0/uuesm2/archive/b.e21.BSSP585cmip6.f09_g17.control.01', 'years': np.arange(2015,2100)}}

cases_sai20 = { # only monthly data
            'b.e21.BSSP585cmip6.f09_g17.2020feedback.01': {
    'directory': '/gpfs/work4/0/uuesm2/archive/b.e21.BSSP585cmip6.f09_g17.2020feedback.01', 'years': np.arange(2020,2045)},
           'b.e21.BSSP585cmip6.f09_g17.2020feedback.02': {
    'directory': '/gpfs/work4/0/uuesm2/archive/b.e21.BSSP585cmip6.f09_g17.2020feedback.02', 'years': np.arange(2045,2100)}}
           
cases_sai80 = { 'b.e21.BSSP585cmip6.f09_g17.feedback.09': { 
    'directory': '/gpfs/work4/0/uuesm2/archive/b.e21.BSSP585cmip6.f09_g17.feedback.09', 'years': np.arange(2080,2100)}}

"""
files for SAI 2080 that are not necessary right now (they only have monthly data)
        'b.e21.BSSP585cmip6.f09_g17.feedback.06': {
'directory': '/gpfs/work4/0/uuesm2/archive/b.e21.BSSP585cmip6.f09_g17.feedback.06', 'years': np.arange(2086,2100)},
        'b.e21.BSSP585cmip6.f09_g17.feedback.05': {
'directory': '/gpfs/work4/0/uuesm2/archive/b.e21.BSSP585cmip6.f09_g17.feedback.05', 'years': np.arange(2080,2086)},
"""

cases_sai20_daily = {
        'b.e21.BSSP585cmip6.f09_g17.2020feedback.02_2090': {
     'directory': '/gpfs/work4/0/uuesm2/archive/b.e21.BSSP585cmip6.f09_g17.2020feedback.02_2090', 'years': np.arange(2090,2100)}}
    
cases = [cases_cnt, cases_sai20, cases_sai80]
cases_daily = [cases_cnt, cases_sai20_daily, cases_sai80] 
scenario_names = ['Control', 'SAI 2020', 'SAI 2080']
scenarios = {name: cesm_scenarios.Scenario(name, case) for name, case in zip(scenario_names, cases)}
scenarios_daily = {name: cesm_scenarios.Scenario(name, case) for name, case in zip(scenario_names, cases_daily)}


""" FUNCTIONS """

def get_area():
    
    scenarios['Control'].get_atm_var('AREA')
    area = scenarios['Control'].var['AREA'].isel(time=0) / 10e6 # units: km2
    
    return area 


def get_T():
    
    # get dictionary with the temperature data
    for name in scenarios:
        scenarios[name].get_atm_var('TREFHT', freq = 'month')

    T1 = { name: (scenarios[name].var['TREFHT'] - 273.15) for name in scenario_names } # units: ºC
    #T1 = { name: center_lon(T0[name]) for name in scenario_names } # shift longitude so the 0 is at the center

    # concatenate SAI 2080 with Control (for plots' continuity) and create dictionary
    T = {
        'Control': T1['Control'],
        'SAI 2020': T1['SAI 2020'], # time: 957 (missing 2020-1-1)
        'SAI 2080': xr.concat((T1['Control'].sel(time = slice('2020', '2079')), T1['SAI 2080']), dim = 'time') } # time: 959
        
    return T


def get_daily_prec():
   
    # get dictionary with the DAILY precipitation data
    for name in scenarios_daily:
        scenarios_daily[name].get_atm_var('PRECC', freq = 'day')
        scenarios_daily[name].get_atm_var('PRECT', freq = 'day')

    precc = { name: (scenarios_daily[name].var['PRECC']*1000*3600*24) for name in scenario_names } # units: mm/day
    prect1 = { name: (scenarios_daily[name].var['PRECT']*1000*3600*24) for name in scenario_names } # units: mm/day
    
    #precc = { name: center_lon(precc0[name]) for name in scenario_names } # shift longitude so the 0 is at the center
    #prect1 = { name: center_lon(prect0[name]) for name in scenario_names } # shift longitude so the 0 is at the center
        
    # concatenate SAI 2080 with Control (for plots' continuity) and create dictionary
    prect = {
    'Control': prect1['Control'],
    'SAI 2020': xr.concat((prect1['Control'].sel(time = slice('2020', '2089')), prect1['SAI 2020']), dim = 'time'),
    'SAI 2080': xr.concat((prect1['Control'].sel(time = slice('2020', '2079')), prect1['SAI 2080']), dim = 'time') } 

    return precc, prect


def get_monthly_prec():

    # get dictionary with the MONTHLY precipitation data
    for name in scenarios:
        scenarios[name].get_atm_var('PRECC', freq = 'month')
        scenarios[name].get_atm_var('PRECL', freq = 'month')

    precc_m = { name: (scenarios[name].var['PRECC']*1000*3600*24) for name in scenario_names } # units: mm/day
    precl = { name: (scenarios[name].var['PRECL']*1000*3600*24) for name in scenario_names } # units: mm/day
    
    #precc_m = { name: center_lon(precc_m0[name]) for name in scenario_names } # shift longitude so the 0 is at the center
    #precl = { name: center_lon(precl0[name]) for name in scenario_names } # shift longitude so the 0 is at the center

    # concatenate SAI 2080 with Control (for plots' continuity) and create dictionary
    prect_m = {
        'Control': precl['Control'] + precc_m['Control'], # time: 960
        'SAI 2020': precl['SAI 2020'] + precc_m['SAI 2020'], # time: 957
        'SAI 2080': xr.concat(
            ((precl['Control'] + precc_m['Control']).sel(time = slice('2020', '2079')), 
             precl['SAI 2080'] + precc_m['SAI 2080']), 
            dim = 'time') } # time: 959   
    
    return precc_m, prect_m


def get_lhflx():
       
    # get dictionary with the latent heat flux
    for name in scenarios:
        scenarios[name].get_atm_var('LHFLX', freq = 'month')

    lh_flux1 = { name: scenarios[name].var['LHFLX'] for name in scenario_names }
    #lh_flux1 = { name: center_lon(lh_flux0[name]) for name in scenario_names } # shift longitude so the 0 is at the center

    # concatenate SAI 2080 with Control (for plots' continuity) and create dictionary
    lh_flux = {
        'Control': lh_flux1['Control'], # time: 960
        'SAI 2020': lh_flux1['SAI 2020'], # time: 957
        'SAI 2080': xr.concat((lh_flux1['Control'].sel(time = slice('2020', '2079')), lh_flux1['SAI 2080']), dim = 'time') } 
            # time: 959

    # negative downwards (condensation), positive upwards (evaporation)
    
    return lh_flux


def get_zonal_wind():
    
    # get dictionary with the zonal wind
    for name in scenarios:
        scenarios[name].get_atm_var('U', freq = 'month')
        
    zonal_wind_all = { name: scenarios[name].var['U'] for name in scenario_names }
    
    # select surface and 250hPa levels
    surface = { name: zonal_wind_all[name].isel(lev = -1) for name in scenario_names }
    top = { name: zonal_wind_all[name].sel(lev = '240', method = 'nearest') for name in scenario_names }
    
    # positive westerly (going to the east) winds, negative easterly (going to the west) winds
    
    return surface, top
        
    
def get_h20():
    
    # get dictionary with the water vapor concentration data
    for name in scenarios:
        scenarios[name].get_atm_var('H2O', freq = 'month')

    # select the last level because is the one that corresponds to the Earth's surface
    h20 = { name: scenarios[name].var['H2O'].isel(lev = -1) for name in scenario_names }
    
    return h20


def get_evap_ocean():
    
    for name in scenarios:
        scenarios[name].get_ocn_var('EVAP_F')

    evap = { name: (scenarios[name].var['EVAP_F'] * 3600 * 24) for name in scenario_names } # units: mm/day
    
    return evap


def get_irri():
    
    for name in scenarios:
        scenarios[name].get_lnd_var('QIRRIG')
        
    irri = { name: (scenarios[name].var['QIRRIG'] * 3600 * 24) for name in scenario_names }
    
    return irri


def get_concentration():
    
    for name in scenarios:
        scenarios[name].get_lnd_var('SO2_CMXF')
        
    con = { name: (scenarios[name].var['SO2_CMXF']) for name in scenario_names }
    
    return con


    