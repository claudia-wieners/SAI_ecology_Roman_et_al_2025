import numpy as np
import xarray as xr
import cesm_scenarios
import matplotlib.pyplot as plt
from matplotlib import gridspec as gs


""" GENERAL MASKS """

def lon_mask(da, lon_initial, lon_final):
    
    """
    Creates a longitudinal mask given a certain condition (1 condition fullfilled, 0 otherwise)

    Parameters:
        da (DataArray): data array containing a field defined over latitudes (lat) and longitudes (lon) 
            and optionally time
        lon_initial (float): initial longitude of the considered region, starting from -180º
        lon_final (float): final longitude of the considered region, ending at 180º

    Returns:
        Data array with the same shape as the original one but containing only 1 or 0 
    """
    
    da_ones = xr.ones_like(da)
    
    if (lon_initial < 0) & (lon_final > 0):
        mask_lon = da_ones.where( (da.lon >= 360 + lon_initial) & (da.lon < 360) | 
                                 (da.lon > 0) & (da.lon <= lon_final), 0 )
    elif (lon_initial < 0) & (lon_final < 0):
        mask_lon = da_ones.where( (da.lon >= 360 + lon_initial) & (da.lon <= 360 + lon_final), 0 )
    else:
        mask_lon = da_ones.where( (da.lon >= lon_initial) & (da.lon <= lon_final), 0 )
        
    return mask_lon


def create_mask(lat_initial, lat_final, mask_lon):
    
    """
    Creates a global mask for a regular rectangle given a certain condition (1 condition fullfilled, 0 otherwise)

    Parameters:
        mask_lon (DataArray): data array containing a longitudinal mask defined over latitudes (lat) 
            and optionally time
        lat_initial (float): initial latitude of the considered region, starting from -90N
        lat_final (float): final latitude of the considered region, ending at 90N

    Returns:
        Data array with the same shape as the original one but containing only 1 or 0
    """
    
    mask = mask_lon.where( (mask_lon.lat >= lat_initial) & (mask_lon.lat <= lat_final), 0 )
        
    return mask


def land_mask():

    """
    Gets the land mask by looking at the LANDFRAC variable in the atm module
    """
    
    cases_cnt = {'b.e21.BSSP585cmip6.f09_g17.control.01': {'directory':
                                                           '/gpfs/work4/0/uuesm2/archive/b.e21.BSSP585cmip6.f09_g17.control.01',
                                                           'years': np.arange(2020,2100) }}

    cases = [cases_cnt]
    scenario_names = ['Control']
    scenarios = {name: cesm_scenarios.Scenario(name, case) for name, case in zip(scenario_names, cases)}

    for name in scenarios:
        scenarios[name].get_atm_var('LANDFRAC', freq = 'month')

    landmask_zero = { name: scenarios[name].var['LANDFRAC'] for name in scenario_names }
    
    # change 0 to nan so they are not plotted
    land_mask = landmask_zero['Control'].where(landmask_zero['Control'].data == 1, np.nan) 
    
    return land_mask


def ocean_mask():
    
    """
    Gets the ocean mask by looking at the OCNFRAC variable in the atm module
    """
    
    cases_cnt = {'b.e21.BSSP585cmip6.f09_g17.control.01': {'directory':
                                                           '/gpfs/work4/0/uuesm2/archive/b.e21.BSSP585cmip6.f09_g17.control.01',
                                                           'years': np.arange(2020,2100) }}

    cases = [cases_cnt]
    scenario_names = ['Control']
    scenarios = {name: cesm_scenarios.Scenario(name, case) for name, case in zip(scenario_names, cases)}

    for name in scenarios:
        scenarios[name].get_atm_var('OCNFRAC', freq = 'month')

    oceanmask_zero = { name: scenarios[name].var['OCNFRAC'] for name in scenario_names }

    # change 0 to nan so they are not plotted
    ocean_mask = oceanmask_zero['Control'].where(oceanmask_zero['Control'].data == 1, np.nan) 
    
    return ocean_mask



""" MEAN """

def global_mean(da, mask = None, ref = None, trend = 'yearly'):

    """
    Takes an xarray DataArray defined over a lat/lon field over time
    and calculates the global/regional annual/seasonal mean computed with correct weight factors.
    Can also work on an exclusively meridional field (only lat coordinates)

    Parameters:
        da (DataArray): data array containing a field defined over latitudes (lat) and longitudes (lon) and time
        mask (DataArray): data array containing the corresponding mask
        ref (DataArray): data array containing a time series for the reference value
            ONLY for the regional mean (remember it has to be the reference over the region!)
        trend (string): trend to calculate, either yearly or seasonal (specify the season)

    Returns:
        Data array with the global or regional annual mean of the field. If plotted, it is a time series
    """
    
    assert trend in ['yearly', 'season_DJF', 'season_MAM', 'season_JJA', 'season_SON'], f'trend must be either yearly, season_DJF, season_MAM, season_JJA, or season_SON and not {trend}'
    
    # apply mask
    if not mask is None: # if mask is not None
        da = da * mask.isel(time = 0)
    
    # detect if data has longitudinal component
    if 'lon' in list(da.coords):
        zonal = False
    else:
        zonal = True
    
    # determine the grid size
    # note: this requires an equally spaced coordinate grid
    dlat = np.deg2rad(np.diff(da.lat.data)[0])
    dlon = np.deg2rad(np.diff(da.lon.data)[0])
    
    # determine normalization factors
    if mask is None:
        norm_factor = 2
    else:
        mask_zonal = mask.sum(dim = "lon") * dlon
        mask_int = ((np.cos(np.deg2rad(da.lat))) * dlat * mask_zonal).sum(dim = "lat")
        norm_factor = float(mask_int[0])
    
    # zonal integration
    if not zonal: # if there is longitudinal component
        da_zonal = da.sum(dim = "lon") * dlon
        if mask is None: # apply zonal normalization if no mask is given
            da_zonal *= 1.0 / (2.0 * np.pi)
        else:
            da_zonal = da_zonal
    else: 
        da_zonal = da

    # meridional integration
    weight = 1.0 * np.cos(np.radians(da.lat)) * dlat / norm_factor
    da_weighted = da_zonal * weight
    da_mean = da_weighted.sum(dim = "lat")
    
    # annual mean (to avoid the yearly cycle)
    if trend == 'yearly':
        da_final_mean = da_mean.groupby('time.year').mean(dim = 'time')
        
    # seasonal mean
    elif trend == 'season_DJF':
        #w_ref = ref.sel(time = slice('2020', '2030')).groupby('time.season').mean(dim = 'time')[0]
        w = da_mean.groupby('time.season')
        da_final_mean = w['DJF'].groupby('time.year').mean(dim = 'time') #- w_ref
        
    elif trend == 'season_MAM':
        #spr_ref = ref.sel(time = slice('2020', '2030')).groupby('time.season').mean(dim = 'time')[2]
        spr = da_mean.groupby('time.season')
        da_final_mean = spr['MAM'].groupby('time.year').mean(dim = 'time') #- spr_ref
        
    elif trend == 'season_JJA':
        #s_ref = ref.sel(time = slice('2020', '2030')).groupby('time.season').mean(dim = 'time')[1]
        s = da_mean.groupby('time.season')
        da_final_mean = s['JJA'].groupby('time.year').mean(dim = 'time') #- s_ref
    
    elif trend == 'season_SON':
        #a_ref = ref.sel(time = slice('2020', '2030')).groupby('time.season').mean(dim = 'time')[3]
        a = da_mean.groupby('time.season')
        da_final_mean = a['SON'].groupby('time.year').mean(dim = 'time') #- a_ref
    
    return da_final_mean


def global_mean1(da, mask = None, ref = None, trend = 'yearly'):

    """
    Takes an xarray DataArray defined over a lat/lon field over time (daily data)
    and calculates the global/regional mean computed with correct weight factors.
    Can also work on an exclusively meridional field (only lat coordinates)

    Parameters:
        da (DataArray): data array containing a field defined over latitudes (lat) and longitudes (lon) and time
        mask (DataArray): data array containing the corresponding mask
        ref (DataArray): data array containing a time series for the reference value
            ONLY for the regional mean (remember it has to be the reference over the region!)
        trend (string): trend to calculate, either yearly or seasonal (specify the season)

    Returns:
        Data array with the global or regional mean of the field. If plotted, it is a time series
    """

    # apply mask
    if not mask is None: # if mask is not None
        da = da * mask.isel(time = 0)
    
    # detect if data has longitudinal component
    if 'lon' in list(da.coords):
        zonal = False
    else:
        zonal = True
    
    # determine the grid size
    # note: this requires an equally spaced coordinate grid
    dlat = np.deg2rad(np.diff(da.lat.data)[0])
    dlon = np.deg2rad(np.diff(da.lon.data)[0])
    
    # determine normalization factors
    if mask is None:
        norm_factor = 2
    else:
        mask_zonal = mask.sum(dim = "lon") * dlon
        mask_int = ((np.cos(np.deg2rad(da.lat))) * dlat * mask_zonal).sum(dim = "lat")
        norm_factor = float(mask_int[0])
    
    # zonal integration
    if not zonal: # if there is longitudinal component
        da_zonal = da.sum(dim = "lon") * dlon
        if mask is None: # apply zonal normalization if no mask is given
            da_zonal *= 1.0 / (2.0 * np.pi)
        else:
            da_zonal = da_zonal
    else: 
        da_zonal = da

    # meridional integration
    weight = 1.0 * np.cos(np.radians(da.lat)) * dlat / norm_factor
    da_weighted = da_zonal * weight
    da_mean = da_weighted.sum(dim = "lat")
    
    return da_mean


def moving_mean(da, X):
       
    """
    Takes an xarray DataArray defined over time and calculates a moving mean over X years and the standard deviation

    Parameters:
        da (DataArray): data array containing a time series
        X (integer): number of years to consider for the window

    Returns:
        Data Array with the calculated mean and its standard deviation
    """
        
    mitjana = da.rolling(year = X, center = 2).mean().compute()
        
    # calculate the standard deviation (measure of the interannual variability)
    std_dev = da.rolling(year = X, center = 2).std().compute()
        
    return mitjana, std_dev



""" AUX FUNCTIONS """

def center_lon(data):
    
    """
    Takes an xarray DataArray and shifts the longitude from 0,360 to -180,180
    
    Parameters:
        data (DataArray): data array defined over a field of lat, lon
        
    Returns: 
        DataArray with the 0 longitude at the center
    """
    
    data_rolled = data.roll(lon = len(data.lon)//2, roll_coords = True)
    
    # update the longitude values to range from -180 to 180
    data_rolled['lon'] = (data_rolled['lon'] + 180) % 360 - 180
    
    return data_rolled


def increase_tantx100(x_initial, x_final):
    
    """
    Takes two xarray DataArray representing the initial and final state 
    and calculates the increase or decrease in percentatge.
    - if x hasn't changed (x_final = x_initial --> increase = 0
    - if x has doubled (x_final = 2 x_intial) --> increase = 100%
    - if x has decrased by half (x_final = 1/2 x_initial) --> increase = -50%

    Parameters:
        x_initial (DataArray): data array containing the initial state
        x_final (DataArray): data array containing the final state

    Returns:
        Data array with the change in percentatge. If plotted, it is a time series
    """

    return ((x_final / x_initial) - 1) * 100


def find_y_final(x, lat_initial, lat_final, lon_initial, lon_final):
    
    """
    Finds the top right latitude as a function of the longitude

    Parameters:
        x: current longitude
        lat_initial: initial latitude of the considered region, starting from -90N
        lat_final: final latitude of the considered region, ending at 90N
        lon_initial: initial longitude of the considered region, starting from the west
        lon_final: final longitude of the considered region, ending at the east

    Returns:
        Data array containing the final latitude depending on the longitude
    """
    
    y_top = ( (lat_final - lat_initial) * (x - lon_initial) / (lon_final - lon_initial)) + lat_initial
    
    # move 360 degrees the longitudes that are on the west
    # (és com moure del darrere al davant les lon>180 per poder calcular bé el triange)
    y_top = y_top.where(y_top.lon < 180, 
                        (lat_final - lat_initial) * (x-360 - lon_initial) / (lon_final - lon_initial) + lat_initial)
    
    return y_top


def maskU(da, y_top1, y_top2, y_top3):
    
    """
    Creates a regional mask for the west of South America

    Parameters:
        da (DataArray): data array containing a field defined over latitudes (lat) and longitudes (lon) and time
        y_top1, y_top2, y_top3 (DataArray): data arrays with the latitude as a function of longitude

    Returns:
        Data array with the same shape as the original one but containing only 1 or 0
    """
    
    mask_u1 = lon_mask(da, -82, -65).where((lon_mask(da, -82, -65).lat >= -19) & 
                                           (lon_mask(da, -82, -65).lat <= y_top1), 0)
    mask_u2 = lon_mask(da, -82, -69).where((lon_mask(da, -82, -69).lat >= -46) & 
                                           (lon_mask(da, -82, -69).lat <= -19), 0)
    mask_u3 = lon_mask(da, -69, -65).where((lon_mask(da, -69, -65).lat >= y_top2) & 
                                           (lon_mask(da, -69, -65).lat <= -19), 0)
    mask_u4 = lon_mask(da, -83, -69).where((lon_mask(da, -83, -69).lat >= -56) & 
                                           (lon_mask(da, -83, -69).lat <= -46), 0)
    mask_u5 = lon_mask(da, -69, -65).where((lon_mask(da, -69, -65).lat >= -56) & 
                                           (lon_mask(da, -69, -65).lat <= y_top3), 0)
    
    mask = mask_u1 + mask_u2 + mask_u3 + mask_u4 + mask_u5
    
    return mask


def maskV(da, y_top1, y_top2):
    
    """
    Creates a regional mask for the south-east of South America

    Parameters:
        da (DataArray): data array containing a field defined over latitudes (lat) and longitudes (lon) and time
        y_top1, y_top2 (DataArray): data arrays with the latitude as a function of longitude

    Returns:
        Data array with the same shape as the original one but containing only 1 or 0
    """
    
    mask_v1 = lon_mask(da, -69, -65).where( (lon_mask(da, -69, -65).lat >= -46) & (lon_mask(da, -69, -65).lat <= y_top1), 0)
    mask_v2 = lon_mask(da, -69, -65).where((lon_mask(da, -69, -65).lat >= y_top2) & (lon_mask(da, -69, -65).lat <= -46), 0)
    mask_v3 = lon_mask(da, -65, -39).where((lon_mask(da, -65, -39).lat >= -56) & (lon_mask(da, -65, -39).lat <= -19), 0)
   
    maskk = mask_v1 + mask_v2 + mask_v3
    mask = np.clip(maskk, 0, 1)
    
    return mask


def maskH(da):
    
    """
    Creates a regional mask for India

    Parameters:
        da (DataArray): data array containing a field defined over latitudes (lat) and longitudes (lon) and time

    Returns:
    Data array with the same shape as the original one but containing only 1 or 0
    """
        
    mask_h1 = lon_mask(da, 61, 94).where(
        (lon_mask(da, 61, 94).lat >= 10) & (lon_mask(da, 61, 94).lat <= 27), 0)
    mask_h2 = lon_mask(da, 61, 100).where(
        (lon_mask(da, 61, 100).lat >= 27) & (lon_mask(da, 61, 100).lat <= 29), 0)
    
    mask = mask_h1 + mask_h2
    
    return mask


def caspian_lake(da):
    
    """
    Creates a mask for the Caspian lake
    """
    
    da_zeros = xr.zeros_like(da)
        
    mask = da_zeros.where( (da.lon > 48) & (da.lon < 53) & (da.lat > 36) & (da.lat < 46), 1 )
    
    return mask



""" REGIONAL STUFF """

def regional_masks(da):
    
    """
    Creates regional masks according to Irvine et al. 2019

    Parameters:
        da (Data Array): data array containing a field defined over latitudes (lat) and longitudes (lon) and time

    Returns:
        Dictionary containing the regional masks
     """
    
    # calculate the final latitude for the non regular regions
    y_topA = find_y_final(da.lon, 43, 72, -12, 39)
    y_topB = find_y_final(da.lon, 43, 72, -12, 39)
    y_topT1 = find_y_final(da.lon, -19, 30, -66, -122)
    y_topT2 = find_y_final(da.lon, 12, 30, -69, -91)
    y_topU1 = find_y_final(da.lon, -21, -5, -65, -81)
    y_topU2 = find_y_final(da.lon, -46, -19, -69, -65)
    y_topU5 = find_y_final(da.lon, -56, -46, -65, -69)
    y_topV1 = find_y_final(da.lon, -46, -19, -69, -65)
    y_topV2 = find_y_final(da.lon, -46, -56, -69, -65)
    
    # create dictionary containing all the regional masks
    masks_zeros = {
        'NEU': lon_mask(da, -12, 39).where((lon_mask(da, -12, 39).lat >= y_topA) & 
                                                    (lon_mask(da, -12, 39).lat <= 72), 0), # a
        'CEU': lon_mask(da, -12, 39).where((lon_mask(da, -12, 39).lat >= 43) & 
                                                     (lon_mask(da, -12, 39).lat <= y_topB), 0), # b
        'MED': create_mask(29, 43, lon_mask(da, -12, 39)), # c
        'SAH': create_mask(16, 29, lon_mask(da, -18, 39)), # d
        'WAF': create_mask(-10, 16, lon_mask(da, -18, 29)), # e
        'EAF': create_mask(-10, 16, lon_mask(da, 30, 51)), # e1
        'SAF': create_mask(-35, -10, lon_mask(da, -6, 51)), # f
        'WAS': create_mask(16, 50, lon_mask(da, 39, 61)) * caspian_lake(da).isel(time = 0), # g
        'SAS': create_mask(10, 29, lon_mask(da, 61, 94)), #maskH(da), # h
        'CAS': create_mask(29, 50, lon_mask(da, 61, 78)), # i
        'TIB': create_mask(29, 50, lon_mask(da, 78, 94)), # j
        'EAS': create_mask(27, 50, lon_mask(da, 94, 146)), # k
        'NAS': create_mask(50, 88, lon_mask(da, 39, 190)), # ñ
        'SEA': create_mask(-10, 27, lon_mask(da, 94, 156)), # l
        'NAU': create_mask(-29, -10, lon_mask(da, 107, 156)), # m
        'SAU': create_mask(-48, -29, lon_mask(da, 107, 190)), # n
        'ALA': create_mask(56, 75, lon_mask(da, -167, -108)), # o
        'CGI': create_mask(50, 88, lon_mask(da, -108, -12)), # p
        'WNA': create_mask(30, 56, lon_mask(da, -135, -108)), # q
        'CNA': create_mask(30, 50, lon_mask(da, -108, -84)), # r
        'ENA': create_mask(25, 50, lon_mask(da, -84, -54)), # s
        'CAM': lon_mask(da, -122, -66).where((lon_mask(da, -122, -66).lat >= y_topT1) & 
                                                         (lon_mask(da, -122, -66).lat <= 30) & 
                                                         (lon_mask(da, -122, -66).lat <= y_topT2), 0), # t
        'WSA': maskU(da, y_topU1, y_topU2, y_topU5), # u
        'SSA': maskV(da, y_topV1, y_topV2), # v
        'AMZ': create_mask(-19, 12, lon_mask(da, -66, -51)), # w
        'NEB': create_mask(-19, 0, lon_mask(da, -51, -34)) } # x   
    
    # apply land mask to the regional masks
    regions = list(masks_zeros.keys())
    var_land_mask = land_mask()
    masks = { region: masks_zeros[region] * var_land_mask.isel(time = 0) for region in regions }
                                                   
    return masks


def regional_plots(da_control, da_sai2020, da_sai2080, ylabel, path, ref_var = None):
    
    """
    Makes regional plots

    Parameters:
        da_control (Data Array): data array containing a time series for the control simulation
        da_sai2020 (Data Array): data array containing a time series for the SAI 2020 simulation
        da_sai2080 (Data Array): data array containing a time series for the SAI 2080 simulation
        ylabel (string): text to write on the y axis
        path: where to save the plots
        ref_var (Data Array): data array containing a time series for the reference value
        
    """

    regions = ['NEU',
               'CEU',
               'MED',
               'SAH',
               'WAF',
               'EAF',
               'SAF',
               'WAS',
               'SAS',
               'CAS',
               'TIB',
               'EAS',
               'NAS',
               'SEA',
               'NAU',
               'SAU',
               'ALA',
               'CGI',
               'WNA',
               'CNA',
               'ENA',
               'CAM',
               'WSA',
               'SSA',
               'AMZ',
               'NEB']

    for i in range(len(regions)):
        print(i)
        fig = plt.figure(figsize = (8, 5))

        # control
        da_control[regions[i]][0].plot(label = 'Control', color = 'blue')
        (da_control[regions[i]][0] - da_control[regions[i]][1]).plot(color = 'blue', alpha = 0.6, linestyle = 'dashed')
        (da_control[regions[i]][0] + da_control[regions[i]][1]).plot(color = 'blue', alpha = 0.6, linestyle = 'dashed')
        plt.fill_between(da_control[regions[i]][0].year, 
                         da_control[regions[i]][0] - da_control[regions[i]][1], 
                         da_control[regions[i]][0] + da_control[regions[i]][1], 
                         alpha = 0.2)

        # sai2020
        da_sai2020[regions[i]][0].plot(label = 'SAI 2020', color = 'orange')
        (da_sai2020[regions[i]][0] - da_sai2020[regions[i]][1]).plot(color = 'orange', alpha = 0.6, linestyle = 'dashed')
        (da_sai2020[regions[i]][0] + da_sai2020[regions[i]][1]).plot(color = 'orange', alpha = 0.6, linestyle = 'dashed')  
        plt.fill_between(da_sai2020[regions[i]][0].year, 
                         da_sai2020[regions[i]][0] - da_sai2020[regions[i]][1], 
                         da_sai2020[regions[i]][0] + da_sai2020[regions[i]][1], 
                         alpha = 0.2)

        # sai2080
        da_sai2080[regions[i]][0].plot(label = 'SAI 2080', color = 'green')
        (da_sai2080[regions[i]][0] - da_sai2080[regions[i]][1]).plot(color = 'green', alpha = 0.6, linestyle = 'dashed')
        (da_sai2080[regions[i]][0] + da_sai2080[regions[i]][1]).plot(color = 'green', alpha = 0.6, linestyle = 'dashed')
        plt.fill_between(da_sai2080[regions[i]][0].year, 
                         da_sai2080[regions[i]][0] - da_sai2080[regions[i]][1], 
                         da_sai2080[regions[i]][0] + da_sai2080[regions[i]][1], 
                         alpha = 0.2)

        if not ref_var is None: # if there is a reference variable:
            plt.axhline(y = ref_var[str(regions[i])], color = 'black', linestyle = 'dashed')
        plt.xticks(fontsize = 11)
        plt.yticks(fontsize = 11)
        plt.xlabel('Year', fontsize = 11)
        plt.ylabel(ylabel, fontsize = 11)
        plt.title(str(regions[i]), fontsize = 12)
        plt.legend()
        plt.grid()
        plt.savefig(path + str(regions[i]) + '.pdf')
        #plt.close()



""" JOINT DISTRIBUTION STUFF """

def last_seasons(da, first_year, X, mask, season = 'winter'):
    
    """
    Takes an xarray DataArray defined over a lat/lon field over time (has to be monthly data!)
    and calculates the X last winters or summers mean

    Parameters:
        da (DataArray): data array containing a field defined over latitudes (lat), longitudes (lon) and time
        first_year (integer): index of the first year from which the average is calculated
        X (integer): number of years to calculate the mean over
        mask (DataArray): data array containing the mask for regional averages 
        season (string): season to calculate, winter or summer

    Returns:
        Data Array containing the last X last seasonal mean
        
    """
    
    assert season in ['winter', 'summer'], f'season must be either winter or summer and not {season}'
    
    if season == 'winter':
        
        # reference winter temperature
        w_ref = (da * mask.isel(time = 0)).sel(time = slice('2020', '2030')).groupby('time.season').mean(
            dim = 'time')[0]   
        w = (da * mask.isel(time = 0)).isel(time = slice(first_year, first_year + (X*12))).groupby('time.season').mean(
            dim = 'time')[0]
        data_final = w - w_ref
        
    elif season == 'summer':
        
        # reference summer temperature
        s_ref = (da * mask.isel(time = 0)).sel(time = slice('2020', '2030')).groupby('time.season').mean(
            dim = 'time')[1] 
        s = (da * mask.isel(time = 0)).isel(time = slice(first_year, first_year + (X*12))).groupby('time.season').mean(
            dim = 'time')[1]
        data_final = s - s_ref
        
    return data_final

def joint_plots(x_data, y_data, x_labels, y_labels, ranges, area, titles, big_title, name, path):
    
    """
    Makes the plots for the joint distribution

    Parameters:
        x_data (list): list containing the data to put in the x axis. Every element must be a data array
        y_data (list): same as x_data but now for the y axis
        x_labels (list): list containing the name of the x axis. Every element must be a string
        y_labels (list): same as x_labels but now for the y axis
        ranges (list): list containing the ranges for every subplot. Every element must be a list
        area (Data Array): data arra ycontaining the area variable
        titles (list): list containing the titles for every subplot. Every element must be a string
        big_title (string): text to write as the main title
        path: where to save the plots

    """
    
    fig = plt.figure(figsize = (20, 6))
    g = gs.GridSpec(nrows=1, ncols=4, width_ratios=[1,1,1, 0.05])
    a = [fig.add_subplot(g[n,m]) for n in range(1) for m in range(3)]
    cax = fig.add_subplot(g[:,3])

    for i, ax in enumerate(a):
        im = ax.hist2d(x_data[i], y_data[i], bins = 80, range = ranges[i], cmin = 1, weights = area.data.flatten()/1400)
        ax.set_title(titles[i], fontsize = 16)
        ax.set_xlabel(x_labels[i], fontsize = 16)
        ax.xaxis.set_tick_params(labelsize = 16)
        ax.set_ylabel(y_labels[i], fontsize = 16)
        ax.yaxis.set_tick_params(labelsize = 16)

        ax.set_aspect('equal')
        
        ax.plot([-60, 60], [60, -60], color = 'gray')
        ax.plot([-60, 60], [-60, 60], color = 'gray')

        # plot lines at x = 0 and y = 0
        ax.axvline(x = 0, color = 'black', linestyle = 'dashed')
        ax.axhline(y = 0, color = 'black', linestyle = 'dashed')
    
    # set colorbar
    cb = fig.colorbar(im[-1], cax = cax, orientation = 'vertical')
    cb.ax.tick_params(labelsize = '14')
    
    plt.suptitle(big_title, fontsize = 18)
    fig.tight_layout()
    plt.savefig(path + name + '.png')
    
    
def regional_plots_joint(da_control, range1, da_sai2020, range2, da_sai2080, range3, regions, area, path):

    """
    Makes regional plots

    Parameters:
        da_control (Data Array): data array containing a time series for the control simulation
        range1 (list): range of the first plot (eg: [[10, 10], [10, 10]])
        da_sai2020 (Data Array): data array containing a time series for the SAI 2020 simulation
        range2 (list): range of the second plot
        da_sai2080 (Data Array): data array containing a time series for the SAI 2080 simulation
        range3 (list): range of the third plot
        regions (list): list of regions to plot
        area (): area variable
        path: where to save the plots

    """

    for i in range(len(regions)):
        print(i)

        # all year
        plt.figure(figsize = (20, 6))

        # sai20 vs control
        plt.subplot(1, 3, 1)
        plt.hist2d(da_control[regions[i]].data.flatten(), da_sai2020[regions[i]].data.flatten(), bins = [85, 100], 
                   range = range1, cmin = 1, weights = area.data.flatten()/1400)

        # plot diagonal lines
        plt.plot([-60, 60], [60, -60], color = 'gray')
        plt.plot([-60, 60], [-60, 60], color = 'gray')
        plt.gca().set_aspect('equal')

        # plot lines at x = 0 and y = 0
        plt.axvline(x = 0, color = 'black', linestyle = 'dashed')
        plt.axhline(y = 0, color = 'black', linestyle = 'dashed')

        plt.xlabel('Control', fontsize = 13)
        plt.ylabel('SAI 2020', fontsize = 13)
        plt.title('SAI 2020 vs Control', fontsize = 14)

        # sai80 vs control
        plt.subplot(1, 3, 2)
        plt.hist2d(da_control[regions[i]].data.flatten(), da_sai2080[regions[i]].data.flatten(), bins = [85, 100],
                   range = range2, cmin = 1, weights = area.data.flatten()/1400)

        # plot diagonal lines
        plt.plot([-60, 60], [60, -60], color = 'gray')
        plt.plot([-60, 60], [-60, 60], color = 'gray')
        plt.gca().set_aspect('equal')

        # plot lines at x = 0 and y = 0
        plt.axvline(x = 0, color = 'black', linestyle = 'dashed')
        plt.axhline(y = 0, color = 'black', linestyle = 'dashed')

        plt.xlabel('Control', fontsize = 13)
        plt.ylabel('SAI 2080', fontsize = 13)
        plt.title('SAI2080 vs Control', fontsize = 14)

        # sai80 vs sai20
        plt.subplot(1, 3, 3)
        plt.hist2d(da_sai2020[regions[i]].data.flatten(), da_sai2080[regions[i]].data.flatten(), bins = [85, 100],
                   range = range3, cmin = 1, weights = area.data.flatten()/1400)

        # plot diagonal lines
        plt.plot([-60, 60], [60, -60], color = 'gray')
        plt.plot([-60, 60], [-60, 60], color = 'gray')
        #plt.gca().set_aspect('equal')

        # plot lines at x = 0 and y = 0
        plt.axvline(x = 0, color = 'black', linestyle = 'dashed')
        plt.axhline(y = 0, color = 'black', linestyle = 'dashed')

        plt.xlabel('SAI 2020', fontsize = 13)
        plt.ylabel('SAI 2080', fontsize = 13)
        plt.title('SAI 2080 vs SAI 2020', fontsize = 14)

        plt.colorbar()        
        plt.xticks(fontsize = 11)
        plt.yticks(fontsize = 11)
        plt.suptitle(str(regions[i]), fontsize = 12)
        plt.savefig(path + str(regions[i]) + '.pdf')



