from pathlib import Path
import numpy as np
import gridit as gi
import pandas as pd

from setup import *
from utils import *



def dis_setup(fn_out):
    grid = gi.Grid.from_vector(GR_SHP, RES)
    aw_stream = grid.array_from_vector(DRAINS)
    # grid_gpd = grid.cell_geodataframe()
    # # grid_gpd.to_file(Path(SPATIAL_DIR, f'{MODEL_NAME}_grid.shp'), driver='ESRI Shapefile')

    # top & bottom
    top = grid.array_from_raster(TOP)
    bottom = grid.array_from_raster(BOTTOM)  # get the bottom elevation
    # top = np.ones_like(bottom) * 10  # flat top at 10mRL
    
    # pp_grid = grid.array_from_vector(GR_SHP).data
    
    # # open domain
    arr = ~grid.array_from_vector(GR_SHP).mask # want the binary version of the domain
    arr = np.where(top.data > 16, 0, arr) # remove areas where elev > 18mRL (mainly at boundaries)
    # arr = np.where((grid == 1) & (arr == 1), 2, arr)  # set pilot point areas to active

    nrow = grid.shape[0]
    ncol = grid.shape[1]
    delr = np.ones(ncol) * RES
    delc = np.ones(nrow) * RES

    # dis
    shallow_bottom = -5 # top of clay /confining layer
    top = top.data  # use the top elevation for the first layer, only where idomain is active
    bottom = np.where(bottom.data < shallow_bottom, shallow_bottom, bottom.data)  # set bottom elevation to top elevation where it is higher
    model_thickness = top.data - bottom  # calculate the thickness of the layers

    if np.any(model_thickness <= 0):
        bottom = np.where(model_thickness <= 0, top - 1, bottom)  # ensure bottom is always below top
        model_thickness = top.data - bottom  # calculate the thickness of the layers

    arr = np.where((model_thickness > 2) & (arr == 1), 1, 0)
    idomain = np.stack([arr] * NLAY, axis=0)
    
    model_thickness = model_thickness * idomain[0]
    
    botm = [bottom]
    # min_b = 1 # min thickness

    # thicknesses = []
    # botm = []
    # b0 = top.copy()
    # for i in range(NLAY):
    #     if i in [0, 1, 2]:  # upper 3 layers with equal thickness
    #         b = 1
    #         ibotm = b0 - b  # calculate the bottom elevation for each layer
    #     else:   
    #         b = np.where(b0 - model_thickness < min_b, min_b, b0 - model_thickness)  # set the thickness of the confining layer
    #         ibotm = b0 - b  # calculate the bottom elevation for each layer
    #     thicknesses.append(b)
    #     botm.append(ibotm)
    #     b0 = ibotm.copy()  # update the bottom elevation for the next layer

    # botm = np.array(botm)

    # save
    fn_out['dis'] = {}
    top_fn = Path(MODEL_DIR, f'{MODEL_NAME}.dis_top.txt').as_posix()
    fn_out['dis']['top'] = tomf6input(top_fn)
    np.savetxt(top_fn, top)  # save top elevation

    fn_out['dis']['btm'] = []
    fn_out['dis']['idomain'] = []
    for i in range(NLAY):
        ilay = i + 1
        idomain_fn = Path(MODEL_DIR, f'{MODEL_NAME}.dis_idomain_layer{ilay}.txt').as_posix()
        btm_fn = Path(MODEL_DIR, f'{MODEL_NAME}.dis_botm_layer{ilay}.txt').as_posix()
        fn_out['dis']['btm'].append(tomf6input(btm_fn))
        fn_out['dis']['idomain'].append(tomf6input(idomain_fn))
        np.savetxt(btm_fn, botm[i])
        # save idomain as int
        idomain = idomain.astype(int)
        np.savetxt(idomain_fn, idomain[i], fmt='%d')

    return grid, idomain, top, nrow, ncol, delr, delc, fn_out, model_thickness


def ghb_aw_setup(grid, idomain, start, end, fn_out, save=True):

    top_min = grid.array_from_raster(TOP, resampling='min')
    ghb_aw_arr = ~grid.array_from_vector(AWANUI_BOUNDARY).mask
    pos_idomain = np.where(idomain[0] > 0, 1, 0)
    in_idomarr, idom_i = get_interior_indices(pos_idomain, layer=0)  # idomainget interior indices of the mbr area
    ghb_aw_mask = np.logical_and(ghb_aw_arr.data, in_idomarr)  # mbr area that is active in the model domain
    
    ghb_aw = ghb_aw_mask * top_min.data
    ghb_aw = np.where(ghb_aw > 8.5, 0, ghb_aw) # remove areas above 8mRL (mainly at boundaries)
    riv_stage_arr = np.where(ghb_aw > 0, AW_WL, 0)  # set max awanui water level to 6.62mRL (surveyed high flow level)
    
    # riv_mask = ~grid.array_from_vector(DRAINS).mask * idomain[0]
    # zones = grid.array_from_vector(DRAIN_ZONES, attribute='elev_ss_0')
    # riv_stage_zones = ~zones.mask
    # riv_rbot = grid.array_from_raster(TOP, resampling='min').data
    # # adjust riv elevations absed on survey data
    # riv_stage_arr = np.where(riv_stage_zones, zones.data, riv_rbot + AWANUI_water_offset) * riv_mask
    
    # save to spatial
    if save:
        riv_stage_fn = Path(SPATIAL_DIR, f'{MODEL_NAME}_riv_stage.dat')
        np.savetxt(riv_stage_fn, riv_stage_arr)
        
    # to dataframe
    riv_kper0 = extract_value_with_indices(
        riv_stage_arr, layer=0, val_col='head', mask_value=0
        )
    
    # riv_kper0['cond'] = GHB_SW['initial_cond_aw']
    # conf_kper0['head'] = GHB_CONF['initial_head']
    cond = grid.array_from_raster(SHALLOW_K).data
    cond = (cond * 2 * 25) / 2 # calculating conductance from K (m/d) and cell length & 2m width (m2) & thickness (m)
    cond_df = extract_value_with_indices(
            cond, layer=0, val_col='cond', mask_value=0
            )
    riv_kper0 = pd.merge(riv_kper0, cond_df, left_on='index', right_on='index')

    # time series for riv elevations ###########################################################
    # riv absolute values during recession
    riv_ts_values = get_awanui_timeseries(start, end)

    riv_kper1 = riv_kper0[['index']].copy()
    riv_kper1['head'] = riv_kper1['index'].apply(lambda x: f's_1_{x[1]+1}_{x[2]+1}')
    riv_kper1['cond'] = riv_kper0['cond']

    # kper2: PAST ############################################################
    ghb_aw_past = ~grid.array_from_vector(AWANUI_PAST).mask
    pos_idomain = np.where(idomain[0] > 0, 1, 0)
    in_idomarr, idom_i = get_interior_indices(pos_idomain, layer=0)  # idomainget interior indices of the mbr area
    ghb_aw_past = np.logical_and(ghb_aw_past.data, in_idomarr)  # mbr area that is active in the model domain
    
    
    riv_mask = ~grid.array_from_vector(DRAINS_PAST).mask + ghb_aw_past
    spring_mask = grid.array_from_vector(SPRING_DRAIN).mask

    past_h = np.where(
            riv_mask * spring_mask * idomain[0],
            AW_WL + PAST_OFFSET,
            0)
    # past_h = np.where(past_h > 8.14, 8.14, past_h)
    
    riv_kper2 = extract_value_with_indices(
        past_h, layer=0, val_col='head', mask_value=0
        )
    riv_kper2 = pd.merge(riv_kper2, cond_df, left_on='index', right_on='index')

    ts_0 = river_stage(riv_kper0, [0, 0], start_t = 0)
    ts_1 = river_stage(riv_kper0, riv_ts_values, start_t = ts_0.index[-1]+1)
    ts_2 = river_stage(riv_kper2, [0], start_t = ts_1.index[-1]+1)
    ts_3 = river_stage(riv_kper2, riv_ts_values*0.5, start_t = ts_2.index[-1]+1)

    riv_kper1 = riv_kper0.copy()
    riv_kper1['head'] = ts_1.columns
    
    riv_kper3 = riv_kper2.copy()
    riv_kper3['head'] = ts_3.columns

    # convert to dataframe
    rivstage_ts = pd.concat([ts_0, ts_1, ts_2, ts_3]).fillna(0)

    # save
    fn_out['ghb_aw'] = {}
    for i, dat in enumerate([riv_kper0, riv_kper1, riv_kper2, riv_kper3]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghbaw_stress_period_data_{i}.txt')
        fn_out['ghb_aw'][i] = tomf6input(fn, list=True)
        savedf2txt(dat, filename=fn, col_order=['head', 'cond'])
    # save ts file
    riv_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_aw_heads.csv')
    rivstage_ts.index.name = '#time'
    rivstage_ts.to_csv(riv_ts_fn, header=False)
    #save headers to a csv
    rivstage_ts.columns.to_series().to_csv(Path(MODEL_DIR, f'{MODEL_NAME}.ghb_aw_head_names.csv'), index=False, header=False)
    fn_out['ghb_aw_ts'] = tomf6tsinput(riv_ts_fn, rivstage_ts)

    return fn_out, riv_ts_values, 9.21, riv_stage_arr, past_h


def get_awanui_timeseries(start, end):
    riv_ts = pd.read_csv(AWANUI_TS)
    riv_ts['DateTime'] = riv_ts['Awanui Stream at Flume (x)'].apply(pd.to_datetime, format = "%Y-%m-%d %H:%M:%S", dayfirst=True)
    riv_ts.set_index('DateTime', inplace=True)
    riv_ts['level_mRL'] = (riv_ts['Awanui Stream at Flume (y)'] * 0.001) - 10

    # get specific dates
    riv_ts = riv_ts[(riv_ts.index >= start) & (riv_ts.index < end)]
    riv_ts['sm_level_mRL'] = riv_ts['level_mRL'].rolling(window=21, center=True, min_periods=1).mean()
    riv_ts['sm_abs_val'] = riv_ts['sm_level_mRL'] - riv_ts['sm_level_mRL'].iloc[0]
    
    # save raw ts
    riv_ts_raw_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ts_awanui_raw.csv')
    riv_ts.index.name = 'DateTime'
    riv_ts.to_csv(riv_ts_raw_fn)

    # riv absolute values during recession
    return riv_ts['sm_abs_val'][1:].values

def get_poukawa_timeseries(start, end):
    pw_ts = pd.read_csv(POUKAWA_TS)
    pw_ts['DateTime'] = pw_ts['Poukawa Stream at Douglas Road (x)'].apply(pd.to_datetime, format = "%Y-%m-%d %H:%M:%S", dayfirst=True)
    pw_ts['level_mRL'] = (pw_ts['Poukawa Stream at Douglas Road (y)'] * 0.001) - 10
    pw_ts.set_index('DateTime', inplace=True)

    # get specific dates
    pw_ts = pw_ts[(pw_ts.index >= start) & (pw_ts.index < end)]
    pw_ts['sm_val'] = pw_ts['level_mRL'].rolling(window=5, center=True, min_periods=1).mean()
    pw_ts['sm_abs_val'] = pw_ts['sm_val'] - pw_ts['sm_val'].iloc[0]

    # save raw ts
    pw_ts_raw_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ts_poukawa_raw.csv')
    pw_ts.index.name = 'DateTime'
    pw_ts.to_csv(pw_ts_raw_fn)
    
    # riv absolute values during recession
    return pw_ts['sm_abs_val'][1:].values


def river_stage(ghb_cells, riv_ts_values, start_t = 1, stper=1):
    """
    ghb_cells: index and initial head values
    riv_ts_values: time series values to add to initial head
    start_t: period start time
    """
    rivstage_dat = {}
    for row in range(len(ghb_cells)):
        index = ghb_cells.at[row, 'index']
        col = f's_1_{index[1]+1}_{index[2]+1}'
        initial_elev = ghb_cells.at[row, 'head']
        kper_list = (initial_elev + riv_ts_values).tolist()
        rivstage_dat[col] = kper_list
    riv_df = pd.DataFrame(
        rivstage_dat, 
        index=np.arange(start_t, len(riv_ts_values) + start_t, stper)
    )

    return riv_df

def conf_stage(ghb_cells, riv_ts_values, start_t = 1):
    """
    ghb_cells: index and initial head values
    riv_ts_values: time series values to add to initial head
    start_t: period start time
    """
    rivstage_dat = {}
    for row in range(len(ghb_cells)):
        index = ghb_cells.at[row, 'index']
        col = f'cf_1_{index[1]+1}_{index[2]+1}'
        initial_elev = ghb_cells.at[row, 'head']
        kper_list = (initial_elev + riv_ts_values).tolist()
        rivstage_dat[col] = kper_list
    riv_df = pd.DataFrame(
        rivstage_dat, 
        index=np.arange(start_t, len(riv_ts_values) + start_t)
    )

    return riv_df


def ghb_spring_setup(grid, idomain, start, end, aw_ts, wetland_WL, aw_present_arr, aw_past_arr, fn_out, save=True):
    spring_mask = ~grid.array_from_vector(SPRING_DRAIN).mask
    spring_mask = np.where((idomain[0] > 0) & (aw_present_arr == 0), spring_mask, 0)
    spring_stage_arr = np.where(spring_mask, SP_WL, 0)
    # spring_stage_arr[spring_stage_arr > 7.59] = 7.59 # correct to no higher than surveyed max
    # spring_stage_arr = np.where(aw_present_arr>0, 0, spring_stage_arr)  # remove areas also covered by awanui ghb

    # past spring_stage
    past_WL = np.max(aw_past_arr)
    past_spring_stage = np.where(spring_mask, past_WL, 0)
    # past_spring_stage = np.where(past_spring_stage > wetland_WL, wetland_WL, past_spring_stage)
    past_spring_stage = np.where(aw_past_arr>0, 0, past_spring_stage)  # remove areas also covered by awanui ghb

    # save to spatial
    if save:
        spring_stage_arr_fn = Path(SPATIAL_DIR, f'{MODEL_NAME}_spring_stage.dat')
        np.savetxt(spring_stage_arr_fn, spring_stage_arr)
    
    # to dataframe
    spring_kper0 = extract_value_with_indices(
        spring_stage_arr, layer=0, val_col='head', mask_value=0
        )
    #conductance
    cond = grid.array_from_raster(SHALLOW_K).data
    cond = (cond * 2 * 25) / 2 # calculating conductance from K (m/d) and cell length & 2m width (m2) & thickness (m)
    cond_df = extract_value_with_indices(
            cond, layer=0, val_col='cond', mask_value=0
            )
    spring_kper0 = pd.merge(spring_kper0, cond_df, left_on='index', right_on='index')

    # spring_kper0['cond'] = GHB_SW['initial_cond_spr']
    spring_kper2 = extract_value_with_indices(
        past_spring_stage, layer=0, val_col='head', mask_value=0
        )
    spring_kper2 = pd.merge(spring_kper2, cond_df, left_on='index', right_on='index')


    # time series for riv elevations ###########################################################
    spring_ts = pd.read_csv(SPRING_TS)
    spring_ts['DateTime'] = spring_ts['DateTime'].apply(pd.to_datetime, format = "%Y-%m-%d")
    spring_ts.set_index('DateTime', inplace=True)
    spring_ts['level_mRL'] = spring_ts['LEVEL']
    spring_ts['sm_level_mRL'] = spring_ts['level_mRL'].rolling(window=21, center=True, min_periods=1).mean()

    # get specific dates
    spring_ts = spring_ts[(spring_ts.index >= start) & (spring_ts.index < end)]
    spring_ts['sm_val'] = spring_ts['level_mRL'].rolling(window=5, center=True, min_periods=1).mean()
    spring_ts['sm_abs_val'] = spring_ts['sm_val'] - spring_ts['sm_val'].iloc[0]
    
    # save raw ts
    spring_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ts_spring_raw.csv')
    spring_ts.reset_index(inplace=True)
    spring_ts.index = np.arange(1, len(spring_ts)+1)
    spring_ts.to_csv(spring_ts_fn)

    # riv absolute values during recession
    spring_ts_values = spring_ts['sm_abs_val'][1:].values
    # correct so values do not exceed 0
    spring_ts_values = np.where(spring_ts_values > 0, 0, spring_ts_values)

    ts_0 = river_stage(spring_kper0, [0, 0], start_t = 0)
    ts_1 = river_stage(spring_kper0, spring_ts_values, start_t = ts_0.index[-1]+1)
    ts_2 = river_stage(spring_kper2, [0], start_t = ts_1.index[-1]+1)
    ts_3 = river_stage(spring_kper2, aw_ts*0.5, start_t = ts_2.index[-1]+1)

    spring_kper1 = spring_kper0.copy()
    spring_kper1['head'] = ts_1.columns

    spring_kper3 = spring_kper2.copy()
    spring_kper3['head'] = ts_3.columns

    # convert to dataframe
    spring_h_ts = pd.concat([ts_0, ts_1, ts_2, ts_3]).fillna(0)

    # save
    fn_out['ghb_sp'] = {}
    for i, dat in enumerate([spring_kper0, spring_kper1, spring_kper2, spring_kper3]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghbspr_stress_period_data_{i}.txt')
        fn_out['ghb_sp'][i] = tomf6input(fn, list=True)
        savedf2txt(dat, filename=fn, col_order=['head', 'cond'])
    # save ts file
    spring_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_spring_heads.csv')
    spring_h_ts.index.name = '#time'
    # spring_h_ts.to_csv(spring_ts_fn)
    spring_h_ts.to_csv(spring_ts_fn, header=False)
    spring_h_ts.columns.to_series().to_csv(Path(MODEL_DIR, f'{MODEL_NAME}.ghb_spring_head_names.csv'), index=False, header=False)
    fn_out['ghb_sp_ts'] = tomf6tsinput(spring_ts_fn, spring_h_ts)

    return fn_out, spring_stage_arr, past_spring_stage

def ghb_pw_setup(grid, idomain, start, end, aw_past_arr, fn_out, save=True):
    # top_min = grid.array_from_raster(TOP, resampling='min')
    ghb_pw_arr = ~grid.array_from_vector(POUKAWA_BOUNDARY).mask
    pos_idomain = np.where(idomain[0] > 0, 1, 0)

    in_idomarr, idom_i = get_interior_indices(pos_idomain, layer=0)  # idomainget interior indices of the mbr area
    ghb_pw_active = np.logical_and(ghb_pw_arr.data, in_idomarr) * idomain[0] # mbr area that is active in the model domain
    ghb_pw_active = np.where(ghb_pw_active, PW_WL, 0)

    # save to spatial
    if save:
        ghb_pw_active_fn = Path(SPATIAL_DIR, f'{MODEL_NAME}_poukawa_stage.dat')
        np.savetxt(ghb_pw_active_fn, ghb_pw_active)

    ghb_pw_indices = get_indices(ghb_pw_active, layer=0, value=True)

    ghb_pw_kper0 = pd.DataFrame({'index': [i[0] for i in ghb_pw_indices]})
    ghb_pw_kper0['head'] = [i[1] for i in ghb_pw_indices]  # head for Poukawa boundary
    # ghb_pw_kper0['cond'] = GHB_SW['initial_cond_pw']
    #conductance
    cond = grid.array_from_raster(SHALLOW_K).data
    cond = (cond * 2 * 25) / 2 # calculating conductance from K (m/d) and cell length & 2m width (m2) & thickness (m)
    cond_df = extract_value_with_indices(
            cond, layer=0, val_col='cond', mask_value=0
            )
    ghb_pw_kper0 = pd.merge(ghb_pw_kper0, cond_df, left_on='index', right_on='index')

    # TS ##########################################################################
    pw_ts_values = get_poukawa_timeseries(start, end)
    aw_ts_values = get_awanui_timeseries(start, end)

    # kper3 - PAST SS ########################################################
    past_wl = np.max(aw_past_arr)

    ghb_pw_arr = ~grid.array_from_vector(PW_PAST).mask
    active_domain = np.where(idomain[0] > 0, 1, 0)
    in_idomarr, idom_i = get_interior_indices(active_domain, layer=0)  # get edges of idomain
    ghb_pw_past = np.logical_and(ghb_pw_arr.data, in_idomarr)
    ghb_pw_past = ghb_pw_past * past_wl
    
    ghb_pw_indices = get_indices(ghb_pw_past, layer=0, value=True)
    ghb_pw_kper2 = pd.DataFrame({'index': [i[0] for i in ghb_pw_indices]})
    ghb_pw_kper2['head'] = [i[1] for i in ghb_pw_indices]  # head for Poukawa boundary
    ghb_pw_kper2 = pd.merge(ghb_pw_kper2, cond_df, left_on='index', right_on='index')

    ts_0 = river_stage(ghb_pw_kper0, [0, 0], start_t = 0)
    ts_1 = river_stage(ghb_pw_kper0, pw_ts_values, start_t = ts_0.index[-1]+1)
    ts_2 = river_stage(ghb_pw_kper2, [0], start_t = ts_1.index[-1]+1)
    ts_3 = river_stage(ghb_pw_kper2, aw_ts_values*0.5, start_t = ts_2.index[-1]+1) #all the same

    ghb_pw_kper1 = ghb_pw_kper0.copy()
    ghb_pw_kper1['head'] = ts_1.columns
    
    ghb_pw_kper3 = ghb_pw_kper2.copy()
    ghb_pw_kper3['head'] = ts_3.columns

    # convert to dataframe
    pwelev_ts = pd.concat([ts_0, ts_1, ts_2, ts_3]).fillna(0)

    #save
    fn_out['ghb_pw'] = {}
    for i, dat in enumerate([ghb_pw_kper0, ghb_pw_kper1, ghb_pw_kper2, ghb_pw_kper3]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_pw_stress_period_data_{i}.txt')
        fn_out['ghb_pw'][i] = tomf6input(fn ,list=True)
        savedf2txt(dat, filename=fn, col_order=['head', 'cond'])
    
    # save ts file
    pwelev_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_pw_heads.csv')
    pwelev_ts.index.name = '#time'
    pwelev_ts.to_csv(pwelev_ts_fn, header=False)
    pwelev_ts.columns.to_series().to_csv(Path(MODEL_DIR, f'{MODEL_NAME}.ghb_pw_head_names.csv'), index=False, header=False)
    fn_out['ghb_pw_ts'] = tomf6tsinput(pwelev_ts_fn, pwelev_ts)
    return fn_out, ghb_pw_kper0, ghb_pw_active, ghb_pw_past


def get_confined_timeseries(start, end):
    conf_ts = pd.read_csv(CONF_TS)
    conf_ts['DateTime'] = conf_ts['DateTime'].apply(pd.to_datetime, format = "%d/%m/%Y %H:%M")
    conf_ts.set_index('DateTime', inplace=True)
    conf_ts['level_mRL'] = conf_ts['LEVEL'] + 10.755
    #resample to daily
    conf_ts.loc[(conf_ts.index >= end), :] = 11.3  # set end value to 11.3

    # fill missing with interpolation
    conf_ts = conf_ts[['level_mRL']].resample('D').mean()
    conf_ts['level_mRL'] = conf_ts['level_mRL'].interpolate(method='linear')
    conf_ts = conf_ts[(conf_ts.index >= start) & (conf_ts.index < end)]
    
    conf_ts['sm_abs_val'] = conf_ts['level_mRL'] - conf_ts['level_mRL'][0]
    
    # save raw ts
    conf_ts_raw_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ts_confined_raw.csv')
    conf_ts.index.name = 'time'
    conf_ts.reset_index(inplace=True)
    conf_ts.index = np.arange(1, len(conf_ts)+1)
    conf_ts.to_csv(conf_ts_raw_fn)

    # ts_0 = pd.Series(
    #     [0], 
    #     index=[1]
    # )
    # ts_0.name = 'heads_ts'
    # # convert to dataframe
    ts_1 = conf_ts['sm_abs_val'][1:].reset_index(drop=True)
    # ts_1.name = 'heads_ts'
    # ts_1.index = np.arange(2, len(ts_1) + 2)

    # past (normalize between 13.95 and 12.95)
    # ts_2 = ts_0.copy()
    # ts_2.index = np.arange(ts_1.index[-1]+1, ts_1.index[-1]+2)
    ts_3 = ts_1.copy()
    ts_3 = (ts_3 - 0) / (ts_3.min() - 0) * -1
    # ts_3.index = np.arange(ts_2.index[-1]+1, ts_2.index[-1]+len(ts_3)+1)

    # conf_ts = pd.concat([ts_0, ts_1, ts_2, ts_3]).fillna(0)

    return ts_1.values, ts_3.values



def ghb_conf_setup(grid, idomain, start, end, model_thickness, fn_out):
    conf_arr = ~grid.array_from_vector(CONF_AREA_ACTIVE).mask * idomain[0]
    conf_arr = np.where(model_thickness < 15, 0, conf_arr)

    # elev_max = grid.array_from_raster(TOP).data
    # conf_h = elev_max + GHB_CONF['initial_head_offset']
    # conf_h = np.where(conf_h < GHB_CONF['initial_head_min'], GHB_CONF['initial_head_min'], conf_h)
    conf_arr = conf_arr * GHB_CONF['initial_head_min']
    

    conf_kper0 = extract_value_with_indices(
        conf_arr, layer=0, val_col='head', mask_value=0
        )
    # conf_kper0['head'] = GHB_CONF['initial_head']
    conf_cond = grid.array_from_raster(CONF_K).data
    conf_cond = (conf_cond * RES * RES) / CONF_THICKNESS # calculating conductance from K (m/d) and cell area (m2) & thickness (m)
    conf_cond = np.where(conf_arr > 0, conf_cond, 0)
    conf_cond_df = extract_value_with_indices(
            conf_cond, layer=0, val_col='cond', mask_value=0
            )

    conf_kper0 = pd.merge(conf_kper0, conf_cond_df, left_on='index', right_on='index')
    conf_kper2 = conf_kper0.copy()

    conf_ts1, conf_ts3 = get_confined_timeseries(start, end)

    #NEW
    ts_0 = conf_stage(conf_kper0, [0, 0], start_t = 0)
    ts_1 = conf_stage(conf_kper0, conf_ts1, start_t = ts_0.index[-1]+1)
    ts_2 = conf_stage(conf_kper2, [0], start_t = ts_1.index[-1]+1)
    ts_3 = conf_stage(conf_kper2, conf_ts3, start_t = ts_2.index[-1]+1)

    conf_kper1 = conf_kper0.copy()
    conf_kper1['head'] = ts_1.columns
    
    conf_kper3 = conf_kper2.copy()
    conf_kper3['head'] = ts_3.columns

    # convert to dataframe
    conf_ts = pd.concat([ts_0, ts_1, ts_2, ts_3]).fillna(0)
    
    #save
    fn_out['ghb_conf'] = {}
    for i, dat in enumerate([conf_kper0, conf_kper1, conf_kper2, conf_kper3]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_conf_stress_period_data_{i}.txt')
        fn_out['ghb_conf'][i] = tomf6input(fn ,list=True)
        savedf2txt(dat, filename=fn, col_order=['head', 'cond'])
    
    # save ts file
    conf_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_conf_heads.csv')
    conf_ts.index.name = '#time'
    conf_ts.to_csv(conf_ts_fn, header=False)
    conf_ts.columns.to_series().to_csv(Path(MODEL_DIR, f'{MODEL_NAME}.ghb_conf_head_names.csv'), index=False, header=False)
    fn_out['ghb_conf_ts'] = tomf6tsinput(conf_ts_fn, conf_ts)
    return fn_out, conf_arr

def wel_mbr_setup(grid, idomain, fn_out, ghb_pw_kper0):
    mbr_arr = ~grid.array_from_vector(MBR).mask
    mbr_indices = []
    for i in range(NLAY):
        in_idomarr, idom_i = get_interior_indices(idomain[i], layer=i)  # idomainget interior indices of the mbr area
        mbr_active = np.logical_and(mbr_arr.data, in_idomarr)  # mbr area that is active in the model domain
        mbr_indices.extend(get_indices(mbr_active, layer=i))
    
    mbr_df = pd.DataFrame({'index': mbr_indices})
    mbr_df['q'] = 138 / len(mbr_df) # m3/d
    mbr_df = mbr_df[~mbr_df['index'].isin(ghb_pw_kper0['index'])]  # remove chd indices from mbr
    
    # save
    fn_out['wel_mbr'] = {}
    for i in range(NPER):
        mbr_fn = Path(MODEL_DIR, f'{MODEL_NAME}.wel_mbr_stress_period_data_{i}.txt')
        fn_out['wel_mbr'][0] = tomf6input(mbr_fn, list=True)
        savedf2txt(mbr_df, filename=mbr_fn, col_order=['q'])
    return fn_out, mbr_df

def wel_inout_setup(grid, idomain, mbr_df, fn_out):
    influx_arr = ~grid.array_from_vector(INFLUX_BOUNDARY).mask
    influx_indices = []
    for i in range(NLAY): # remove influx from bottom 2 layers
        in_idomarr, idom_i = get_interior_indices(idomain[i], layer=i)  # idomainget interior indices of the mbr area
        influx_active = np.logical_and(influx_arr.data, in_idomarr)  # mbr area that is active in the model domain
        influx_indices.extend(get_indices(influx_active, layer=i))
    influx_df = pd.DataFrame({'index': influx_indices})

    # outflux boundaries
    outflux_arr = ~grid.array_from_vector(OUTFLUX_BOUNDARY).mask
    outflux_indices = []
    for i in range(NLAY):
        in_idomarr, idom_i = get_interior_indices(idomain[i], layer=i)  # idomainget interior indices of the mbr area
        outflux_active = np.logical_and(outflux_arr.data, in_idomarr)  # mbr area that is active in the model domain
        outflux_indices.extend(get_indices(outflux_active, layer=i))
    outflux_df = pd.DataFrame({'index': outflux_indices})

    # influx/outfluc
    influx_df['q'] = 1  # m3/d
    outflux_df['q'] = -1  # m3/d

    # remove where mbr is active
    influx_df = influx_df[~influx_df['index'].isin(mbr_df['index'])]
    outflux_df = outflux_df[~outflux_df['index'].isin(mbr_df['index'])]

    # save
    fn_out['wel_influx'] = {}
    fn_out['wel_outflux'] = {}
    for i in range(NPER):
        influx_fn = Path(MODEL_DIR, f'{MODEL_NAME}.wel_influx_stress_period_data_{i}.txt')
        outflux_fn = Path(MODEL_DIR, f'{MODEL_NAME}.wel_outflux_stress_period_data_{i}.txt')
        fn_out['wel_influx'][0] = tomf6input(influx_fn, list=True)
        fn_out['wel_outflux'][0] = tomf6input(outflux_fn, list=True)
        savedf2txt(influx_df, filename=influx_fn, col_order=['q'])
        savedf2txt(outflux_df, filename=outflux_fn, col_order=['q'])
    
    return fn_out

def npf_setup(grid, idomain, fn_out):
    all_k = []
    for i in range(NLAY):
        k = grid.array_from_raster(SHALLOW_K).data
        # k = np.ones_like(idomain[0]) * KH_PRIOR['initial']  # horizontal hydraulic conductivity in m/day
        # if i == 1:
        #     k = np.ones_like(idomain[0]) * 0.001  # horizontal hydraulic conductivity in m/day
        # else:
        #     k = np.ones_like(idomain[0]) * 100  # horizontal hydraulic conductivity in m/day
        all_k.append(k)
    k_hor = np.array(all_k)  # horizontal hydraulic conductivity

    # save
    fn_out['npf_k'] = []
    for i in range(NLAY):
        ilay = i + 1  # layer number starts from 1
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.npf_k_layer{ilay}.txt').as_posix()
        fn_out['npf_k'].append(tomf6input(fn))
        np.savetxt(fn, k_hor[i])
    return fn_out, k_hor

def sto_ss_setup(idomain, fn_out):
    sto_ss = np.ones_like(idomain) * SS_PRIOR['initial']

    # save
    fn_out['sto_ss'] = []
    for i in range(NLAY):
        ilay = i + 1  # layer number starts from 1
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.sto_ss_layer{ilay}.txt').as_posix()
        fn_out['sto_ss'].append(tomf6input(fn))
        np.savetxt(fn, sto_ss[i])
    return fn_out

def rch_setup(idomain, fn_out):
    fn_out['rch'] = {}
    for i in range(NPER):
        if i in [0, 2]:  # steady state periods
            recharge_rate = RCH['initial_rf'] + RCH['initial_mbr']
        else:
            recharge_rate = RCH['initial_mbr']
        recharge = np.ones_like(idomain[0]) * recharge_rate * idomain[0]
        rch_fn = Path(MODEL_DIR, f'{MODEL_NAME}.rcha_recharge_{i}.txt').as_posix()
        fn_out['rch'][i] = tomf6input(rch_fn)
        np.savetxt(rch_fn, recharge) 
    return fn_out

def pk4_ts(start, end):
    # time series for riv elevations ###########################################################
    # get 1 day before start
    # start = start - pd.Timedelta(days=1)
    gw_ts = pd.read_csv(PK4_TS)
    gw_ts['DateTime'] = gw_ts['DateTime'].apply(pd.to_datetime, format = "%Y-%m-%d")
    gw_ts.set_index('DateTime', inplace=True)
    gw_ts['level_mRL'] = gw_ts['LEVEL']
    gw_ts['sm_level_mRL'] = gw_ts['level_mRL'].rolling(window=21, center=True, min_periods=1).mean()

    # get specific dates
    gw_ts = gw_ts[(gw_ts.index >= start) & (gw_ts.index < end)]
    gw_ts['abs_val'] = gw_ts['level_mRL'] - gw_ts['level_mRL'].iloc[0]
    gw_ts['sm_abs_val'] = gw_ts['abs_val'].rolling(window=5, center=True, min_periods=1).mean()

    gw_ts.index = np.arange(1, len(gw_ts['level_mRL']) + 1)

    # save ts file
    gw_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.pk4_level.csv')
    gw_ts.to_csv(gw_ts_fn)
    return gw_ts

def poukawa_ts(start, end):
    # time series for riv elevations ###########################################################
    sw_ts = pd.read_csv(POUKAWA_Q)
    sw_ts['DateTime'] = sw_ts['Poukawa Stream at Douglas Road (x)'].apply(pd.to_datetime, format = "%Y-%m-%d %H:%M:%S", dayfirst=True)
    sw_ts.set_index('DateTime', inplace=True)
    sw_ts['m3d'] = (sw_ts['Poukawa Stream at Douglas Road (y)'] * 84.6)

    # get specific dates
    sw_ts = sw_ts[(sw_ts.index > start) & (sw_ts.index < end)]
    sw_ts['sm_m3d'] = sw_ts['m3d'].rolling(window=5, center=True, min_periods=1).mean()
    sw_ts.index = np.arange(2, len(sw_ts['m3d']) + 2)

    # save ts file
    sw_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.poukawa_flow_m3d.csv')
    sw_ts.to_csv(sw_ts_fn)
    return sw_ts

def awanui_ts(start, end):
    # time series for riv elevations ###########################################################
    sw_ts = pd.read_csv(AWANUI_Q)
    sw_ts['DateTime'] = sw_ts['Awanui Stream at Flume (x)'].apply(pd.to_datetime, format = "%Y-%m-%d %H:%M:%S", dayfirst=True)
    sw_ts.set_index('DateTime', inplace=True)
    sw_ts['m3d'] = (sw_ts['Awanui Stream at Flume (y)'] * 84.6)

    # get specific dates
    sw_ts = sw_ts[(sw_ts.index > start) & (sw_ts.index < end)]
    sw_ts['sm_m3d'] = sw_ts['m3d'].rolling(window=21, center=True, min_periods=1).mean()
    sw_ts.index = np.arange(2, len(sw_ts['m3d']) + 2)

    # save ts file
    sw_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.awanui_flow_m3d.csv')
    sw_ts.to_csv(sw_ts_fn)
    return sw_ts