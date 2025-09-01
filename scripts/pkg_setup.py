from pathlib import Path
import numpy as np
import gridit as gi
import pandas as pd

from setup import *
from utils import *



def dis_setup(fn_out):
    grid = gi.Grid.from_vector(DOMAIN, RES)
    # grid_gpd = grid.cell_geodataframe()
    # # grid_gpd.to_file(Path(SPATIAL_DIR, f'{MODEL_NAME}_grid.shp'), driver='ESRI Shapefile')

    # top & bottom
    top = grid.array_from_raster(TOP)
    bottom = grid.array_from_raster(BOTTOM)  # get the bottom elevation

    # open domain
    arr = ~grid.array_from_vector(DOMAIN).mask # want the binary version of the domain
    # arr = np.where(top.data > 15, 0, arr)

    idomain = np.stack([arr] * NLAY, axis=0)

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

    return grid, idomain, top, nrow, ncol, delr, delc, fn_out

def ghb_aw_setup(grid, idomain, start, end, fn_out):
    riv_mask = ~grid.array_from_vector(DRAINS).mask * idomain[0]
    zones = grid.array_from_vector(DRAIN_ZONES, attribute='elev_ss_0')
    riv_stage_zones = ~zones.mask
    riv_rbot = grid.array_from_raster(TOP, resampling='min').data
    # adjust riv elevations absed on survey data
    riv_stage_arr = np.where(riv_stage_zones, zones.data, riv_rbot + 0.3) * riv_mask
    
    # to dataframe
    riv_kper0 = extract_value_with_indices(
        riv_stage_arr, layer=0, val_col='head', mask_value=0
        )
    riv_kper0['cond'] = RES**2 * 1

    # time series for riv elevations ###########################################################
    riv_ts = pd.read_csv(AWANUI_TS)
    riv_ts['DateTime'] = riv_ts['Awanui Stream at Flume (x)'].apply(pd.to_datetime, format = "%Y-%m-%d %H:%M:%S", dayfirst=True)
    riv_ts.set_index('DateTime', inplace=True)
    riv_ts['level_mRL'] = (riv_ts['Awanui Stream at Flume (y)'] * 0.001) - 10
    riv_ts['sm_level_mRL'] = riv_ts['level_mRL'].rolling(window=21, center=True, min_periods=1).mean()

    # get specific dates
    riv_ts = riv_ts[(riv_ts.index >= start) & (riv_ts.index <= end)]
    riv_ts['abs_val'] = riv_ts['level_mRL'] - riv_ts['level_mRL'].iloc[0]
    riv_ts['sm_abs_val'] = riv_ts['abs_val'].rolling(window=21, center=True, min_periods=1).mean()

    # riv absolute values during recession
    riv_ts_values = riv_ts['sm_abs_val'].values

    rivstage_dat = {}
    for row in range(len(riv_kper0)):
        index = riv_kper0.at[row, 'index']
        col = f's_1_{index[1]+1}_{index[2]+1}'
        base_elev = riv_kper0.at[row, 'head']
        rivstage_dat[col] = (base_elev + riv_ts_values).tolist()
    # convert to dataframe
    rivstage_ts = pd.DataFrame(
        rivstage_dat, 
        index=np.arange(1, len(riv_ts_values) + 1)
        )

    riv_kper1 = riv_kper0[['index']].copy()
    riv_kper1['head'] = riv_kper1['index'].apply(lambda x: f's_1_{x[1]+1}_{x[2]+1}')
    riv_kper1['cond'] = riv_kper0['cond']

    # kper2: PAST ############################################################
    riv_mask = ~grid.array_from_vector(DRAINS_PAST).mask
    wetland_mask = ~grid.array_from_vector(WETLANDA_PAST).mask
    past_arr = np.where(wetland_mask + riv_mask > 0, 1, 0)
    top = grid.array_from_raster(TOP).data
    wetland_top = top * wetland_mask
    top_med = np.median(wetland_top[wetland_top > 0])
    past_h = np.where(wetland_mask, top_med, np.where(riv_mask, top, 0)) * idomain[0]
    
    # to dataframe
    riv_kper2 = extract_value_with_indices(
        past_h, layer=0, val_col='head', mask_value=0
        )
    riv_kper2['cond'] = RES**2 * 1
    # same as kper2 (for parameterization)
    riv_kper3 = riv_kper2.copy()

    # save
    fn_out['ghb_aw'] = {}
    for i, dat in enumerate([riv_kper0, riv_kper1, riv_kper2, riv_kper3]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghbaw_stress_period_data_{i}.txt')
        fn_out['ghb_aw'][i] = tomf6input(fn, list=True)
        savedf2txt(dat, filename=fn, col_order=['head', 'cond'])
    # save ts file
    riv_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.awanui_stage.csv')
    rivstage_ts.index.name = '#time'
    rivstage_ts.to_csv(riv_ts_fn)
    fn_out['ghb_aw_ts'] = tomf6tsinput(riv_ts_fn, rivstage_ts)

    return fn_out

def ghb_spring_setup(grid, idomain, start, end, fn_out):
    spring_mask = ~grid.array_from_vector(SPRING_DRAIN).mask * idomain[0]
    zones = grid.array_from_vector(DRAIN_ZONES, attribute='elev_ss_0')
    spring_stage_zones = ~zones.mask
    spring_rbot = grid.array_from_raster(TOP, resampling='min').data
    # adjust riv elevations absed on survey data
    spring_stage_arr = np.where(spring_stage_zones, zones.data, spring_rbot + 0.3) * spring_mask
    
    # to dataframe
    spring_kper0 = extract_value_with_indices(
        spring_stage_arr, layer=0, val_col='head', mask_value=0
        )
    spring_kper0['cond'] = RES**2 * 1

    # time series for riv elevations ###########################################################
    spring_ts = pd.read_csv(SPRING_TS)
    spring_ts['DateTime'] = spring_ts['DateTime'].apply(pd.to_datetime, format = "%Y-%m-%d")
    spring_ts.set_index('DateTime', inplace=True)
    spring_ts['level_mRL'] = spring_ts['LEVEL']
    spring_ts['sm_level_mRL'] = spring_ts['level_mRL'].rolling(window=21, center=True, min_periods=1).mean()

    # get specific dates
    spring_ts = spring_ts[(spring_ts.index >= start) & (spring_ts.index <= end)]
    spring_ts['abs_val'] = spring_ts['level_mRL'] - spring_ts['level_mRL'].iloc[0]
    spring_ts['sm_abs_val'] = spring_ts['abs_val'].rolling(window=5, center=True, min_periods=1).mean()

    # assert len(riv_ts_resampled) == NSTEPS

    # riv absolute values during recession
    spring_ts_values = spring_ts['sm_abs_val'].values

    spring_dat = {}
    for row in range(len(spring_kper0)):
        index = spring_kper0.at[row, 'index']
        col = f's_1_{index[1]+1}_{index[2]+1}'
        base_elev = spring_kper0.at[row, 'head']
        spring_dat[col] = (base_elev + spring_ts_values).tolist()
    # convert to dataframe
    spring_h_ts = pd.DataFrame(
        spring_dat, 
        index=np.arange(1, len(spring_ts_values) + 1)
        )
    # make sure spring heads are like those observed
    spring_arr = grid.array_from_vector(SPRING).data
    spring_loc = extract_value_with_indices(
        spring_arr, layer=0, val_col='head', mask_value=0
        )['index'].values[0]
    spring_idx = f's_1_{spring_loc[1]+1}_{spring_loc[2]+1}'
    spring_h_ts[spring_idx] = spring_ts['level_mRL'].values
    
    # create spring_kper1
    spring_kper1 = spring_kper0[['index']].copy()
    spring_kper1['head'] = spring_kper1['index'].apply(lambda x: f's_1_{x[1]+1}_{x[2]+1}')
    spring_kper1['cond'] = spring_kper0['cond']

    # make empty kper2 and kper3 (for parameterization)
    spring_kper2 = pd.DataFrame()
    spring_kper3 = pd.DataFrame()

    # save
    fn_out['ghb_sp'] = {}
    for i, dat in enumerate([spring_kper0, spring_kper1, spring_kper2, spring_kper3]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghbspr_stress_period_data_{i}.txt')
        fn_out['ghb_sp'][i] = tomf6input(fn, list=True)
        savedf2txt(dat, filename=fn, col_order=['head', 'cond'])
    # save ts file
    spring_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.spring_stage.csv')
    spring_h_ts.index.name = '#time'
    spring_h_ts.to_csv(spring_ts_fn)
    # spring_h_ts.to_csv(spring_ts_fn, header=False)
    fn_out['ghb_sp_ts'] = tomf6tsinput(spring_ts_fn, spring_h_ts)

    return fn_out

def ghb_pw_setup(grid, idomain, start, end, fn_out):
    top_min = grid.array_from_raster(TOP, resampling='min')
    ghb_pw_arr = ~grid.array_from_vector(POUKAWA_BOUNDARY).mask
    in_idomarr, idom_i = get_interior_indices(idomain[0], layer=0)  # idomainget interior indices of the mbr area
    ghb_pw_active = np.logical_and(ghb_pw_arr.data, in_idomarr)  # mbr area that is active in the model domain
    ghb_pw_active = ghb_pw_active * top_min.data
    ghb_pw_indices = get_indices(ghb_pw_active, layer=0, value=True)

    ghb_pw_kper0 = pd.DataFrame({'index': [i[0] for i in ghb_pw_indices]})
    ghb_pw_kper0['head'] = [i[1] for i in ghb_pw_indices]  # head for Poukawa boundary
    ghb_pw_kper0['cond'] = RES**2 * 1

    # TS ##########################
    pw_ts = pd.read_csv(POUKAWA_TS)
    pw_ts['DateTime'] = pw_ts['Poukawa Stream at Douglas Road (x)'].apply(pd.to_datetime, format = "%Y-%m-%d %H:%M:%S", dayfirst=True)
    pw_ts['level_mRL'] = (pw_ts['Poukawa Stream at Douglas Road (y)'] * 0.001) - 10
    pw_ts['sm_level_mRL'] = pw_ts['level_mRL'].rolling(window=21, center=True, min_periods=1).mean()
    pw_ts.set_index('DateTime', inplace=True)

    # get specific dates
    pw_ts = pw_ts[(pw_ts.index >= start) & (pw_ts.index <= end)]
    pw_ts['abs_val'] = pw_ts['level_mRL'] - pw_ts['level_mRL'].iloc[0]
    pw_ts['sm_abs_val'] = pw_ts['abs_val'].rolling(window=5, center=True, min_periods=1).mean()

    # resample to match number of time steps
    # pw_ts_resampled = riv_ts.resample(f'{DAYS_P_STEP}D', origin=start).first()

    # riv absolute values during recession
    pw_ts_values = pw_ts['sm_abs_val'].values

    pwelev_dat = {}
    for row in range(len(ghb_pw_kper0)):
        index = ghb_pw_kper0.at[row, 'index'] 
        col = f'h_1_{index[1]+1}_{index[2]+1}'
        base_elev = ghb_pw_kper0.at[row, 'head']
        pwelev_dat[col] = (base_elev + pw_ts_values).tolist()
    # convert to dataframe
    pwelev_ts = pd.DataFrame(
        pwelev_dat,
        index=np.arange(1, len(pw_ts_values) + 1),
        )

    ghb_pw_kper1 = ghb_pw_kper0[['index']].copy()
    ghb_pw_kper1['head'] = ghb_pw_kper1['index'].apply(lambda x: f'h_1_{x[1]+1}_{x[2]+1}')
    ghb_pw_kper1['cond'] = ghb_pw_kper0['cond']

    # kper3 - PAST SS
    top_elev = grid.array_from_raster(TOP) # heads slightly higher than min_top (no drain)
    ghb_pw_arr = ~grid.array_from_vector(PW_PAST).mask
    in_idomarr, idom_i = get_interior_indices(idomain[0], layer=0)  # get edges of idomain
    ghb_pw_active = np.logical_and(ghb_pw_arr.data, in_idomarr)
    ghb_pw_active = ghb_pw_active * top_elev.data
    ghb_pw_indices = get_indices(ghb_pw_active, layer=0, value=True)

    ghb_pw_kper2 = pd.DataFrame({'index': [i[0] for i in ghb_pw_indices]})
    ghb_pw_kper2['head'] = [i[1] for i in ghb_pw_indices]  # head for Poukawa boundary
    ghb_pw_kper2['cond'] = RES**2 * 1

    # kper4 - PAST TS (same: keep heads same as kper3)
    ghb_pw_kper3 = ghb_pw_kper2.copy()

    #save
    fn_out['ghb_pw'] = {}
    for i, dat in enumerate([ghb_pw_kper0, ghb_pw_kper1, ghb_pw_kper2, ghb_pw_kper3]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_pw_stress_period_data_{i}.txt')
        fn_out['ghb_pw'][i] = tomf6input(fn ,list=True)
        savedf2txt(dat, filename=fn, col_order=['head', 'cond'])
    
    # save ts file
    pwelev_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_pw_heads.csv')
    pwelev_ts.index.name = '#time'
    pwelev_ts.to_csv(pwelev_ts_fn)
    fn_out['ghb_pw_ts'] = tomf6tsinput(pwelev_ts_fn, pwelev_ts)
    return fn_out, ghb_pw_kper0

def ghb_conf_setup(grid, idomain, start, end, fn_out):
    conf_arr = ~grid.array_from_vector(CONF_AREA_ACTIVE).mask * idomain[0]
    conf_kper0 = extract_value_with_indices(
        conf_arr, layer=0, val_col='head', mask_value=0
        )
    conf_kper0['head'] = 13.5  # head for confining layer
    conf_kper0['cond'] = RES**2 * 1

    # TS ##########################
    conf_ts = pd.read_csv(CONF_TS)
    conf_ts['DateTime'] = conf_ts['DateTime'].apply(pd.to_datetime, format = "%d/%m/%Y %H:%M")
    conf_ts.set_index('DateTime', inplace=True)
    conf_ts['level_mRL'] = conf_ts['LEVEL'] + 10.755
    #resample to daily
    conf_ts = conf_ts[(conf_ts.index <= end)]
    conf_ts = conf_ts[['level_mRL']].resample('D').mean()
    conf_ts['level_mRL'] = conf_ts['level_mRL'].interpolate(method='time')
    conf_ts = conf_ts[(conf_ts.index >= start) & (conf_ts.index <= end)]
    # conf_ts['abs_val'] = conf_ts['level_mRL'] - conf_ts['level_mRL'].iloc[0]
    # conf_ts['sm_abs_val'] = conf_ts['abs_val'].rolling(window=5, center=True, min_periods=1).mean()

    # convert to dataframe
    conf_ts_df = conf_ts['level_mRL'].reset_index(drop=True)
    conf_ts_df.name = 'heads_ts'
    conf_ts_df.index = conf_ts_df.index + 1

    conf_kper1 = conf_kper0[['index']].copy()
    conf_kper1['head'] = 'heads_ts'
    conf_kper1['cond'] = conf_kper0['cond']

    # PAST
    conf_kper2 = conf_kper0.copy()
    conf_kper3 = conf_kper0.copy()
    
    #save
    fn_out['ghb_conf'] = {}
    for i, dat in enumerate([conf_kper0, conf_kper1, conf_kper2, conf_kper3]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_conf_stress_period_data_{i}.txt')
        fn_out['ghb_conf'][i] = tomf6input(fn ,list=True)
        savedf2txt(dat, filename=fn, col_order=['head', 'cond'])
    
    # save ts file
    conf_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.ghb_conf_heads.csv')
    conf_ts_df.index.name = '#time'
    conf_ts_df.to_csv(conf_ts_fn)
    fn_out['ghb_conf_ts'] = tomf6tsinput(conf_ts_fn, conf_ts_df)
    return fn_out

def wel_mbr_setup(grid, idomain, fn_out, ghb_pw_kper0):
    mbr_arr = ~grid.array_from_vector(MBR).mask
    mbr_indices = []
    for i in range(NLAY):
        in_idomarr, idom_i = get_interior_indices(idomain[i], layer=i)  # idomainget interior indices of the mbr area
        mbr_active = np.logical_and(mbr_arr.data, in_idomarr)  # mbr area that is active in the model domain
        mbr_indices.extend(get_indices(mbr_active, layer=i))
    
    mbr_df = pd.DataFrame({'index': mbr_indices})
    mbr_df['q'] = 2 # m3/d
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

def npf_setup(idomain, fn_out):
    all_k = []
    for i in range(NLAY):
        if i == 1:
            k = np.ones_like(idomain[0]) * 0.001  # horizontal hydraulic conductivity in m/day
        else:
            k = np.ones_like(idomain[0]) * 100  # horizontal hydraulic conductivity in m/day
        all_k.append(k)
    k_hor = np.array(all_k)  # horizontal hydraulic conductivity

    # save
    fn_out['npf_k'] = []
    for i in range(NLAY):
        ilay = i + 1  # layer number starts from 1
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.npf_k_layer{ilay}.txt').as_posix()
        fn_out['npf_k'].append(tomf6input(fn))
        np.savetxt(fn, k_hor[i])
    return fn_out

def sto_ss_setup(idomain, fn_out):
    sto_ss = np.ones_like(idomain) * 1e-4

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
        if i in [0, 2]:
            recharge = np.ones_like(idomain[0]) * 0.0001 * idomain[0]
        else:
            recharge = np.zeros_like(idomain[0])
        rch_fn = Path(MODEL_DIR, f'{MODEL_NAME}.rcha_recharge_{i}.txt').as_posix()
        fn_out['rch'][i] = tomf6input(rch_fn)
        np.savetxt(rch_fn, recharge)  # save bottom elevation for each layer
    return fn_out