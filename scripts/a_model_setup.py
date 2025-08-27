import matplotlib.pyplot as plt
import os
from pathlib import Path
import shutil
import gridit as gi
import numpy as np
import flopy as fp
import pandas as pd
from scipy.ndimage import uniform_filter

import helpers as helpers

from utils import *
from setup import *

def main():

    # print working directory
    print(f'Working directory: {os.getcwd()}')
    print(f'Model directory: {os.path.abspath(MODEL_DIR)}')

    # create directories if they do not exist
    for d in [SPATIAL_DIR, FIG_DIR]:
        if not os.path.exists(d):
            os.makedirs(d)
    # if dir exists delete it
    if os.path.exists(MODEL_DIR):
        shutil.rmtree(MODEL_DIR)
    # create model directory
    os.makedirs(MODEL_DIR)  # create model directory
    
    # copy all the contents of bin into model directory
    if os.path.exists(BIN_DIR):
        if os.name == 'nt':  # if on Windows, copy files
            os_bin = Path(BIN_DIR, 'windows')
        elif os.name == 'posix':  # if on Linux or MacOS, copy files
            os_bin = Path(BIN_DIR, 'linux')
        else:
            raise ValueError(f'Unsupported OS: {os.name}. Please check the BIN_DIR path.')

        shutil.copytree(os_bin, MODEL_DIR, dirs_exist_ok=True)  # copy bin directory to model directory
        # Make executable files executable
        import stat
        for root, dirs, files in os.walk(MODEL_DIR):
            for file in files:
                # Make all files that could be executables executable
                if file in ['mf6']:
                    file_path = Path(root, file)
                    current_permissions = os.stat(file_path).st_mode
                    os.chmod(file_path, current_permissions | stat.S_IEXEC | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
                    print(f"Made executable: {file_path}")
    else:
        print(f'Bin directory {BIN_DIR} does not exist. Please check the path.')

    fn_out = {}
    # TDIS ####################################################################
    start = pd.to_datetime(START_DATE)
    end = pd.to_datetime(END_DATE)
    period_days = (end - start).days
    PERIOD_DAYS = period_days
    NSTEPS = period_days  # number of time steps
    print(f'Simulation period: {PERIOD_DAYS} days with {NSTEPS} time steps.')
    

    # DIS ####################################################################
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


    # RIV #################################################################################
    riv_mask = ~grid.array_from_vector(DRAINS).mask * idomain[0]
    zones = grid.array_from_vector(DRAIN_ZONES, attribute='elev_ss_0')
    riv_stage_zones = ~zones.mask * riv_mask
    riv_rbot = grid.array_from_raster(TOP, resampling='min').data * riv_mask

    # adjust riv elevations absed on survey data
    riv_stage_arr = np.where(riv_stage_zones, zones.data, riv_rbot + 0.3)
    riv_rbot_arr = np.where(riv_rbot < riv_stage_arr, riv_rbot, riv_stage_arr - 0.3)
    
    # extract non-NaN values from the riv elevation array
    riv_stage = extract_value_with_indices(
        riv_stage_arr, layer=0, val_col='stage', mask_value=0
        )
    riv_rbot = extract_value_with_indices(
        riv_rbot_arr, layer=0, val_col='rbot', mask_value=0
        )
    riv_kper0 = pd.merge(riv_stage, riv_rbot, on='index', how='inner')  # merge the two dataframes on the index column
    riv_kper0['cond'] = RES**2 * 1

    # # where rbot > stage, set rbot to stage - 0.1
    # riv_kper0['rbot'].clip(lower=riv_kper0['stage']-0.1, inplace=True)

    # time series for riv elevations ###########################################################
    riv_ts = pd.read_csv(AWANUI_TS)
    riv_ts['DateTime'] = riv_ts['Awanui Stream at Flume (x)'].apply(pd.to_datetime, format = "%Y-%m-%d %H:%M:%S", dayfirst=True)
    riv_ts['level_mRL'] = (riv_ts['Awanui Stream at Flume (y)'] * 0.001) - 10

    riv_ts['sm_level_mRL'] = riv_ts['level_mRL'].rolling(window=21, center=True, min_periods=1).mean()

    # get specific dates
    riv_ts = riv_ts[(riv_ts['DateTime'] >= start) & (riv_ts['DateTime'] <= end)].reset_index(drop=True)
    assert len(riv_ts) == NSTEPS, f"riv time series length {len(riv_ts)} does not match number of time steps {NSTEPS}."
    riv_ts['abs_val'] = riv_ts['level_mRL'] - riv_ts['level_mRL'].iloc[0]
    riv_ts['sm_abs_val'] = riv_ts['abs_val'].rolling(window=21, center=True, min_periods=1).mean()

    # riv absolute values during recession
    riv_ts_values = riv_ts['sm_abs_val'].values

    rivstage_dat = {}
    for row in range(len(riv_kper0)):
        index = riv_kper0.at[row, 'index']
        col = f's_1_{index[1]+1}_{index[2]+1}'
        base_elev = riv_kper0.at[row, 'stage']
        rivstage_dat[col] = (base_elev + riv_ts_values).tolist()
    # convert to dataframe
    rivstage_ts = pd.DataFrame(
        rivstage_dat, 
        index=np.arange(1, NSTEPS + 1)
        )

    riv_kper1 = riv_kper0[['index']].copy()
    riv_kper1['stage'] = riv_kper1['index'].apply(lambda x: f's_1_{x[1]+1}_{x[2]+1}')
    riv_kper1['cond'] = riv_kper0['cond']
    riv_kper1['rbot'] = riv_kper0['rbot']

    # manual adjustment to rbot where rbot > stage
    for i, row in riv_kper1.iterrows():
        if np.any(row['rbot'] >= rivstage_ts[row['stage']]):
            riv_kper1.at[i, 'old_rbot'] = riv_kper1.at[i, 'rbot']
            riv_kper1.at[i, 'rbot'] = rivstage_ts[row['stage']].min() - 0.1

    # save
    fn_out['riv_aw'] = {}
    for i, dat in enumerate([riv_kper0, riv_kper1]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.riv_stress_period_data_{i}.txt')
        fn_out['riv_aw'][i] = tomf6input(fn, list=True)
        savedf2txt(dat, filename=fn, col_order=['stage', 'cond', 'rbot'])
    # save ts file
    riv_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.riv_stage.csv')
    rivstage_ts.to_csv(riv_ts_fn, header=False)
    fn_out['riv_aw_ts'] = tomf6tsinput(riv_ts_fn, rivstage_ts)
    
    # CHD (POUKAWA) ##########################################################################
    top_min = grid.array_from_raster(TOP, resampling='min')
    chd_pw_arr = ~grid.array_from_vector(POUKAWA_BOUNDARY).mask
    in_idomarr, idom_i = get_interior_indices(idomain[0], layer=0)  # idomainget interior indices of the mbr area
    chd_pw_active = np.logical_and(chd_pw_arr.data, in_idomarr)  # mbr area that is active in the model domain
    chd_pw_active = chd_pw_active * top_min.data
    chd_pw_indices = get_indices(chd_pw_active, layer=0, value=True)

    chd_pw_kper0 = pd.DataFrame({'index': [i[0] for i in chd_pw_indices]})
    chd_pw_kper0['head'] = [i[1] for i in chd_pw_indices]  # head for Poukawa boundary

    # TS ##########################
    pw_ts = pd.read_csv(POUKAWA_TS)
    pw_ts['DateTime'] = pw_ts['Poukawa Stream at Douglas Road (x)'].apply(pd.to_datetime, format = "%Y-%m-%d %H:%M:%S", dayfirst=True)
    pw_ts['level_mRL'] = (pw_ts['Poukawa Stream at Douglas Road (y)'] * 0.001) - 10

    pw_ts['sm_level_mRL'] = pw_ts['level_mRL'].rolling(window=21, center=True, min_periods=1).mean()

    # get specific dates
    pw_ts = pw_ts[(pw_ts['DateTime'] >= start) & (pw_ts['DateTime'] <= end)].reset_index(drop=True)
    assert len(pw_ts) == NSTEPS, f"riv time series length {len(pw_ts)} does not match number of time steps {NSTEPS}."
    pw_ts['abs_val'] = pw_ts['level_mRL'] - pw_ts['level_mRL'].iloc[0]
    pw_ts['sm_abs_val'] = pw_ts['abs_val'].rolling(window=5, center=True, min_periods=1).mean()

    # riv absolute values during recession
    pw_ts_values = pw_ts['sm_abs_val'].values

    pwelev_dat = {}
    for row in range(len(chd_pw_kper0)):
        index = chd_pw_kper0.at[row, 'index'] 
        col = f'h_1_{index[1]+1}_{index[2]+1}'
        base_elev = chd_pw_kper0.at[row, 'head']
        pwelev_dat[col] = (base_elev + pw_ts_values).tolist()
    # convert to dataframe
    pwelev_ts = pd.DataFrame(pwelev_dat, index=np.arange(1, NSTEPS + 1))

    chd_pw_kper1 = chd_pw_kper0[['index']].copy()
    chd_pw_kper1['head'] = chd_pw_kper1['index'].apply(lambda x: f'h_1_{x[1]+1}_{x[2]+1}')

    #save
    fn_out['chd_pw'] = {}
    for i, dat in enumerate([chd_pw_kper0, chd_pw_kper1]):
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.chd_pw_stress_period_data_{i}.txt')
        fn_out['chd_pw'][i] = tomf6input(fn ,list=True)
        savedf2txt(dat, filename=fn, col_order=['head'])
    
    # save ts file
    pwelev_ts_fn = Path(MODEL_DIR, f'{MODEL_NAME}.chd_heads.csv')
    pwelev_ts.to_csv(pwelev_ts_fn, header=False)
    fn_out['chd_pw_ts'] = tomf6tsinput(pwelev_ts_fn, pwelev_ts)

    
    # WEL (MBR) ##########################################################################
    mbr_arr = ~grid.array_from_vector(MBR).mask
    mbr_indices = []
    for i in range(NLAY-2):  # remove mbr from bottom 2 layers
        in_idomarr, idom_i = get_interior_indices(idomain[i], layer=i)  # idomainget interior indices of the mbr area
        mbr_active = np.logical_and(mbr_arr.data, in_idomarr)  # mbr area that is active in the model domain
        mbr_indices.extend(get_indices(mbr_active, layer=i))
    
    mbr_df = pd.DataFrame({'index': mbr_indices})
    mbr_df['q'] = 2 # m3/d
    mbr_df = mbr_df[~mbr_df['index'].isin(chd_pw_kper0['index'])]  # remove chd indices from mbr
    
    # save
    fn_out['wel_mbr'] = {}
    for i in range(NPER):
        mbr_fn = Path(MODEL_DIR, f'{MODEL_NAME}.wel_mbr_stress_period_data_{i}.txt')
        fn_out['wel_mbr'][0] = tomf6input(mbr_fn, list=True)
        savedf2txt(mbr_df, filename=mbr_fn, col_order=['q'])

    # WEL (INFLUX/OUTFLUX) ###################################################################

    # influx boundaries
    influx_arr = ~grid.array_from_vector(INFLUX_BOUNDARY).mask
    influx_indices = []
    for i in range(NLAY-2): # remove influx from bottom 2 layers
        in_idomarr, idom_i = get_interior_indices(idomain[i], layer=i)  # idomainget interior indices of the mbr area
        influx_active = np.logical_and(influx_arr.data, in_idomarr)  # mbr area that is active in the model domain
        influx_indices.extend(get_indices(influx_active, layer=i))
    influx_df = pd.DataFrame({'index': influx_indices})

    # outflux boundaries
    outflux_arr = ~grid.array_from_vector(OUTFLUX_BOUNDARY).mask
    outflux_indices = []
    for i in range(NLAY-2):
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

    # NPF: K ###################################################################
    all_k = []
    for i in range(NLAY):
        if i == 1:
            k = np.ones((nrow, ncol)) * 0.001  # horizontal hydraulic conductivity in m/day
        else:
            k = np.ones((nrow, ncol)) * 100  # horizontal hydraulic conductivity in m/day
        all_k.append(k)
    k_hor = np.array(all_k)  # horizontal hydraulic conductivity

    # save
    fn_out['npf_k'] = []
    for i in range(NLAY):
        ilay = i + 1  # layer number starts from 1
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.npf_k_layer{ilay}.txt').as_posix()
        fn_out['npf_k'].append(tomf6input(fn))
        np.savetxt(fn, k_hor[i])

    # STO: SS ###################################################################
    sto_ss = np.ones_like(idomain) * 1e-4

    # save
    fn_out['sto_ss'] = []
    for i in range(NLAY):
        ilay = i + 1  # layer number starts from 1
        fn = Path(MODEL_DIR, f'{MODEL_NAME}.sto_ss_layer{ilay}.txt').as_posix()
        fn_out['sto_ss'].append(tomf6input(fn))
        np.savetxt(fn, sto_ss[i])
    
    # RCH: Rainfall ################################################################
    fn_out['rch'] = {}
    for i in range(NPER):
        if i == 0:
            recharge = np.ones_like(idomain[0]) * 0.0001 * idomain[0]
        else:
            recharge = np.zeros_like(idomain[0])
        rch_fn = Path(MODEL_DIR, f'{MODEL_NAME}.rcha_recharge_{i}.txt').as_posix()
        fn_out['rch'][i] = tomf6input(rch_fn)
        np.savetxt(rch_fn, recharge)  # save bottom elevation for each layer

    # --------------------------------------------------------------------------
    # other model parameters
    init_h = np.stack([top] * NLAY)  # initial head, based on average riv elevation

    # 2 BUILD A MODEL -------------------------------------------------------

    sim = fp.mf6.MFSimulation(
        sim_name=MODEL_NAME, # name of simulation
        version='mf6', # version of MODFLOW
        exe_name=f'{MODEL_DIR}/mf6',
        sim_ws=MODEL_DIR, # path to workspace where all files are stored
        )

    tdis = fp.mf6.ModflowTdis(
        simulation=sim, # add to the simulation called sim (defined in prevous code cell)
        time_units='DAYS', 
        nper=NPER, # number of stress periods
        perioddata=[
            (1, 1, 1),
            (PERIOD_DAYS, NSTEPS, 1),
            ], # period length, number of steps, timestep multiplier
        )

    ims = fp.mf6.ModflowIms(
        simulation=sim, 
        complexity='COMPLEX',
        print_option='ALL'
        )

    gwf = fp.mf6.ModflowGwf(
        simulation=sim, 
        modelname=MODEL_NAME,
        model_nam_file=f"{MODEL_NAME}.nam",
        save_flows=True,
        )

    dis = fp.mf6.ModflowGwfdis(
        model=gwf, # add to groundwater flow model called gwf
        length_units='METERS', 
        nlay=NLAY, 
        nrow=nrow, 
        ncol=ncol,
        delr=delr, 
        delc=delc, 
        top=fn_out['dis']['top'],
        botm=fn_out['dis']['btm'],
        idomain=fn_out['dis']['idomain'],
        xorigin=grid.bounds[0],
        yorigin=grid.bounds[1],
        )
    
    sto = fp.mf6.ModflowGwfsto(
        model=gwf,
        ss=fn_out['sto_ss'], 
        steady_state={0: True},
        transient={1: True},
        )

    npf = fp.mf6.ModflowGwfnpf(
        model=gwf, #node property flow package
        save_specific_discharge=True, # save the specific discharge for every cell
        icelltype=0, # 0 means constant saturated thickness
        k=fn_out['npf_k'], # horizontal k value
        )

    ic = fp.mf6.ModflowGwfic(
        model=gwf, 
        strt=init_h, # initial head, only used for iterative solution in steady model (arbitrary)
        )

    rch = fp.mf6.ModflowGwfrcha(
        model=gwf,
        recharge=fn_out['rch'], # recharge file names for each layer
        pname='rch' # package name
    )

    riv_aw = fp.mf6.ModflowGwfriv(
        model=gwf, # add riv package to model gwf (created in previous code cell)
        timeseries=fn_out['riv_aw_ts'], # time series file for riv stages
        stress_period_data=fn_out['riv_aw'],
        pname='riv_aw', # package name
        save_flows=True,
        )

    chd_pw = fp.mf6.ModflowGwfchd(
        model=gwf,
        timeseries=fn_out['chd_pw_ts'], # time series file for chd heads
        stress_period_data=fn_out['chd_pw'],
        pname='chd_pw', # package name
        save_flows=True, # save flows for this package 
        )

    wel_mbr = fp.mf6.ModflowGwfwel(
        model=gwf,
        stress_period_data=fn_out['wel_mbr'],
        pname='wel_mbr' # package name
        )
    
    influx = fp.mf6.ModflowGwfwel(
        model=gwf,
        stress_period_data=fn_out['wel_influx'],
        pname='influx' # package name
        )
    outflux = fp.mf6.ModflowGwfwel(
        model=gwf,
        stress_period_data=fn_out['wel_outflux'],
        pname='outflux' # package name
        )

    oc = fp.mf6.ModflowGwfoc(
        model=gwf, # add output control to model gwf (created in previous code cell)
        budget_filerecord=f"{MODEL_NAME}.cbc", # file name where all budget output is stored
        head_filerecord=f"{MODEL_NAME}.hds", # file name where all head output is stored
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

    # --------------------------------------------------------

    print('Writing model files...')
    sim.write_simulation()  # write all model files to disk
    print('Running model...')
    success, _ = sim.run_simulation()  # run the model

    if not success:
        raise Exception("Model did not run successfully. Please check the output files for errors.")

    # OUTPUT PLOTS --------------------------------------------------------
    # visual check
    pmv = fp.plot.PlotMapView(model=gwf, layer=0) # create view of layer 0
    pmv.plot_array(top *idomain[0], masked_values=[1e30], alpha=0.5, cmap='viridis') # plot top elevation

    pmv.plot_bc(name='chd_pw', color='lightblue') # add 'chd' cells
    # pmv.plot_bc(name='mbr', color='orange') # add 'wells' cells
    pmv.plot_bc(name='riv_aw', color='blue') # add 'wells' cells
    pmv.plot_bc(name='influx', color='green') # add 'influx' cells
    pmv.plot_bc(name='outflux', color='red') # add 'outflux' cells
    # pmv.plot_inactive(color='lightgray', alpha=0.5) # plot inactive cells
    # pmv.plot_grid(colors='silver', lw=0.01); # add grid

    # save to figures
    plt.savefig(Path(FIG_DIR, f'{MODEL_NAME}_domain.png'), dpi=300, bbox_inches='tight') # save figure

    # --------------------------------------------------------------

    # plot heads
    idom_plt = np.where(idomain[0] == 1, np.nan, idomain[0])  # create a mask for the active domain
    head = gwf.output.head().get_data()
    for i in range(NLAY):
        
        pmv = fp.plot.PlotMapView(model=gwf, layer=i)

        pmv.plot_bc(name='chd_pw', color='purple') # add 'chd' cells
        pmv.plot_bc(name='mbr', color='orange') # add 'wells' cells
        pmv.plot_bc(name='drn_r', color='blue') # add 'wells' cells
        pmv.plot_bc(name='chd_conf', color='orange') # add 'chd' cells
        # pmv.plot_bc(name='chd_conf_out', color='gold') # add 'chd' cells
        pmv.plot_bc(name='influx', color='green') # add 'influx' cells
        pmv.plot_bc(name='outflux', color='red') # add 'outflux' cells
        
        # pmv.plot_inactive(color='lightgray', alpha=0.5) # plot inactive cells
        # pmv.plot_array(idom_plt, cmap='gray', alpha=0.2) # plot head array for layer 0
        # pmv.plot_grid(colors='silver', lw=0.1); # add grid
        cs = pmv.contour_array(head[i], linewidths=1, colors='k') # contour plot of heads
        plt.clabel(cs, fmt='%1.1f'); # add contour labels with one decimal place

        # save to figures
        plt.savefig(Path(FIG_DIR, f'{MODEL_NAME}_heads{i}.png'), dpi=300, bbox_inches='tight') # save figure
        plt.close()

    # -------------------------------------------------------------
    # create a cross section grid & plot heads

    cross_col = 15 #which row to show cross section
    crossview = fp.plot.PlotCrossSection(model=gwf, line={'column': cross_col})
    strtArray = crossview.plot_array(head, masked_values=[1e30], alpha=0.8) # plot the active domain
    crossview.plot_grid(colors='black', lw=1); # add grid
    crossview.plot_inactive(color='lightgray') # plot inactive cells
    cb = plt.colorbar(strtArray, shrink=0.5) # add color bar
    # strtArray = crossview.plot_array(head, masked_values=[1e30], alpha = 0.5) # plot the array of heads in cross section
    plt.savefig(Path(FIG_DIR, f'{MODEL_NAME}_xsection_col{cross_col}.png'), dpi=300, bbox_inches='tight') # save figure
    plt.close()  # close the figure to avoid memory issues

    # TEST OBS ------------------------------------------------------------
    os.chdir(MODEL_DIR)  # change directory to model directory
    sample_path = os.path.relpath(SAMPLES, MODEL_DIR)
    helpers.extract_heads_and_budget(model_name=f'{MODEL_NAME}')  # extract heads and budget from model output
    helpers.extract_spring_obs(gwf=None, model_name=f'{MODEL_NAME}', samples_path=f'{sample_path}')  # extract spring 

if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations