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
from pkg_setup import *

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
    PERLEN = (end - start).days

    print(f'Simulation period: {PERLEN} days with {NSTEPS} time steps.')
    
    # dis
    grid, idomain, top, nrow, ncol, delr, delc, fn_out = dis_setup(fn_out)
    # ghbs
    fn_out = ghb_aw_setup(grid, idomain, start, end, fn_out)
    fn_out = ghb_spring_setup(grid, idomain, start, end, fn_out)
    fn_out, ghb_pw_kper0 = ghb_pw_setup(grid, idomain, start, end, fn_out)
    fn_out = ghb_conf_setup(grid, idomain, start, end, fn_out)
    # wel
    fn_out, mbr_df = wel_mbr_setup(grid, idomain, fn_out, ghb_pw_kper0)
    fn_out = wel_inout_setup(grid, idomain, mbr_df, fn_out)
    # npf
    fn_out = npf_setup(idomain, fn_out)
    # sto
    fn_out = sto_ss_setup(idomain, fn_out)
    # rch
    # fn_out = rch_setup(idomain, fn_out)
    #obs - save file
    pk4_h = pk4_ts(start, end)
    aw_Q = awanui_ts(start, end)
    pw_Q = poukawa_ts(start, end)

    
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
            (PERLEN, NSTEPS, 1),
            (1, 1, 1),
            (PERLEN, NSTEPS, 1)
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
        steady_state={
            0: True,
            2: True,
            },
        transient={
            1: True,
            3: True,
            },
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

    # rch = fp.mf6.ModflowGwfrcha(
    #     model=gwf,
    #     recharge=fn_out['rch'], # recharge file names for each layer
    #     save_flows=True, # save flows for this package
    #     pname='rch' # package name
    # )

    ghb_aw = fp.mf6.ModflowGwfghb(
        model=gwf, # add riv package to model gwf (created in previous code cell)
        timeseries=fn_out['ghb_aw_ts'], # time series file for riv stages
        stress_period_data=fn_out['ghb_aw'],
        pname='ghb_aw', # package name
        save_flows=True,
        )
    
    ghb_conf = fp.mf6.ModflowGwfghb(
        model=gwf,
        timeseries=fn_out['ghb_conf_ts'], # time series file for chd heads
        stress_period_data=fn_out['ghb_conf'],
        pname='ghb_conf', # package name
        save_flows=True, # save flows for this package 
        )

    ghb_pw = fp.mf6.ModflowGwfghb(
        model=gwf,
        timeseries=fn_out['ghb_pw_ts'], # time series file for chd heads
        stress_period_data=fn_out['ghb_pw'],
        pname='ghb_pw', # package name
        save_flows=True, # save flows for this package 
        )
    
    ghb_sp = fp.mf6.ModflowGwfghb(
        model=gwf,
        timeseries=fn_out['ghb_sp_ts'], # time series file for spring heads
        stress_period_data=fn_out['ghb_sp'],
        pname='ghb_sp', # package name
        save_flows=True, # save flows for this package 
        )

    wel_mbr = fp.mf6.ModflowGwfwel(
        model=gwf,
        stress_period_data=fn_out['wel_mbr'],
        pname='wel_mbr', # package name
        save_flows=True,
        )
    
    influx = fp.mf6.ModflowGwfwel(
        model=gwf,
        stress_period_data=fn_out['wel_influx'],
        pname='influx', # package name
        save_flows=True,
        )
    outflux = fp.mf6.ModflowGwfwel(
        model=gwf,
        stress_period_data=fn_out['wel_outflux'],
        pname='outflux', # package name
        save_flows=True,
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
    # pmv = fp.plot.PlotMapView(model=gwf, layer=0) # create view of layer 0
    # pmv.plot_array(top *idomain[0], masked_values=[1e30], alpha=0.5, cmap='viridis') # plot top elevation

    # pmv.plot_bc(name='chd_pw', color='lightblue') # add 'chd' cells
    # # pmv.plot_bc(name='mbr', color='orange') # add 'wells' cells
    # pmv.plot_bc(name='ghb_aw', color='blue') # add 'wells' cells
    # pmv.plot_bc(name='influx', color='green') # add 'influx' cells
    # pmv.plot_bc(name='outflux', color='red') # add 'outflux' cells
    # # pmv.plot_inactive(color='lightgray', alpha=0.5) # plot inactive cells
    # # pmv.plot_grid(colors='silver', lw=0.01); # add grid

    # # save to figures
    # plt.savefig(Path(FIG_DIR, f'{MODEL_NAME}_domain.png'), dpi=300, bbox_inches='tight') # save figure

    # --------------------------------------------------------------

    # plot heads
    # idom_plt = np.where(idomain[0] == 1, np.nan, idomain[0])  # create a mask for the active domain
    # head = gwf.output.head().get_data()
    # for i in range(NLAY):
        
    #     pmv = fp.plot.PlotMapView(model=gwf, layer=i)

    #     pmv.plot_bc(name='ghb_pw', color='purple') # add 'chd' cells
    #     pmv.plot_bc(name='mbr', color='orange') # add 'wells' cells
    #     pmv.plot_bc(name='drn_r', color='blue') # add 'wells' cells
    #     # pmv.plot_bc(name='chd_conf_out', color='gold') # add 'chd' cells
    #     pmv.plot_bc(name='influx', color='green') # add 'influx' cells
    #     pmv.plot_bc(name='outflux', color='red') # add 'outflux' cells
        
    #     # pmv.plot_inactive(color='lightgray', alpha=0.5) # plot inactive cells
    #     # pmv.plot_array(idom_plt, cmap='gray', alpha=0.2) # plot head array for layer 0
    #     # pmv.plot_grid(colors='silver', lw=0.1); # add grid
    #     cs = pmv.contour_array(head[i], linewidths=1, colors='k') # contour plot of heads
    #     plt.clabel(cs, fmt='%1.1f'); # add contour labels with one decimal place

    #     # save to figures
    #     plt.savefig(Path(FIG_DIR, f'{MODEL_NAME}_heads{i}.png'), dpi=300, bbox_inches='tight') # save figure
    #     plt.close()

    # -------------------------------------------------------------
    # create a cross section grid & plot heads

    # cross_col = 15 #which row to show cross section
    # crossview = fp.plot.PlotCrossSection(model=gwf, line={'column': cross_col})
    # strtArray = crossview.plot_array(head, masked_values=[1e30], alpha=0.8) # plot the active domain
    # crossview.plot_grid(colors='black', lw=1); # add grid
    # crossview.plot_inactive(color='lightgray') # plot inactive cells
    # cb = plt.colorbar(strtArray, shrink=0.5) # add color bar
    # # strtArray = crossview.plot_array(head, masked_values=[1e30], alpha = 0.5) # plot the array of heads in cross section
    # plt.savefig(Path(FIG_DIR, f'{MODEL_NAME}_xsection_col{cross_col}.png'), dpi=300, bbox_inches='tight') # save figure
    # plt.close()  # close the figure to avoid memory issues
    TRUTHREL_DIR = os.path.abspath(TRUTH_DIR)
    # TEST OBS ------------------------------------------------------------
    os.chdir(MODEL_DIR)  # change directory to model directory
    sample_path = os.path.relpath(SAMPLES, MODEL_DIR)
    helpers.extract_heads_and_budget(model_name=f'{MODEL_NAME}')  # extract heads and budget from model output
    helpers.extract_spring_obs(
        gwf=None, 
        model_name=f'{MODEL_NAME}', 
        samples_path=f'{sample_path}',
        conf_h_path=f'{fn_out['ghb_conf_ts']['timeseries']['filename']}',
        aw_h_path=f'{fn_out['ghb_aw_ts']['timeseries']['filename']}',
        pw_h_path=f'{fn_out['ghb_pw_ts']['timeseries']['filename']}',
        spring_h_path=f'{fn_out['ghb_sp_ts']['timeseries']['filename']}',
        pw_q_path=f'{MODEL_NAME}.poukawa_flow_m3d.csv',
        )  # extract spring 

    # CREATE TRUTH -------------------------------------------------------
    # copy file to truth directory
    obs_results  = pd.read_csv(f'obs_results.csv')
    obs_results.loc[:, 'pk4head'] = 7.8  # set all to 7.8 initially
    mask = (obs_results['time'] <= 53) & (obs_results['time'] > 1)
    # pk4 head
    obs_results.loc[mask, 'pk4head'] = (
        obs_results.loc[mask, 'time']
        .apply(lambda x: pk4_h.loc[int(x), 'level_mRL'])
        .values
    )
    obs_results['pk4confdiff'] = obs_results['pk4head'] - obs_results['confhead']
    obs_results['pk4springdiff'] = obs_results['pk4head'] - obs_results['springhead']
    obs_results['awswgwdiff'] = 0 # known positive flux
    obs_results['pwswgwdiff'] = 0 # known positive flux
    
    obs_results.loc[:, 'awtot'] = aw_Q['sm_m3d'].iloc[0] * -1  # mean monthly 10/11
    obs_results.loc[mask, 'awtot'] = (
        obs_results.loc[mask, 'time']
        .apply(lambda x: aw_Q.loc[int(x), 'sm_m3d'] * -1)
        .values
    )
    
    # inflow from poukawa + local flow
    # local flow assume 0.0013 m/d per m2 of catchment area
    # poukawa area = 830m * 1m = 830m2
    # local flow (pw_cum) = 830m2 * 0.0013m/d = 1.079 m3/d
    # total inflow = poukawa + local flow
    obs_results.loc[:, 'pwtot'] = pw_Q['sm_m3d'].iloc[0] * -1  # mean monthly 10/11
    obs_results.loc[mask, 'pwtot'] = (
        obs_results.loc[mask, 'time']
        .apply(lambda x: pw_Q.loc[int(x), 'sm_m3d'] * -1)
        .values
    )

    mask = (obs_results['pk4springdiff'] < 0)
    obs_results['springq'] = -3
    obs_results.loc[mask, 'springq'] = 0  # no spring flow if head diff is negative
    obs_results['springcum'] = -150
    obs_results.loc[mask, 'springcum'] = 0  # no spring flow if head diff is negative
    obs_results['pwcum'] = obs_results['pwtot'] * 0.0002
    obs_results['pwtot'] = obs_results['pwtot'] + obs_results['pwcum']
    obs_results['awcum'] = obs_results['awtot'] - obs_results['pwtot'] - obs_results['springcum']
    
    # save to truth directory
    obs_results.to_csv(Path(TRUTHREL_DIR, 'obs_results_truth.csv'), index=False)

    

if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations