import matplotlib.pyplot as plt
import os
from pathlib import Path
import shutil
import numpy as np
import flopy as fp
import pandas as pd

import helpers as helpers

from utils import *
from setup import *
from pkg_setup import *

def main():

    # print working directory
    print(f'Working directory: {os.getcwd()}')
    print(f'Model directory: {os.path.abspath(MODEL_DIR)}')

    # create directories if they do not exist
    for d in [SPATIAL_DIR, FIG_DIR, POST_DIR]:
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

    # save the setup.py file as a text file in MODEL_DIR
    setup_file_path = Path(__file__).parent / 'setup.py'
    destination_path = Path(MODEL_DIR) / 'model_setup.txt'
    shutil.copy2(setup_file_path, destination_path)
    print(f'Setup file copied to: {destination_path}')

    fn_out = {}
    # TDIS ####################################################################
    start = pd.to_datetime(START_DATE)
    end = pd.to_datetime(END_DATE)
    PERLEN = (end - start).days
    NSTEPS = PERLEN  # 1 day time step

    print(f'Simulation period: {PERLEN} days with {NSTEPS} time steps.')
    
    # dis
    grid, idomain, top, nrow, ncol, delr, delc, fn_out, model_thickness = dis_setup(fn_out)
    elev_constraint = grid.array_from_raster(TOP)
    grid_gpd = grid.cell_geodataframe()
    grid_gpd['idomain'] = idomain[0].flatten()
    grid_gpd.to_file(Path(SPATIAL_DIR, f'{MODEL_NAME}_grid.shp'), driver='ESRI Shapefile')
    # ghbs
    active_domain = np.where(idomain > 0, 1, 0)
    fn_out, aw_ts, wetland_WL = ghb_aw_setup(grid, active_domain, start, end, fn_out)
    fn_out = ghb_spring_setup(grid, active_domain, start, end, aw_ts, wetland_WL, fn_out)
    fn_out, ghb_pw_kper0 = ghb_pw_setup(grid, active_domain, start, end, fn_out)
    fn_out = ghb_conf_setup(grid, active_domain, start, end, model_thickness, fn_out)
    # wel
    # fn_out, mbr_df = wel_mbr_setup(grid, idomain, fn_out, ghb_pw_kper0)
    # fn_out = wel_inout_setup(grid, idomain, mbr_df, fn_out)
    # npf
    fn_out = npf_setup(grid, active_domain, fn_out)
    # sto
    fn_out = sto_ss_setup(active_domain, fn_out)
    # rch
    fn_out = rch_setup(active_domain, fn_out)
    #obs - save file
    pk4_h = pk4_ts(start, end)
    aw_Q = awanui_ts(start, end)
    pw_Q = poukawa_ts(start, end)

    # copy PP_SHP into model directory; need to copy all files associated with shapefile
    pp_shp_path = Path(PP_SHP)
    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg']:
        src_file = pp_shp_path.with_suffix(ext)
        if src_file.exists():
            shutil.copy2(src_file, MODEL_DIR)

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
        start_date_time=start.strftime('%Y-%m-%d'), #YYYY-MM-DDThh:mmTZD
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
        ss = fn_out['sto_ss'],
        # sy=SS_PRIOR['sy'], 
        # ss=SS_PRIOR['ss'],
        # iconvert=SS_PRIOR['iconvert'],
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

    rch = fp.mf6.ModflowGwfrcha(
        model=gwf,
        recharge=fn_out['rch'], # recharge file names for each layer
        save_flows=True, # save flows for this package
        pname='rch' # package name
    )
    
    ghb_aw = fp.mf6.ModflowGwfghb(
        print_input=True,
        model=gwf, # add riv package to model gwf (created in previous code cell)
        timeseries=fn_out['ghb_aw_ts'],
        stress_period_data=fn_out['ghb_aw'],
        pname='ghb_aw', # package name
        save_flows=True,
        )
    # ghb_aw.ts.initialize(fn_out['ghb_aw_ts'])
    
    ghb_conf = fp.mf6.ModflowGwfghb(
        print_input=True,
        # print_flows=True,
        save_flows=True, # save flows for this package 
        model=gwf,
        timeseries=fn_out['ghb_conf_ts'],
        stress_period_data=fn_out['ghb_conf'],
        pname='ghb_conf', # package name
        )
    # ghb_conf.ts.initialize(fn_out['ghb_conf_ts'])

    ghb_pw = fp.mf6.ModflowGwfghb(
        print_input=True,
        model=gwf,
        timeseries=fn_out['ghb_pw_ts'],
        stress_period_data=fn_out['ghb_pw'],
        pname='ghb_pw', # package name
        save_flows=True, # save flows for this package 
        )
    
    ghb_sp = fp.mf6.ModflowGwfghb(
        print_input=True,
        model=gwf,
        timeseries=fn_out['ghb_sp_ts'],
        stress_period_data=fn_out['ghb_sp'],
        pname='ghb_sp', # package name
        save_flows=True, # save flows for this package 
        )

    # wel_mbr = fp.mf6.ModflowGwfwel(
    #     model=gwf,
    #     stress_period_data=fn_out['wel_mbr'],
    #     pname='wel_mbr', # package name
    #     save_flows=True,
    #     )
    
    # influx = fp.mf6.ModflowGwfwel(
    #     model=gwf,
    #     stress_period_data=fn_out['wel_influx'],
    #     pname='influx', # package name
    #     save_flows=True,
    #     )
    # outflux = fp.mf6.ModflowGwfwel(
    #     model=gwf,
    #     stress_period_data=fn_out['wel_outflux'],
    #     pname='outflux', # package name
    #     save_flows=True,
    #     )

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

    
    TRUTHREL_DIR = os.path.abspath(TRUTH_DIR)
    # TEST OBS ------------------------------------------------------------
    os.chdir(MODEL_DIR)  # change directory to model directory
    sample_path = os.path.relpath(SAMPLES, MODEL_DIR)
    
    budget = helpers.extract_budget(model_name=f'{MODEL_NAME}')  # extract heads and budget from model output    
    samples = helpers.extract_model_heads(model_name=f'{MODEL_NAME}', sample_path=sample_path)
    ghb_dfs = helpers.extract_ghb_fluxes(model_name=f'{MODEL_NAME}', sample_path=sample_path)

    # CREATE TRUTH -------------------------------------------------------
    # tdis_data = sim.tdis.perioddata.get_data()
    
    # awanui flux truth
    helpers.create_awghb_truth(ghb_dfs, TRUTHREL_DIR)
    helpers.create_sprghb_truth(ghb_dfs, TRUTHREL_DIR)
    helpers.create_confghb_truth(ghb_dfs, TRUTHREL_DIR)
    # top model constraints (head can't be greater than 1m above top)
    helpers.create_head_truth(elev_constraint, TRUTHREL_DIR)
    # sample location truth
    helpers.samples_truth(gwf, TRUTHREL_DIR)
    helpers.create_budget_truth(budget, TRUTHREL_DIR)


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations