import os
import glob
import re
import stat
import shutil
import pyemu
import numpy as np
import flopy as fp
import pandas as pd
import matplotlib.pyplot as plt

from setup import *
from utils_pest import *
from helpers import *

def main():

    if os.path.exists(PEST_DIR):
        shutil.rmtree(PEST_DIR)
    shutil.copytree(MODEL_DIR, PEST_DIR)
    
    # copy all the contents of bin into model directory
    if os.path.exists(BIN_DIR):
        if os.name == 'nt':  # if on Windows, copy files
            os_bin = os.path.join(BIN_DIR, 'windows')

        elif os.name == 'posix':  # if on Linux or MacOS, copy files
            os_bin = os.path.join(BIN_DIR, 'linux')
        else:
            raise ValueError(f'Unsupported OS: {os.name}. Please check the BIN_DIR path.')

        files = os.listdir(os_bin)
        for f in files:
            if os.path.exists(os.path.join(PEST_DIR, f)):
                file_path = os.path.join(PEST_DIR, f)
                os.remove(file_path)
                shutil.copy2(os.path.join(os_bin, f), file_path)
                current_permissions = os.stat(file_path).st_mode
                os.chmod(file_path, current_permissions | stat.S_IEXEC | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
                print(f"Made executable: {file_path}")
                        
    else:
        print(f'Bin directory {BIN_DIR} does not exist. Please check the path.')
# -----------------------------------------------------------------

    # load simulation
    sim = fp.mf6.MFSimulation.load(sim_ws=PEST_DIR)
    # load flow model
    gwf = sim.get_model()

    # run the model once to make sure it works
    # pyemu.os_utils.run("mf6", cwd=PEST_DIR)
    # run modpath7
    # pyemu.os_utils.run(f'mp7 {MODEL_NAME}.mpsim', cwd=PEST_DIR)

    # -----------------------------------------------------------------

    sr = pyemu.helpers.SpatialReference.from_namfile(
            os.path.join(PEST_DIR, f"{MODEL_NAME}.nam"),
            delr=gwf.dis.delr.array, delc=gwf.dis.delc.array)

    # instantiate PstFrom
    pf = pyemu.utils.PstFrom(original_d=PEST_DIR, # where the model is stored
                            new_d=HPTEMP_DIR, # the PEST template folder
                            remove_existing=True, # ensures a clean start
                            longnames=False, # set False if using PEST/PEST_HP
                            spatial_reference=sr, #the spatial reference we generated earlier
                            zero_based=False, # does the MODEL use zero based indices? For example, MODFLOW does NOT
                            # start_datetime=start_datetime, # required when specifying temporal correlation between parameters
                            echo=False) # to stop PstFrom from writting lots of infromation to the notebook; experiment by setting it as True to see the difference; usefull for troubleshooting


    # PARAMETERIZATION --------------------------------------------------
    # exponential variogram for spatially varying parameters
    v_space = pyemu.geostats.ExpVario(
        name='main_gs', contribution=1.0, a=1000, anisotropy=1.0, bearing=0.0)

    v_fine = pyemu.geostats.ExpVario(
        name='fine_gs', contribution=1.0, a=100, anisotropy=1.0, bearing=0.0)

    # geostatistical structure for spatially varying parameters
    fine_gs = pyemu.geostats.GeoStruct(variograms=v_fine, transform='log')
    grid_gs = pyemu.geostats.GeoStruct(variograms=v_space, transform='log')

    # plot the gs if you like:
    # _ = grid_gs.plot()

    # get the IDOMAIN array
    ib = gwf.dis.idomain.array

    # setup pilot points

    # # set up pst file
    # K ----------------------------------------------
    # shallow aquifer
    define_mult_array(
        pf, HPTEMP_DIR,
        tag=f'{MODEL_NAME}.npf_k_layer',
        sr=sr,
        ib=ib,
        grid_gs=grid_gs,
        fine_gs=fine_gs,
        lb=0.01, ub=1e2, 
        # ulb=1e-4, uub=1e3,
        add_coarse=True,
        pp_space=5,
        lays=np.arange(NLAY).tolist()
        )
    
    define_mult_array(
        pf, HPTEMP_DIR,
        tag=f'{MODEL_NAME}.sto_ss_layer',
        sr=sr,
        ib=ib,
        grid_gs=grid_gs,
        fine_gs=fine_gs,
        lb=1e-2, ub=1e2, 
        # ulb=1e-4, uub=1e3,
        add_coarse=True,
        pp_space=5,
        lays=np.arange(NLAY).tolist()
        )

    # RECHARGE ------------------------------------------------------
    rcha(
        pf, HPTEMP_DIR, ib,
        tag=f'{MODEL_NAME}.rcha_recharge',
        grid_gs=grid_gs,
        fine_gs=fine_gs,
        constant_gs=None,
        rch_bounds=[1e-2, 1e2],
        pp_space=5,
        )

    ghb(pf, HPTEMP_DIR,
        name='ghb-aw',
        tag=f'{MODEL_NAME}.ghbaw_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        cond_bounds=[1e-4, 1e4],
        # cond_ultbounds=[1e-6, 1e6],
        head_bounds=[-2, 2],
        # head_ultbounds=[5, 20]
        )
    
    ghb(pf, HPTEMP_DIR,
        name='ghb-pw',
        tag=f'{MODEL_NAME}.ghb_pw_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        cond_bounds=[1e-4, 1e4],
        # cond_ultbounds=[1e-6, 1e6],
        head_bounds=[-2, 2],
        # head_ultbounds=[5, 20]
        )
    
    ghb(pf, HPTEMP_DIR,
        name='ghb-spring',
        tag=f'{MODEL_NAME}.ghbspr_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        cond_bounds=[1e-4, 1e4],
        # cond_ultbounds=[1e-6, 1e6],
        head_bounds=[-2, 2],
        # head_ultbounds=[5, 15]
        )
    
    ghb(pf, HPTEMP_DIR,
        name='ghb-conf',
        tag=f'{MODEL_NAME}.ghb_conf_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        cond_bounds=[1e-4, 1e4],
        # cond_ultbounds=[1e-6, 1e6],
        head_bounds=[-2, 2],
        # head_ultbounds=[5, 20]
        )
    
    # add ts data
    # add_ts_parameters(pf, 'spring', HPTEMP_DIR, f'{MODEL_NAME}.spring_stage.csv', -1, 1)
    # add_ts_parameters(pf, 'pw', HPTEMP_DIR, f'{MODEL_NAME}.ghb_pw_heads.csv', -1, 1)
    # add_ts_parameters(pf, 'conf', HPTEMP_DIR, f'{MODEL_NAME}.ghb_conf_heads.csv', -1, 1)
    # add_ts_parameters(pf, 'aw', HPTEMP_DIR, f'{MODEL_NAME}.awanui_stage.csv', -1, 1)
    
    
    # OBSERVATIONS ------------------------------------------------------
    # budget = pd.read_csv(os.path.join(TRUTH_DIR, "incremental_budget.csv"))
    # index_cols = [list(budget.columns.values)[0]]
    # use_cols = list(budget.columns.values)[1:]
    # obsgp='budget'
    # pf.add_observations(
    #     'incremental_budget.csv',
    #     index_cols=index_cols,
    #     use_cols=use_cols,
    #     prefix=obsgp,
    #     obsgp=obsgp)
    
    # pst = pf.build_pst()

    # budget.set_index('time', inplace=True)
    # for col in use_cols:
    #     for time in budget.index:
    #         t = int(time)
    #         pst.observation_data.loc[
    #             f'oname:{obsgp}_otype:lst_usecol:{col}_tim:{t}','obsval'] = budget[col].iloc[t]
            
    #         pst.observation_data.loc[
    #             f'oname:{obsgp}_otype:lst_usecol:{col}_totim:{t}','weight'] = budget[col].iloc[2]


    # add stress period head/fluxes observations ------------------------------
    # std_obs_results  = pd.read_csv(os.path.join(TRUTH_DIR, f'obs_results_std.csv'))
    index_cols = ['time']
    obs_use_cols = ['pk4', 'pk4-spr-diff', 'pk4-conf-diff']
    pf.add_observations(
        'output.sample_heads.csv', # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=obs_use_cols, # skip the index column
        prefix='ts-heads',
        obsgp='ts-heads',
    )

    # std_obs_results  = pd.read_csv(os.path.join(TRUTH_DIR, f'obs_results_std.csv'))
    index_cols = ['time']
    obs_use_cols = ['flux']
    pf.add_observations(
        'output.spring_fluxes.csv', # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=obs_use_cols, # skip the index column
        prefix='ts-flux',
        obsgp='ts-flux',
    )

    # add budget
    index_cols = ['kper', 'kstp', 'k', 'i', 'j']
    use_cols = ['AWq']
    pf.add_observations(
        "output.GHB_AW_fluxes.csv", # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols, # skip the index column
        prefix='AWq',
        obsgp='AWq',
    )

    # add confined flux
    index_cols = ['kper', 'kstp', 'k', 'i', 'j']
    use_cols = ['CONFq']
    pf.add_observations(
        "output.GHB_CONF_fluxes.csv", # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols, # skip the index column
        prefix='CONFq',
        obsgp='CONFq',
    )

    # add budget
    budget = pd.read_csv(os.path.join(HPTEMP_DIR, "output.budget.csv"))
    index_cols = ['kper']
    use_cols = list(budget.iloc[:, 1:].columns)
    pf.add_observations(
        "output.budget.csv", # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols, # skip the index column
        prefix='budget',
        obsgp='budget',
    )

    # add head arrays
    pattern = os.path.join(HPTEMP_DIR, "output.heads_lyr*.dat")
    head_files = glob.glob(pattern)

    # Extract layer, kper, kstp from filenames
    file_info = []
    for file_path in head_files:
        filename = os.path.basename(file_path)
        # Use regex to extract numbers
        match = re.search(r'output\.heads_lyr(\d+)_kper(\d+)_kstp(\d+)\.dat', filename)
        if match:
            lyr = int(match.group(1))
            kper = int(match.group(2))
            kstp = int(match.group(3))
            file_info.append({
                'file_path': file_path,
                'filename': filename,
                'layer': lyr,
                'kper': kper,
                'kstp': kstp
            })
    for info in file_info:
        file_path = info['file_path']
        lyr = info['layer']
        kper = info['kper']
        kstp = info['kstp']
        
        # Add as observation to pyemu
        pf.add_observations(
            info['filename'],
            prefix=f'h-lyr{lyr}-kper{kper}-kstp{kstp}',
            obsgp=f'head-arr',
            zone_array=ib[0,:,:], # assuming layer 1 for zone array
        )


    # FORWARD RUN SCRIPT --------------------------------------------------
    pf.mod_sys_cmds.append("mf6") #do this only once

    sample_rel = os.path.relpath(SAMPLES, HPTEMP_DIR)
    # post-processing to get observations
    pf.add_py_function(
        f"{SCRIPTS_DIR}/helpers.py",
        f"extract_budget(model_name=f'{MODEL_NAME}')", is_pre_cmd=False)
    pf.add_py_function(
        f"{SCRIPTS_DIR}/helpers.py",
        f"extract_model_heads(model_name=f'{MODEL_NAME}', sample_path='{sample_rel}')", is_pre_cmd=False)
    pf.add_py_function(
        f"{SCRIPTS_DIR}/helpers.py",
        f"extract_ghb_fluxes(model_name=f'{MODEL_NAME}', sample_path='{sample_rel}')", is_pre_cmd=False)

    pst = pf.build_pst()

    ###############################################################################
    ts_heads = pd.read_csv(os.path.join(TRUTH_DIR, 'output.sample_heads.truth.csv'), index_col=0)
    AWq = pd.read_csv(os.path.join(TRUTH_DIR, "output.GHB_AW_fluxes.truth.csv"), index_col=[0, 1, 2, 3, 4])
    heads = np.loadtxt(
        os.path.join(TRUTH_DIR, "output.heads.truth.dat"))
    heads_std = np.loadtxt(
        os.path.join(TRUTH_DIR, "output.heads.std.dat"))
    heads_weight = np.loadtxt(
        os.path.join(TRUTH_DIR, "output.heads.weight.dat"))

    # adjust obgnme to main groups (for some reason not working above)
    pst.observation_data['obgnme'] = pst.observation_data.apply(
        lambda x: x['obgnme'].split('_')[0].split(':')[-1], axis=1)
    pst.observation_data['standard_deviation'] = 0
    pst.observation_data['weight'] = 0

    # add truth values for ts-heads
    for _, row in pst.observation_data.iterrows():
        obgnme = row['obgnme']

        if obgnme == 'ts-heads':
            time = int(row['time'])
            col = row['usecol']
            
            try:
                pst.observation_data.at[row.name, 'obsval'] = ts_heads.loc[time, col]
                pst.observation_data.at[row.name, 'standard_deviation'] = ts_heads.loc[time, 'std']
                pst.observation_data.at[row.name, 'weight'] = ts_heads.loc[time, 'weight']

            except:
                continue
        
        elif obgnme == 'AWq':
            kper = int(row['kper'])
            kstp = int(row['kstp'])
            k = int(row['k'])
            i = int(row['i'])
            j = int(row['j'])
            
            pst.observation_data.at[row.name, 'obsval'] = AWq.loc[(kper, kstp, k, i, j), 'AWq']
            pst.observation_data.at[row.name, 'standard_deviation'] = AWq.loc[(kper, kstp, k, i, j), 'std']
            pst.observation_data.at[row.name, 'weight'] = AWq.loc[(kper, kstp, k, i, j), 'weight']
        
        elif obgnme == 'head-arr':
            i = int(row['i'])
            j = int(row['j'])
            
            pst.observation_data.at[row.name, 'obsval'] = heads[i, j]
            pst.observation_data.at[row.name, 'standard_deviation'] = heads_std[i, j]
            pst.observation_data.at[row.name, 'weight'] = heads_weight[i, j]
            pst.observation_data.at[row.name, 'obgnme'] = 'less_' + row['obgnme']
    
    # create phi factor file
    phi_dict = {
        'head-arr': 20,
        'ts-heads': 50,
        'AWq': 5,
    }
    # to dataframe with no column names or index
    df = pd.DataFrame(list(phi_dict.items()))
    df.to_csv(os.path.join(HPTEMP_DIR, 'phi_factors.csv'), index=False, header=False)
    
    
    ## ADD FORECASTS ------------------------------------------------------
    forecast_obgnme = ['ts-flux', 'budget']
    
    mask = pst.observation_data.loc[:,'obgnme'].isin(forecast_obgnme)
    forecasts = pst.observation_data[mask].index.tolist()
    
    pst.pestpp_options['forecasts'] = forecasts


    # WRITE PEST -------------------------------------------------------
    print("Writing PEST template file...")
    pst.control_data.noptmax = 0 # just run parameter values in the file
    pst.observation_data = pst.observation_data.fillna(0)
    pst_file = f'{MODEL_NAME}.pst'
    pst.write(os.path.join(HPTEMP_DIR, pst_file), version=2)


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations