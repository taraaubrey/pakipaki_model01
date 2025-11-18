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

    start_datetime = pd.to_datetime(START_DATE)

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
                            new_d=TEMP_DIR, # the PEST template folder
                            remove_existing=True, # ensures a clean start
                            longnames=False, # set False if using PEST/PEST_HP
                            spatial_reference=sr, #the spatial reference we generated earlier
                            zero_based=False, # does the MODEL use zero based indices? For example, MODFLOW does NOT
                            # start_datetime=start_datetime, # required when specifying temporal correlation between parameters
                            echo=False) # to stop PstFrom from writting lots of infromation to the notebook; experiment by setting it as True to see the difference; usefull for troubleshooting


    # PARAMETERIZATION --------------------------------------------------
    # exponential variogram for spatially varying parameters
    v_coarse = pyemu.geostats.ExpVario(
        name='main_gs', contribution=1.0, a=1000, anisotropy=1.0, bearing=0.0)
    v_fine = pyemu.geostats.ExpVario(
        name='fine_gs', contribution=1.0, a=100, anisotropy=1.0, bearing=0.0)
    
    vh_coarse = pyemu.geostats.GauVario(
        name='head_gs', contribution=1.0, a=1000, anisotropy=1.0, bearing=0.0)
    vh_fine = pyemu.geostats.GauVario(
        name='head_gs', contribution=1.0, a=100, anisotropy=1.0, bearing=0.0)

    v_time = pyemu.geostats.ExpVario(
        name='time_gs', contribution=1.0, a=14, anisotropy=1.0, bearing=0.0)

    # geostatistical structure for spatially varying parameters
    fine_gs = pyemu.geostats.GeoStruct(variograms=v_fine, transform='log')
    coarse_gs = pyemu.geostats.GeoStruct(variograms=v_coarse, transform='log')
    
    h_fine_gs = pyemu.geostats.GeoStruct(variograms=vh_fine, transform='none')
    h_coarse_gs = pyemu.geostats.GeoStruct(variograms=vh_coarse, transform='none')
    
    time_gs = pyemu.geostats.GeoStruct(variograms=v_time, transform='none')

    # plot the gs if you like:
    # _ = grid_gs.plot()

    # get the IDOMAIN array
    ib = gwf.dis.idomain.array


    # TIMESERIES data ------------------------------------
    aw_df = pd.read_csv(os.path.join(TEMP_DIR, f'{MODEL_NAME}.ghb_aw_head_names.csv'), header=None)
    sp_df = pd.read_csv(os.path.join(TEMP_DIR, f'{MODEL_NAME}.ghb_spring_head_names.csv'), header=None)
    pw_df = pd.read_csv(os.path.join(TEMP_DIR, f'{MODEL_NAME}.ghb_pw_head_names.csv'), header=None)

    dts = pd.to_datetime(start_datetime) + pd.to_timedelta(np.cumsum(sim.tdis.perioddata.array["perlen"]),unit='d')
    dts = dts.to_list()

    # add ts data
    # add_ts_parameters(
    #     pf, TEMP_DIR,
    #     datetimes=dts,
    #     name='ghb-ts-aw',
    #     f=f'{MODEL_NAME}.ghb_aw_heads.csv',
    #     bounds=[-1e-1, 1e-1],
    #     parnames=list(aw_df[0]),
    #     use_cols=np.arange(1, aw_df.shape[0]+1).tolist(),
    #     time_gs={'tgs': time_gs},
    #     )
    # add_ts_parameters(
    #     pf, TEMP_DIR,
    #     datetimes=dts,
    #     name='ghb-ts-pw',
    #     f=f'{MODEL_NAME}.ghb_pw_heads.csv',
    #     bounds=[-1e-1, 1e-1],
    #     parnames=list(pw_df[0]),
    #     use_cols=np.arange(1, pw_df.shape[0]+1).tolist(),
    #     time_gs={'tgs': time_gs},
    #     )
    # add_ts_parameters(
    #     pf, TEMP_DIR,
    #     datetimes=dts,
    #     name='ghb-ts-spring',
    #     f=f'{MODEL_NAME}.ghb_spring_heads.csv',
    #     bounds=[-1e-1, 1e-1],
    #     parnames=list(sp_df[0]),
    #     use_cols=np.arange(1, sp_df.shape[0]+1).tolist(),
    #     time_gs={'tgs': time_gs},
    #     )
    # add_ts_parameters(
    #     pf, TEMP_DIR,
    #     datetimes=dts,
    #     name='ghb-ts-conf',
    #     f=f'{MODEL_NAME}.ghb_conf_heads.csv',
    #     bounds=[-1e-1, 1e-1],
    #     parnames=['head_ts'],
    #     use_cols=[1],
    #     time_gs={'tgs': time_gs},
    #     )
    # # set up pst file
    # K ----------------------------------------------
    # shallow aquifer
    K_ins = [
        {
            'type': 'grid',
            'gs': fine_gs,
            'zone_array': np.where(ib[0] == 2, 1, 0),
            'name_suffix': '-gr'
        },
        {
            'type': 'pilotpoints',
            'gs': coarse_gs,
            'pp_space': 'pilot_points.shp',
            'zone_array': np.where(ib[0] > 0, 1, 0),
            'name_suffix': '-pp'
        },
        {
            'type': 'constant',
            'zone_array': np.where(ib[0] > 0, 1, 0),
            'name_suffix': '-cn'
        },
    ]
    
    define_mult_array(
        pf, TEMP_DIR,
        tag=f'{MODEL_NAME}.npf_k_layer',
        ib=ib,
        p_ins = K_ins,
        lb=KH_PRIOR['lb'], ub=KH_PRIOR['ub'],
        ulb=KH_PRIOR['ulb'], uub=KH_PRIOR['uub'], 
        lays=np.arange(NLAY).tolist()
        )
    # S ----------------------------------------------
    define_mult_array(
        pf, TEMP_DIR,
        tag=f'{MODEL_NAME}.sto_ss_layer',
        ib=ib,
        p_ins = K_ins,
        lb=SS_PRIOR['lb'], ub=SS_PRIOR['ub'],
        ulb=SS_PRIOR['ulb'], uub=SS_PRIOR['uub'],
        lays=np.arange(NLAY).tolist()
        )

    # RECHARGE ------------------------------------------------------
    # define_mult_array(
    #     pf, TEMP_DIR,
    #     tag=f'{MODEL_NAME}.rcha_recharge',
    #     ib=ib,
    #     p_ins = K_ins,
    #     lb=RCH['lb'], ub=RCH['ub'],
    #     ulb=RCH['ulb'], uub=RCH['uub'],
    #     lays=[0, 2]
    #     )
    
    # GHBS ------------------------------------------------------
    ghb_cond(pf, TEMP_DIR,
        name='ghb-aw',
        tag=f'{MODEL_NAME}.ghbaw_stress_period_data',
        grid_gs={'gr': fine_gs},
        cond_bounds=[GHB_SW['cond_lb'], GHB_SW['cond_ub']],
        ult_bounds=[GHB_SW['cond_ulb'], GHB_SW['cond_uub']],
        )
    ghb_cond(pf, TEMP_DIR,
        name='ghb-pw',
        tag=f'{MODEL_NAME}.ghb_pw_stress_period_data',
        grid_gs={'gr': fine_gs},
        cond_bounds=[GHB_SW['cond_lb'], GHB_SW['cond_ub']],
        ult_bounds=[GHB_SW['cond_ulb'], GHB_SW['cond_uub']],
        )    
    ghb_cond(pf, TEMP_DIR,
        name='ghb-spring',
        tag=f'{MODEL_NAME}.ghbspr_stress_period_data',
        grid_gs={'gr': fine_gs},
        cond_bounds=[GHB_SW['cond_lb'], GHB_SW['cond_ub']],
        ult_bounds=[GHB_SW['cond_ulb'], GHB_SW['cond_uub']],
        )
    ghb_cond(pf, TEMP_DIR,
        name='ghb-conf',
        tag=f'{MODEL_NAME}.ghb_conf_stress_period_data',
        grid_gs={'gr': coarse_gs},
        cond_bounds=[GHB_CONF['cond_lb'], GHB_CONF['cond_ub']],
        ult_bounds=[GHB_CONF['cond_ulb'], GHB_CONF['cond_uub']],
        )
    
    ghb_heads(pf, TEMP_DIR,
        name='ghb-aw',
        tag=f'{MODEL_NAME}.ghbaw_stress_period_data',
        head_gs={'hg': h_fine_gs},
        head_bounds=[GHB_SW['head_lb'], GHB_SW['head_ub']],
        ult_bounds=[GHB_SW['head_ulb'], GHB_SW['head_uub']],
        )
    ghb_heads(pf, TEMP_DIR,
        name='ghb-pw',
        tag=f'{MODEL_NAME}.ghb_pw_stress_period_data',
        head_gs={'hg': h_fine_gs},
        head_bounds=[GHB_SW['head_lb'], GHB_SW['head_ub']],
        ult_bounds=[GHB_SW['head_ulb'], GHB_SW['head_uub']],
        )
    ghb_heads(pf, TEMP_DIR,
        name='ghb-spring',
        tag=f'{MODEL_NAME}.ghbspr_stress_period_data',
        head_gs={'hg': h_fine_gs},
        head_bounds=[GHB_SW['head_lb'], GHB_SW['head_ub']],
        ult_bounds=[GHB_SW['head_ulb'], GHB_SW['head_uub']],
        )
    ghb_heads(pf, TEMP_DIR,
        name='ghb-conf',
        tag=f'{MODEL_NAME}.ghb_conf_stress_period_data',
        head_gs={'hg': h_coarse_gs},
        head_bounds=[GHB_CONF['head_lb'], GHB_CONF['head_ub']],
        ult_bounds=[GHB_CONF['head_ulb'], GHB_CONF['head_uub']]
        )

    print('*' * 40)
    par_info = [(df['pargp'].unique(), df.shape[0]) for df in pf.par_dfs]
    total_pars = sum([npar for grp, npar in par_info])
    for grp, npar in par_info:
        print(f'Parameter group: {grp[0]}, number of parameters: {npar}')
    print(f'Total number of parameters: {total_pars}')
    print('*' * 40)

    # OBSERVATIONS ------------------------------------------------------
    print('SETTING OBSERVATIONS...')

    # add stress period head/fluxes observations
    index_cols = ['kper', 'kstp', 'time']
    obs_use_cols = ['pk4', 'pk4-spr-diff', 'pk4-conf-diff']
    pf.add_observations(
        'output.sample_heads.csv', # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=obs_use_cols, # skip the index column
        prefix='ts-heads',
        obsgp='ts-heads',
    )

    index_cols = ['kper', 'kstp', 'time']
    obs_use_cols = ['flux']
    pf.add_observations(
        'output.spring_fluxes.csv', # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=obs_use_cols, # skip the index column
        prefix='ts-sprflux',
        obsgp='ts-sprflux',
    )

    # add budget
    index_cols = ['kper', 'kstp', 'time', 'i', 'j']
    use_cols = ['AWq']
    pf.add_observations(
        "output.GHB_AW_fluxes.csv", # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols, # skip the index column
        prefix='arr-awq',
        obsgp='arr-awq',
    )

    # add budget
    index_cols = ['kper', 'kstp', 'time', 'i', 'j']
    use_cols = ['SPq']
    pf.add_observations(
        "output.GHB_SP_fluxes.csv", # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols, # skip the index column
        prefix='arr-spq',
        obsgp='arr-spq',
    )

    # add confined flux
    index_cols = ['kper', 'kstp', 'time', 'i', 'j']
    use_cols = ['CONFq']
    pf.add_observations(
        "output.GHB_CONF_fluxes.csv", # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols, # skip the index column
        prefix='arr-confq',
        obsgp='arr-confq',
    )

    # add recession
    index_cols = ['kper']
    use_cols=['pk4', 'pk4-spr-diff', 'pk4-conf-diff', 'fliptime']
    pf.add_observations(
        "output.sample_recession_rates.csv", # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols,
        prefix='recession',
        obsgp='recession',
    )

    # add budget
    index_cols = ['kper']
    use_cols = ['awanui','confined','poukawa','spring','inout','percentdiscrepancy','recharge','stoss','total','sw' ,'inflow']
    pf.add_observations(
        "output.budget.csv", # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols, # skip the index column
        prefix='budget',
        obsgp='budget',
    )

    # add head arrays
    pattern = os.path.join(TEMP_DIR, "output.heads_*.dat")
    head_files = glob.glob(pattern)

    # Extract layer, kper, kstp from filenames
    file_info = []
    for file_path in head_files:
        filename = os.path.basename(file_path)
        # Use regex to extract numbers
        match = re.search(r'output.heads_kper(\d+)_kstp(\d+)_time(\d+).(\d+).dat', filename)
        if match:
            kper = int(match.group(1))
            kstp = int(match.group(2))
            time = int(match.group(3))
            file_info.append({
                'file_path': file_path,
                'filename': filename,
                'kper': kper,
                'kstp': kstp,
                'time': time
            })
    for info in file_info:
        file_path = info['file_path']
        time = info['time']
        kper = info['kper']
        kstp = info['kstp']

        # Add as observation to pyemu
        pf.add_observations(
            info['filename'],
            prefix=f'arr-h_kper:{kper}_kstp:{kstp}_time:{time}',
            obsgp=f'arr-h',
            zone_array=np.where(ib[0] > 0, 1, 0),
        )

    # FORWARD RUN SCRIPT --------------------------------------------------
    pf.mod_sys_cmds.append("mf6") #do this only once

    sample_rel = os.path.relpath(SAMPLES, TEMP_DIR)
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

    # # Add preferred-value regularization
    # pyemu.helpers.zero_order_tikhonov(
    #     pst, parbounds=True, reset=False
    # )

    ###############################################################################
    print('ADDING OBSERVATION VALUES AND WEIGHTS FROM TRUTH...')
    # load truth data
    ts_heads = pd.read_csv(
        os.path.join(
            TRUTH_DIR, 
            'output.sample_heads.truth.csv'), 
        index_col=0)
    AWq = pd.read_csv(
        os.path.join(
            TRUTH_DIR, 
            "output.GHB_AW_fluxes.truth.csv"), 
        index_col=['kper', 'kstp', 'i', 'j'])
    SPq = pd.read_csv(
        os.path.join(
            TRUTH_DIR, 
            "output.GHB_SP_fluxes.truth.csv"), 
        index_col=['kper', 'kstp', 'i', 'j'])
    CONFq = pd.read_csv(
        os.path.join(
            TRUTH_DIR, 
            "output.GHB_CONF_fluxes.truth.csv"), 
        index_col=['kper', 'kstp', 'i', 'j'])
    recession_truth = pd.read_csv(
        os.path.join(
            TRUTH_DIR, 
            'output.sample_recession_rates.truth.csv'), 
        index_col='kper')
    budget = pd.read_csv(
        os.path.join(
            TRUTH_DIR, 
            "output.budget.truth.csv"), 
        index_col='kper')

    heads = np.loadtxt(
        os.path.join(
            TRUTH_DIR, 
            "output.heads.truth.dat"))
    heads_std = np.loadtxt(
        os.path.join(
            TRUTH_DIR, 
            "output.heads.std.dat"))
    heads_weight = np.loadtxt(
        os.path.join(
            TRUTH_DIR, 
            "output.heads.weight.dat"))
    
    get_oname = lambda x: x['obgnme'].split('_')[0].split(':')[-1]

    # adjust obgnme to main groups (for some reason not working above)
    # pst.observation_data['obgnme'] = pst.observation_data.apply(get_oname, axis=1)
    pst.observation_data['standard_deviation'] = 0.
    pst.observation_data['weight'] = 0.

    print(f'Number of non-zero observations before adding weighting (check): {pst.nnz_obs}')

    # add truth values
    for _, row in pst.observation_data.iterrows():
        oname = row['oname']

        if oname == 'ts-heads':
            time = int(row['time'])
            col = row['usecol']

            if time < 54: # only first 54 time steps have truth data
                pst.observation_data.at[row.name, 'obsval'] = ts_heads.loc[time, col]
                pst.observation_data.at[row.name, 'standard_deviation'] = ts_heads.loc[time, 'std']
                pst.observation_data.at[row.name, 'weight'] = ts_heads.loc[time, 'weight']
        
        elif oname == 'arr-awq':
            kper = int(row['kper'])
            kstp = int(row['kstp'])
            i = int(row['i'])
            j = int(row['j'])
            
            if kper == 1:
                pst.observation_data.at[row.name, 'obsval'] = AWq.loc[(kper, kstp, i, j), 'AWq']
                pst.observation_data.at[row.name, 'standard_deviation'] = AWq.loc[(kper, kstp, i, j), 'std']
                pst.observation_data.at[row.name, 'weight'] = AWq.loc[(kper, kstp, i, j), 'weight']
        
        elif oname == 'arr-spq':
            kper = int(row['kper'])
            kstp = int(row['kstp'])
            i = int(row['i'])
            j = int(row['j'])
            
            if kper == 1:
                pst.observation_data.at[row.name, 'obsval'] = SPq.loc[(kper, kstp, i, j), 'SPq']
                pst.observation_data.at[row.name, 'standard_deviation'] = SPq.loc[(kper, kstp, i, j), 'std']
                pst.observation_data.at[row.name, 'weight'] = SPq.loc[(kper, kstp, i, j), 'weight']
        
        elif oname == 'arr-h':
            i = int(row['i'])
            j = int(row['j'])
            kper = int(row['kper'])
            kstp = int(row['kstp'])
            idom = ib[0, i-1, j-1]

            if (kper < 3) and (kstp in np.arange(1,53, 10)):
                # if idom == 1 & (i % 10 == 0) & (j % 10 == 0):
                #     pst.observation_data.at[row.name, 'obsval'] = heads[i, j] - (HEAD_offset*2) # 1m below ground level
                #     pst.observation_data.at[row.name, 'standard_deviation'] = heads_std[i, j]*20
                #     pst.observation_data.at[row.name, 'weight'] = 1/(heads_std[i, j]*20)
                # else:
                pst.observation_data.at[row.name, 'obsval'] = heads[i, j]
                pst.observation_data.at[row.name, 'standard_deviation'] = heads_std[i, j]
                pst.observation_data.at[row.name, 'weight'] = heads_weight[i, j]
                pst.observation_data.at[row.name, 'obgnme'] = 'less_' + row['obgnme']
            
                
            # else: # elsewhere in the model domain (less penalty on heads here; more uncertain about heads here)
            #     if (kper < 3) & (i % 10) & (j % 10):
            #         pst.observation_data.at[row.name, 'obsval'] = heads[i, j]
            #         pst.observation_data.at[row.name, 'standard_deviation'] = heads_std[i, j]
            #         pst.observation_data.at[row.name, 'weight'] = heads_weight[i, j]
            #         # pst.observation_data.at[row.name, 'obgnme'] = 'less_' + row['obgnme']
            
        elif oname == 'arr-confq':
            kper = int(row['kper'])
            kstp = int(row['kstp'])
            i = int(row['i'])
            j = int(row['j'])
            
            pst.observation_data.at[row.name, 'obsval'] = CONFq.loc[(kper, kstp, i, j), 'CONFq']
            pst.observation_data.at[row.name, 'standard_deviation'] = CONFq.loc[(kper, kstp, i, j), 'std']
            pst.observation_data.at[row.name, 'weight'] = CONFq.loc[(kper, kstp, i, j), 'weight']
            pst.observation_data.at[row.name, 'obgnme'] = 'greater_' + row['obgnme']
        
        elif oname == 'budget':
            kper = int(row['kper'])
            col = row['usecol']

            if (col in ['confined']): # needs to be greater than recharge because negative
                pst.observation_data.at[row.name, 'obsval'] = budget.loc[kper, col]
                pst.observation_data.at[row.name, 'obgnme'] = 'greater_' + row['obgnme']
                pst.observation_data.at[row.name, 'standard_deviation'] = budget.loc[kper, 'conf-std']
                pst.observation_data.at[row.name, 'weight'] = budget.loc[kper, 'conf-weight']
                
            elif (col in ['inflow']) & (kper < 3): # recharge plus confined less than awanui max flows
                pst.observation_data.at[row.name, 'obsval'] = budget.loc[kper, col]
                pst.observation_data.at[row.name, 'obgnme'] = 'less_' + row['obgnme']
                pst.observation_data.at[row.name, 'standard_deviation'] = budget.loc[kper, 'inflow-std']
                pst.observation_data.at[row.name, 'weight'] = budget.loc[kper, 'inflow-weight']

        elif oname == 'recession':
            kper = int(row['kper'])
            col = row['usecol']
            
            if (kper == 2) & (col in ['pk4', 'pk4-spr-diff', 'pk4-conf-diff']):
                pst.observation_data.at[row.name, 'obsval'] = recession_truth.loc[kper, col]
                pst.observation_data.at[row.name, 'standard_deviation'] = recession_truth.loc[kper, 'std']
                pst.observation_data.at[row.name, 'weight'] = recession_truth.loc[kper, 'weight']
    
    print('*' * 40)
    print(f'Number of observations: {pst.nobs}')
    print(f'Number of non-zero observations by group:')
    print(pst.observation_data[pst.observation_data['weight'] > 0].groupby('obgnme').size())
    print('Total number of non-zero observations:', pst.observation_data[pst.observation_data['weight'] > 0].shape[0])
    print('*' * 40)

    # RENAME OBSERVATION GROUPS TO AVOID SUBSTRING CONFLICTS --------------
    # PEST++ does pattern matching on phi factors, so obg1 matches obg1, obg10, obg11, etc.
    # Add unique numeric prefix to each group (e.g., 1obg0, 2obg1, 3obg2, 4obg10, etc.)
    # For less_/greater_ groups, map them to the same base group
    print('\nRenaming observation groups to avoid PEST++ pattern matching conflicts...')

    # First, identify base groups (without less_/greater_ prefix)
    base_groups = set()
    for obgnme in pst.observation_data['obgnme'].unique():
        if obgnme.startswith('less_') or obgnme.startswith('greater_'):
            # Extract base group name
            base = obgnme.split('_', 1)[1]
            base_groups.add(base)
        else:
            base_groups.add(obgnme)

    # Create mapping for base groups with unique prefixes
    base_to_new = {}
    for idx, base in enumerate(sorted(base_groups), start=1):
        base_to_new[base] = f"{idx}{base}"

    # Create mapping of old to new names
    old_to_new = {}
    for obgnme in pst.observation_data['obgnme'].unique():
        if obgnme.startswith('less_'):
            base = obgnme.split('_', 1)[1]
            old_to_new[obgnme] = f"less_{base_to_new[base]}"
        elif obgnme.startswith('greater_'):
            base = obgnme.split('_', 1)[1]
            old_to_new[obgnme] = f"greater_{base_to_new[base]}"
        else:
            old_to_new[obgnme] = base_to_new[obgnme]

    # Apply renaming
    pst.observation_data['obgnme'] = pst.observation_data['obgnme'].map(old_to_new)

    print('\nRenamed observation groups to avoid PEST++ pattern matching conflicts:')
    for old, new in old_to_new.items():
        print(f'  {old:20s} -> {new}')

    # PHI FACTORS ------------------------------------------------------
    # map dict for phi factors
    phi_obgnme = {}
    for obgnme in pst.observation_data['obgnme'].unique():
        oname = pst.observation_data[pst.observation_data['obgnme'] == obgnme]['oname'].values[0]
        weights = pst.observation_data[pst.observation_data['obgnme'] == obgnme]['weight'].values
        # only consider non-zero weights
        if np.sum(weights) == 0:
            continue

        phi_factor = PHI_OBS.get(oname, None)
        if phi_factor:
            # Keep the full obgnme including any less_/greater_ prefix
            phi_obgnme[obgnme] = phi_factor

    df = pd.DataFrame(list(phi_obgnme.items()), columns=['obgnme', 'phi_factor'])
    # Sort by obgnme for consistency
    df = df.sort_values('obgnme')
    df.to_csv(os.path.join(TEMP_DIR, 'phi_factors.csv'), index=False, header=False)
    pst.pestpp_options['ies_phi_factor_file'] = 'phi_factors.csv'


    ## ADD FORECASTS ------------------------------------------------------
    forecast_obgnme = ['ts-flux', 'budget']
    
    mask = pst.observation_data.loc[:,'oname'].isin(forecast_obgnme)
    forecasts = pst.observation_data[mask].index.tolist()
    
    pst.pestpp_options['forecasts'] = forecasts

    # PARAMETER COVARIANCE ---------------------------------------------
    cov_file = os.path.join(TEMP_DIR, 'prior_cov.jcb')
    cov = pf.build_prior(
        fmt='coo', 
        filename=cov_file)
    cov.to_binary(cov_file)
    print(f"Saved covariance matrix to {cov_file}")

    # WRITE PEST -------------------------------------------------------
    # chnage onames to integers
    # pst.observation_data.index = [int(n) for n in list(np.arange(0,len(pst.observation_data['obgnme'])))]

    print("Writing PEST template file...")
    pst.control_data.noptmax = 0 # just run parameters in file
    pst.observation_data = pst.observation_data.fillna(0)
    pst_file = f'{MODEL_NAME}.pst'
    final_pst = os.path.join(TEMP_DIR, pst_file)

    pst.write(final_pst, version=2)

    # OBSERVATION COVARIANCE ---------------------------------------------

    def update_time_covariance_direct(grp_data, v_time, var_dict, cov_matrix, name_to_idx):
        """Add temporal covariance directly to matrix"""
        onames = grp_data.obsnme.tolist()
        times = [int(time) for time in grp_data.time.values]

        if len(onames) < 2:
            # Single observation - just diagonal
            variance = var_dict[onames[0]]
            idx = name_to_idx[onames[0]]
            cov_matrix[idx, idx] = variance
            return

        # temporal coordinates
        coords = np.zeros((len(times), 3))
        coords[:, 0] = times  # Time goes in x-coordinate

        cov_array = v_time.covariance_matrix(
            x=coords[:, 0],
            y=coords[:, 1],  # All zeros
            names=onames
        )
        
        # Get variances for this group from the dictionary
        variances = np.array([var_dict[name] for name in onames])
        std_devs = np.sqrt(variances)

        # Scale correlation matrix by standard deviations
        cov_scaled = cov_array.x * np.outer(std_devs, std_devs)

        # Write directly to matrix (full symmetric matrix)
        for i, name_i in enumerate(onames):
            idx_i = name_to_idx[name_i]
            for j, name_j in enumerate(onames):
                idx_j = name_to_idx[name_j]
                cov_matrix[idx_i, idx_j] = cov_scaled[i, j]

    def update_spatial_covariance_direct(grp_data, v_space, var_dict, cov_matrix, name_to_idx):
        """Add spatial covariance directly to matrix"""
        onames = grp_data.obsnme.tolist()
        
        if len(onames) < 2:
            # Single observation - just diagonal
            variance = var_dict[onames[0]]
            idx = name_to_idx[onames[0]]
            cov_matrix[idx, idx] = variance
            return
        
        i_list = [int(i)-1 for i in grp_data['i'].values]
        j_list = [int(j)-1 for j in grp_data['j'].values]

        # Get spatial model coordinates from cell indices
        x_coords = sr.xcentergrid[i_list, j_list]
        y_coords = sr.ycentergrid[i_list, j_list]

        cov_array = v_space.covariance_matrix(
            x=x_coords,
            y=y_coords,
            names=onames
        )
        
        # Get variances for this group from the dictionary
        variances = np.array([var_dict[name] for name in onames])
        std_devs = np.sqrt(variances)

        # Scale correlation matrix by standard deviations
        cov_scaled = cov_array.x * np.outer(std_devs, std_devs)

        # Write directly to matrix (full symmetric matrix)
        for i, name_i in enumerate(onames):
            idx_i = name_to_idx[name_i]
            for j, name_j in enumerate(onames):
                idx_j = name_to_idx[name_j]
                cov_matrix[idx_i, idx_j] = cov_scaled[i, j]

    # ===================================================================
    # BUILD OBSERVATION COVARIANCE MATRIX
    # ===================================================================

    # Filter to non-zero weighted observations only
    nz_obs = pst.observation_data[pst.observation_data.weight > 0].copy()
    print(f"\n{'='*60}")
    print(f"BUILDING OBSERVATION COVARIANCE MATRIX")
    print(f"{'='*60}")
    print(f"Total observations: {len(pst.observation_data)}")
    print(f"Non-zero weighted observations: {len(nz_obs)}")
    print(f"Reduction: {100 * (1 - len(nz_obs)/len(pst.observation_data)):.1f}%")

    # Get variance lookup for non-zero weighted obs only
    v = pyemu.Cov.from_obsweights(final_pst)
    all_variances = v.x[:, 0]
    all_names = v.row_names
    var_dict = dict(zip(all_names, all_variances))

    # Pre-allocate the full covariance matrix
    nz_obs_names = nz_obs.obsnme.tolist()
    n = len(nz_obs_names)
    cov_matrix = np.zeros((n, n))
    name_to_idx = {name: i for i, name in enumerate(nz_obs_names)}

    print(f"Pre-allocated {n}x{n} covariance matrix")

    # Process each observation group
    obs_groups = nz_obs.obgnme.unique()

    for obgnme in obs_groups:
        group_mask = nz_obs.obgnme == obgnme
        grp_data = nz_obs[group_mask].copy()
        oname = grp_data['oname'].unique()[0]
        
        if len(grp_data) == 0:
            continue

        print(f"\nProcessing group: {oname} ({obgnme}) ({len(grp_data)} obs)")

        # Array observations - temporal and spatial correlation with subsampling
        if oname in ['arr-awq', 'arr-spq']:
            # create (i, j) subsampling lists
            i_j_list = list(zip(grp_data['i'], grp_data['j']))
            # i_list = [int(i) for i in grp_data['i'].unique()]
            # j_list = [int(j) for j in grp_data['j'].unique()]
            times_unique = [int(t) for t in grp_data['time'].unique()]
            
            print(f"  Original grid: {len(i_j_list)} locations")
            # print(f"  Original grid: {len(i_list)} x {len(j_list)} = {len(i_list)*len(j_list)} cells")
            
            # Temporal correlation within each subsampled (i,j) location
            n_temporal = 0
            for i, j in i_j_list:
                mask_ij = (grp_data['i'] == i) & (grp_data['j'] == j)
                grp_data_ij = grp_data[mask_ij]
                
                if len(grp_data_ij) > 0:
                    update_time_covariance_direct(
                        grp_data_ij, v_time, var_dict, cov_matrix, name_to_idx
                    )
                    n_temporal += 1

            
            print(f"  Processed temporal correlation for {n_temporal} locations")
            
            # Spatial correlation at each time step (using subsampled cells)
            n_spatial = 0
            for time in times_unique:
                mask_time = grp_data['time'] == str(time)
                grp_data_time = grp_data[mask_time]
                
                if len(grp_data_time) > 0:
                    update_spatial_covariance_direct(
                        grp_data_time, v_coarse, var_dict, cov_matrix, name_to_idx
                    )
                    n_spatial += 1
            
            print(f"  Processed spatial correlation for {n_spatial} time steps")
        
        # Time series observations - temporal correlation only
        elif oname in ['ts-heads', 'ts-sprflux']:
            usecols = grp_data['usecol'].unique()
            n_cols = 0
            for col in usecols:
                mask_col = grp_data['usecol'] == col
                grp_data_col = grp_data[mask_col]
                
                if len(grp_data_col) > 0:
                    update_time_covariance_direct(
                        grp_data_col, v_time, var_dict, cov_matrix, name_to_idx
                    )
                    n_cols += 1
            
            print(f"  Processed temporal correlation for {n_cols} locations")
        
        # Other observation groups - just diagonal
        else:
            n_diag = 0
            for oname in grp_data.obsnme:
                if oname in var_dict and oname in name_to_idx:
                    idx = name_to_idx[oname]
                    cov_matrix[idx, idx] = var_dict[oname]
                    n_diag += 1
            
            print(f"  Processed diagonal for {n_diag} observations")

    # ===================================================================
    # CREATE COV OBJECT AND SAVE
    # ===================================================================

    print(f"\nCreating pyemu.Cov object...")
    obscov = pyemu.Cov(x=cov_matrix, names=nz_obs_names, isdiagonal=False)

    print(f"Saving to binary format...")
    obscov.to_binary(os.path.join(TEMP_DIR, 'obscov.jcb'))
    pst.pestpp_options['observation_covariance'] = 'obscov.jcb'
    pst.write(final_pst, version=2)

    # Calculate matrix statistics
    n_nonzero = np.count_nonzero(cov_matrix)
    matrix_density = 100 * n_nonzero / (n**2)

    print(f"\n{'='*60}")
    print(f"OBSERVATION COVARIANCE SUMMARY")
    print(f"{'='*60}")
    print(f"Matrix size: {n} x {n}")
    print(f"Non-zero entries: {n_nonzero:,}")
    print(f"Matrix density: {matrix_density:.2f}%")
    print(f"Saved to: obscov.jcb")
    print(f"{'='*60}\n")
    
    # RUN PESTPP-IES --------------------------------------------------
    pyemu.os_utils.run("pestpp-ies {0}".format(pst_file), cwd=TEMP_DIR)
    
    # PRIOR PARAMETER COVARIANCE --------------------------------------------------
    print("Adding parameter covariance...")

    pe_full = pf.draw(num_reals=NREALS, use_specsim=True)
    pe_full.enforce() # enforces parameter bounds

    full_pe = os.path.join(TEMP_DIR, 'prior_pe.jcb')
    pe_full.to_binary(full_pe) #writes the parameter ensemble to binary file
    
    pst.write(final_pst, version=2)

    print('done')
    

if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations