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
                            longnames=True, # set False if using PEST/PEST_HP
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
    
    vh_space = pyemu.geostats.GauVario(
        name='head_gs', contribution=1.0, a=1000, anisotropy=1.0, bearing=0.0)
    vh_fine = pyemu.geostats.GauVario(
        name='head_gs', contribution=1.0, a=100, anisotropy=1.0, bearing=0.0)

    v_time = pyemu.geostats.ExpVario(
        name='time_gs', contribution=1.0, a=10, anisotropy=1.0, bearing=0.0)

    # geostatistical structure for spatially varying parameters
    fine_gs = pyemu.geostats.GeoStruct(variograms=v_fine, transform='log')
    grid_gs = pyemu.geostats.GeoStruct(variograms=v_space, transform='log')
    
    h_fine_gs = pyemu.geostats.GeoStruct(variograms=vh_fine, transform='none')
    h_grid_gs = pyemu.geostats.GeoStruct(variograms=vh_space, transform='none')
    
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
    add_ts_parameters(
        pf, TEMP_DIR,
        datetimes=dts,
        name='ghb-ts-aw',
        f=f'{MODEL_NAME}.ghb_aw_heads.csv',
        bounds=[-1e-1, 1e-1],
        parnames=list(aw_df[0]),
        use_cols=np.arange(1, aw_df.shape[0]+1).tolist(),
        time_gs={'tgs': time_gs},
        )
    add_ts_parameters(
        pf, TEMP_DIR,
        datetimes=dts,
        name='ghb-ts-pw',
        f=f'{MODEL_NAME}.ghb_pw_heads.csv',
        bounds=[-1e-1, 1e-1],
        parnames=list(pw_df[0]),
        use_cols=np.arange(1, pw_df.shape[0]+1).tolist(),
        time_gs={'tgs': time_gs},
        )
    add_ts_parameters(
        pf, TEMP_DIR,
        datetimes=dts,
        name='ghb-ts-spring',
        f=f'{MODEL_NAME}.ghb_spring_heads.csv',
        bounds=[-1e-1, 1e-1],
        parnames=list(sp_df[0]),
        use_cols=np.arange(1, sp_df.shape[0]+1).tolist(),
        time_gs={'tgs': time_gs},
        )
    add_ts_parameters(
        pf, TEMP_DIR,
        datetimes=dts,
        name='ghb-ts-conf',
        f=f'{MODEL_NAME}.ghb_conf_heads.csv',
        bounds=[-1e-1, 1e-1],
        parnames=['head_ts'],
        use_cols=[1],
        time_gs={'tgs': time_gs},
        )
    # # set up pst file
    # K ----------------------------------------------
    # shallow aquifer
    define_mult_array(
        pf, TEMP_DIR,
        tag=f'{MODEL_NAME}.npf_k_layer',
        ib=ib,
        grid_gs=fine_gs,
        grid_suffix='-fngr',
        lb=KH_PRIOR['lb'], ub=KH_PRIOR['ub'], 
        pp_space=5,
        pp_gs = grid_gs,
        lays=np.arange(NLAY).tolist()
        )
    # S ----------------------------------------------
    define_mult_array(
        pf, TEMP_DIR,
        tag=f'{MODEL_NAME}.sto_ss_layer',
        ib=ib,
        grid_gs=grid_gs,
        grid_suffix='-fngr',
        lb=SS_PRIOR['lb'], ub=SS_PRIOR['ub'], 
        pp_space=5,
        pp_gs = grid_gs,
        lays=np.arange(NLAY).tolist()
        )

    # RECHARGE ------------------------------------------------------
  
    rcha(
        pf, TEMP_DIR, ib,
        tag=f'{MODEL_NAME}.rcha_recharge',
        grid_gs=fine_gs,
        grid_suffix='-fngr',
        rch_bounds=[RCH['lb'], RCH['ub']],
        pp_gs=grid_gs,
        pp_space=5,
        )
    
    # GHBS ------------------------------------------------------
    h_ghb_gs = {
        'fngr': h_fine_gs,
        # 'crgr': h_grid_gs,
    }
    
    ghb(pf, TEMP_DIR,
        name='ghb-aw',
        tag=f'{MODEL_NAME}.ghbaw_stress_period_data',
        grid_gs={'fngr': fine_gs},
        head_gs=h_ghb_gs,
        cond_bounds=[GHB_SW['cond_lb'], GHB_SW['cond_ub']],
        head_bounds=[GHB_SW['head_lb'], GHB_SW['head_ub']],
        )
    
    ghb(pf, TEMP_DIR,
        name='ghb-pw',
        tag=f'{MODEL_NAME}.ghb_pw_stress_period_data',
        grid_gs={'fngr': fine_gs},
        head_gs=h_ghb_gs,
        cond_bounds=[GHB_SW['cond_lb'], GHB_SW['cond_ub']],
        head_bounds=[GHB_SW['head_lb'], GHB_SW['head_ub']],
        )
    
    ghb(pf, TEMP_DIR,
        name='ghb-spring',
        tag=f'{MODEL_NAME}.ghbspr_stress_period_data',
        grid_gs={'fngr': fine_gs},
        head_gs=h_ghb_gs,
        cond_bounds=[GHB_SW['cond_lb'], GHB_SW['cond_ub']],
        head_bounds=[GHB_SW['head_lb'], GHB_SW['head_ub']],
        )
    
    conf_ghb(pf, TEMP_DIR,
        name='ghb-conf',
        tag=f'{MODEL_NAME}.ghb_conf_stress_period_data',
        cond_gs=grid_gs,
        head_gs=h_grid_gs,
        cond_bounds=[GHB_CONF['cond_lb'], GHB_CONF['cond_ub']],
        head_bounds=[GHB_CONF['head_lb'], GHB_CONF['head_ub']],
        )

    # OBSERVATIONS ------------------------------------------------------
    print('SETTING OBSERVATIONS...')

    # add stress period head/fluxes observations
    index_cols = ['time']
    obs_use_cols = ['pk4', 'pk4-spr-diff', 'pk4-conf-diff']
    pf.add_observations(
        'output.sample_heads.csv', # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=obs_use_cols, # skip the index column
        prefix='ts-heads',
        obsgp='ts-heads',
    )

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

    # add recession
    recession = pd.read_csv(os.path.join(TEMP_DIR, "output.sample_recession_rates.csv"))
    index_cols = ['idx']
    use_cols=list(recession.columns[1:])
    pf.add_observations(
        "output.sample_recession_rates.csv", # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols, # skip the index column
        prefix='recession',
        obsgp='recession',
    )

    # add budget
    budget = pd.read_csv(os.path.join(TEMP_DIR, "output.budget.csv"))
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
    pattern = os.path.join(TEMP_DIR, "output.heads_lyr*.dat")
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

    ###############################################################################
    ts_heads = pd.read_csv(os.path.join(TRUTH_DIR, 'output.sample_heads.truth.csv'), index_col=0)
    AWq = pd.read_csv(os.path.join(TRUTH_DIR, "output.GHB_AW_fluxes.truth.csv"), index_col=[0, 1, 2, 3, 4])
    recession_truth = pd.read_csv(os.path.join(TRUTH_DIR, 'output.sample_recession_rates.truth.csv'), index_col=0)
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

    budget['awanui'] = 0 
    budget['confined'] = 0
    budget['poukawa'] = 0
    budget['std'] = 3
    budget['weight'] = 1/budget['std']

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
            
            if kper == 1:
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
        
        elif obgnme == 'budget':
            kper = int(row['kper'])
            col = row['usecol']
            
            if col in ['awanui', 'confined', 'poukawa']:
                pst.observation_data.at[row.name, 'obsval'] = budget.loc[kper-1, col]
                pst.observation_data.at[row.name, 'standard_deviation'] = budget.loc[kper-1, 'std']
                pst.observation_data.at[row.name, 'weight'] = budget.loc[kper-1, 'weight']

                if col in ['awanui','poukawa']:
                    pst.observation_data.at[row.name, 'obgnme'] = 'less_' + row['obgnme']

                elif col in ['confined']:
                    pst.observation_data.at[row.name, 'obgnme'] = 'greater_' + row['obgnme']
        elif obgnme == 'recession':
            idx = int(row['idx'])
            col = row['usecol']
            
            try:
                pst.observation_data.at[row.name, 'obsval'] = recession_truth.loc[idx, col]
                pst.observation_data.at[row.name, 'standard_deviation'] = recession_truth.loc[idx, 'std']
                pst.observation_data.at[row.name, 'weight'] = recession_truth.loc[idx, 'weight']

            except:
                continue
    

    # create phi factor file
    phi_dict = {
        'head-arr': 100,
        'ts-heads': 200,
        'recession': 50,
        'AWq': 1,
        'budget': 10,
    }
    # to dataframe with no column names or index
    df = pd.DataFrame(list(phi_dict.items()))
    df.to_csv(os.path.join(TEMP_DIR, 'phi_factors.csv'), index=False, header=False)

    ## ADD FORECASTS ------------------------------------------------------
    forecast_obgnme = ['ts-flux', 'budget']
    
    mask = pst.observation_data.loc[:,'obgnme'].isin(forecast_obgnme)
    forecasts = pst.observation_data[mask].index.tolist()
    
    pst.pestpp_options['forecasts'] = forecasts
    pst.pestpp_options['ies_phi_factor_file'] = 'phi_factors.csv'

    # PARAMETER COVARIANCE ---------------------------------------------
    cov = pf.build_prior(fmt='coo', filename=os.path.join(TEMP_DIR, "prior_cov.jcb"))

    # cov = pyemu.Cov.from_parameter_data(pst)

    # Add cross-covariance between ghb-conf and ghb-aw HEAD parameters
    # (conductance parameters remain uncorrelated)
    print("Adding cross-covariance between ghb-conf and ghb-aw head parameters...")

    # Get all parameter names
    par_names = pst.parameter_data.index.tolist()

    # Identify ghb-conf head parameters (exclude conductance)
    conf_head_pars = [p for p in par_names if 'ghb-conf' in p and 'head' in p]
    conf_head_pars += [p for p in par_names if 'ghb-ts-conf' in p]

    # Identify ghb-aw head parameters (exclude conductance)
    aw_head_pars = [p for p in par_names if 'ghb-aw' in p and 'head' in p]
    aw_head_pars += [p for p in par_names if 'ghb-ts-aw' in p]

    print(f"Found {len(conf_head_pars)} conf head parameters and {len(aw_head_pars)} aw head parameters")

    # Set correlation coefficient
    correlation = 0.6

    cov_copy = cov.x.copy()
    # Get parameter standard deviations from the covariance matrix diagonal
    par_std = {}
    for i, par_name in enumerate(cov.row_names):
        par_std[par_name] = np.sqrt(cov_copy[i, i])

    # Create index mapping for covariance matrix
    name_to_idx = {name: i for i, name in enumerate(cov.row_names)}

    # Helper function to extract indices from parameter name
    def extract_indices(par_name):
        """Extract idx0:idx1:idx2 from parameter name, or return None if not present"""
        # Look for pattern like idx0:idx1:idx2 in the parameter name
        match = re.search(r'idx(\d+):idx(\d+):idx(\d+)', par_name)
        if match:
            return (match.group(1), match.group(2), match.group(3))
        return None

    # Create dictionaries mapping indices to parameter names for faster lookup
    conf_by_indices = {}
    for par in conf_head_pars:
        indices = extract_indices(par)
        if indices and par in name_to_idx:
            conf_by_indices[indices] = par

    aw_by_indices = {}
    for par in aw_head_pars:
        indices = extract_indices(par)
        if indices and par in name_to_idx:
            aw_by_indices[indices] = par

    # Set cross-covariance only between parameters with matching indices
    n_cross_cov = 0
    for indices, conf_par in conf_by_indices.items():
        # Check if there's a matching aw parameter with same indices
        if indices in aw_by_indices:
            aw_par = aw_by_indices[indices]

            if conf_par not in par_std or aw_par not in par_std:
                continue

            idx_conf = name_to_idx[conf_par]
            idx_aw = name_to_idx[aw_par]

            # Calculate cross-covariance: cov[i,j] = correlation * std[i] * std[j]
            cross_cov_value = correlation * par_std[conf_par] * par_std[aw_par]

            # Set symmetric entries in covariance matrix
            cov.x[idx_conf, idx_aw] = cross_cov_value
            cov.x[idx_aw, idx_conf] = cross_cov_value
            n_cross_cov += 1

    print(f"Set {n_cross_cov} cross-covariance terms with correlation = {correlation}")
    print(f"(Only correlated parameters with matching spatial indices)")

    # Save covariance matrix to file
    cov_file = os.path.join(TEMP_DIR, 'prior_cov.jcb')
    cov.to_binary(cov_file)
    print(f"Saved covariance matrix to {cov_file}")


    # WRITE PEST -------------------------------------------------------
    
    print("Writing PEST template file...")
    pst.control_data.noptmax = 0 # just run parameters in file
    pst.observation_data = pst.observation_data.fillna(0)
    pst_file = f'{MODEL_NAME}.pst'
    final_pst = os.path.join(TEMP_DIR, pst_file)

    pst.write(final_pst, version=2)

    # RUN PESTPP-IES --------------------------------------------------
    pyemu.os_utils.run("pestpp-ies {0}".format(pst_file), cwd=TEMP_DIR)
    
    # PRIOR PARAMETER COVARIANCE --------------------------------------------------
    print("Adding parameter covariance...")
    
    pe = pf.draw(num_reals=NREALS_PRIOR, use_specsim=True)
    pe.enforce() # enforces parameter bounds
    
    pe_f = os.path.join(TEMP_DIR, 'prior_pe.jcb')
    pe.to_binary(pe_f) #writes the parameter ensemble to binary file
    
    pst.write(final_pst, version=2)

    print('done')
    

if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations