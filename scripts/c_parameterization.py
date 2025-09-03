import os
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
        pf, TEMP_DIR,
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
        pf, TEMP_DIR,
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
    # define_mult_array(pf, TEMP_DIR,
    #         tag=f'{MODEL_NAME}.rcha_recharge',
    #         sr=sr,
    #         ib=ib,
    #         grid_gs=grid_gs,
    #         lb=0.01, ub=1, ulb=0, uub=1e-1,
    #         add_coarse=True)

    wel(pf, TEMP_DIR,
        name='mbr',
        tag=f'{MODEL_NAME}.wel_mbr_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        q_bounds=[1e-2, 1e2],
        # q_ultbounds=[0.01, 10]
        )

    wel(pf, TEMP_DIR,
        name='influx',
        tag=f'{MODEL_NAME}.wel_influx_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        q_bounds=[0.1, 10],
        # q_ultbounds=[0.01, 10]
        )

    wel(pf, TEMP_DIR,
        name='outflux',
        tag=f'{MODEL_NAME}.wel_outflux_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        q_bounds=[0.1, 10],
        # q_ultbounds=[0.1, 90]
        )

    ghb(pf, TEMP_DIR,
        name='ghb_aw',
        tag=f'{MODEL_NAME}.ghbaw_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        cond_bounds=[1e-4, 1e4],
        # cond_ultbounds=[1e-6, 1e6],
        head_bounds=[-2, 2],
        # head_ultbounds=[5, 20]
        )
    
    ghb(pf, TEMP_DIR,
        name='ghb_pw',
        tag=f'{MODEL_NAME}.ghb_pw_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        cond_bounds=[1e-4, 1e4],
        # cond_ultbounds=[1e-6, 1e6],
        head_bounds=[-2, 2],
        # head_ultbounds=[5, 20]
        )
    
    ghb(pf, TEMP_DIR,
        name='ghb_spring',
        tag=f'{MODEL_NAME}.ghbspr_stress_period_data',
        fine_gs=fine_gs,
        grid_gs=grid_gs,
        constant_gs=None,
        cond_bounds=[1e-4, 1e4],
        # cond_ultbounds=[1e-6, 1e6],
        head_bounds=[-2, 2],
        # head_ultbounds=[5, 15]
        )
    
    ghb(pf, TEMP_DIR,
        name='ghb_conf',
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
    # add_ts_parameters(pf, 'spring', TEMP_DIR, f'{MODEL_NAME}.spring_stage.csv', -1, 1)
    # add_ts_parameters(pf, 'pw', TEMP_DIR, f'{MODEL_NAME}.ghb_pw_heads.csv', -1, 1)
    # add_ts_parameters(pf, 'conf', TEMP_DIR, f'{MODEL_NAME}.ghb_conf_heads.csv', -1, 1)
    # add_ts_parameters(pf, 'aw', TEMP_DIR, f'{MODEL_NAME}.awanui_stage.csv', -1, 1)
    
    
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
    obs_results = pd.read_csv(os.path.join(TRUTH_DIR, 'obs_results_truth.csv'))
    std_obs_results  = pd.read_csv(os.path.join(TRUTH_DIR, f'obs_results_std.csv'))
    index_cols = ['kper', 'kstp']
    use_cols = list(obs_results.columns.values)[5:]
    pf.add_observations(
        'obs_results.csv', # this is being read from the template file not the above truth file
        index_cols=index_cols,
        use_cols=use_cols, # skip the index column
        prefix='ts',
        obsgp='ts',
    )

    # FORWARD RUN SCRIPT --------------------------------------------------
    # pst = pf.build_pst()
    pf.mod_sys_cmds.append("mf6") #do this only once
    # pf.mod_sys_cmds.append(f"mp7 {MODEL_NAME}.mpsim") #do this only once
    # pst = pf.build_pst()

    sample_rel = os.path.relpath(SAMPLES, TEMP_DIR)
    # post-processing to get observations
    pf.add_py_function(
        f"{SCRIPTS_DIR}/helpers.py",
        f"extract_heads_and_budget(model_name='{MODEL_NAME}')", is_pre_cmd=False)
    pf.add_py_function(
        f"{SCRIPTS_DIR}/helpers.py",
        f"extract_spring_obs(gwf=None, model_name='{MODEL_NAME}', samples_path=r'{sample_rel}', conf_h_path='{MODEL_NAME}.ghb_conf_heads.csv', aw_h_path='{MODEL_NAME}.awanui_stage.csv', pw_h_path='{MODEL_NAME}.ghb_pw_heads.csv', spring_h_path='{MODEL_NAME}.spring_stage.csv', pw_q_path='{MODEL_NAME}.poukawa_flow_m3d.csv')", is_pre_cmd=False)

    pst = pf.build_pst()
    # pst_file = f'{MODEL_NAME}.pst'
    # pst.write(os.path.join(TEMP_DIR, pst_file),version=2)

    ###############################################################################

    # build weights from standard deviations
    init_values = obs_results.set_index(['kper', 'kstp']).to_dict(orient='index')
    std_values = std_obs_results.set_index(['kper', 'kstp']).to_dict(orient='index')
    for col in use_cols:
        for kper_kstp in list(std_values.keys()):
            kper = int(kper_kstp[0])
            kstp = int(kper_kstp[1])
            if col in std_values[kper_kstp].keys():
                initial_val = float(init_values[kper_kstp][col])
                std = float(std_values[kper_kstp][col])
                if std > 0:
                    weight = 1.0 / std
                    pst.observation_data.loc[
                        f'oname:ts_otype:lst_usecol:{col}_kper:{kper}_kstp:{kstp}','obsval'] = initial_val
                    pst.observation_data.loc[
                        f'oname:ts_otype:lst_usecol:{col}_kper:{kper}_kstp:{kstp}','standard_deviation'] = std
                    pst.observation_data.loc[
                        f'oname:ts_otype:lst_usecol:{col}_kper:{kper}_kstp:{kstp}','weight'] = weight
                else:
                    pst.observation_data.loc[
                        f'oname:ts_otype:lst_usecol:{col}_kper:{kper}_kstp:{kstp}','obsval'] = initial_val
                    pst.observation_data.loc[
                        f'oname:ts_otype:lst_usecol:{col}_kper:{kper}_kstp:{kstp}','standard_deviation'] = std
                    pst.observation_data.loc[
                        f'oname:ts_otype:lst_usecol:{col}_kper:{kper}_kstp:{kstp}','weight'] = 0
            else:
                pst.observation_data.loc[
                        f'oname:ts_otype:lst_usecol:{col}_kper:{kper}_kstp:{kstp}','obsval'] = initial_val
                pst.observation_data.loc[
                        f'oname:ts_otype:lst_usecol:{col}_kper:{kper}_kstp:{kstp}','standard_deviation'] = std
                pst.observation_data.loc[
                    f'oname:ts_otype:lst_usecol:{col}_kper:{kper}_kstp:{kstp}','weight'] = 0

    # adjust the obs which are supposed to be greater than constraints
    val = pst.observation_data.loc[
                    f'oname:ts_otype:lst_usecol:awswgwdiff_kper:1_kstp:1','obgnme']
    pst.observation_data.loc[
                    f'oname:ts_otype:lst_usecol:awswgwdiff_kper:1_kstp:1','obgnme'] = 'greater_' + val
    
    val = pst.observation_data.loc[
                    f'oname:ts_otype:lst_usecol:pwswgwdiff_kper:1_kstp:1','obgnme']
    pst.observation_data.loc[
                    f'oname:ts_otype:lst_usecol:pwswgwdiff_kper:1_kstp:1','obgnme'] = 'greater_' + val
    
    ## ADD FORECASTS ------------------------------------------------------
    forecasts = [
        'oname:ts_otype:lst_usecol:pk4head_kper:3_kstp:1',
        'oname:ts_otype:lst_usecol:pk4head_kper:4_kstp:1',
        'oname:ts_otype:lst_usecol:pk4head_kper:4_kstp:2',
        'oname:ts_otype:lst_usecol:pk4head_kper:4_kstp:3',
        'oname:ts_otype:lst_usecol:pk4head_kper:4_kstp:4',
        'oname:ts_otype:lst_usecol:pk4head_kper:4_kstp:5',
        'oname:ts_otype:lst_usecol:pk4springdiff_kper:3_kstp:1',
        'oname:ts_otype:lst_usecol:pk4springdiff_kper:4_kstp:1',
        'oname:ts_otype:lst_usecol:pk4springdiff_kper:4_kstp:2',
        'oname:ts_otype:lst_usecol:pk4springdiff_kper:4_kstp:3',
        'oname:ts_otype:lst_usecol:pk4springdiff_kper:4_kstp:4',
        'oname:ts_otype:lst_usecol:pk4springdiff_kper:4_kstp:5',
        'oname:ts_otype:lst_usecol:springq_kper:1_kstp:1',
        'oname:ts_otype:lst_usecol:springq_kper:2_kstp:1',
        'oname:ts_otype:lst_usecol:springq_kper:2_kstp:2',
        'oname:ts_otype:lst_usecol:springq_kper:2_kstp:3',
        'oname:ts_otype:lst_usecol:springq_kper:2_kstp:4',
        'oname:ts_otype:lst_usecol:springq_kper:2_kstp:5',
        'oname:ts_otype:lst_usecol:springq_kper:3_kstp:1',
        'oname:ts_otype:lst_usecol:springq_kper:4_kstp:1',
        'oname:ts_otype:lst_usecol:springq_kper:4_kstp:2',
        'oname:ts_otype:lst_usecol:springq_kper:4_kstp:3',
        'oname:ts_otype:lst_usecol:springq_kper:4_kstp:4',
        'oname:ts_otype:lst_usecol:springq_kper:4_kstp:5',
    ]
    pst.pestpp_options['forecasts'] = forecasts


    # WRITE PEST -------------------------------------------------------
    print("Writing PEST template file...")
    pst.control_data.noptmax = 0 # just run parameter values in the file
    pst_file = f'{MODEL_NAME}.pst'
    pst.write(os.path.join(TEMP_DIR, pst_file), version=2)

    # RUN PESTPP-IES --------------------------------------------------
    pyemu.os_utils.run("pestpp-ies {0}".format(pst_file), cwd=TEMP_DIR)
    
    # PRIOR PARAMETER COVARIANCE --------------------------------------------------
    print("Adding parameter covariance...")
    pe_f = os.path.join(TEMP_DIR, 'prior_pe.jcb')

    pe = pf.draw(num_reals=NREALS_PRIOR, use_specsim=True)
    pe.enforce() # enforces parameter bounds
    pe.to_binary(pe_f) #writes the parameter ensemble to binary file
    
    final_pst = os.path.join(TEMP_DIR, pst_file)
    pst.write(final_pst, version=2)
    
    
    
    # TEST RUN PESTPP-IES --------------------------------------------------
    # test_pst = os.path.join(TEMP_DIR,"test.pst")
    # # grab the first realization from the ensemble
    # pst.parameter_data.loc[:,"parval1"] = pe.loc[pe.index[20],pst.par_names].values
    # pst.control_data.noptmax = 0
    # pst.write(test_pst, version=2)
    
    # pyemu.os_utils.run("pestpp-glm test.pst", cwd=TEMP_DIR)
    
    # PLOT TEST REALIZATION --------------------------------------------------
    # df = pd.read_csv(os.path.join(TEMP_DIR,"mult2model_info.csv"))
    # kh1_df = df.loc[df.model_file.str.contains("npf_k_layer1"),:]

    # org_arr = np.loadtxt(os.path.join(TEMP_DIR,kh1_df.org_file.iloc[0]))
    # inp_arr = np.loadtxt(os.path.join(TEMP_DIR,kh1_df.model_file.iloc[1]))
    # mlt_arrs = [np.loadtxt(os.path.join(TEMP_DIR,afile)) for afile in kh1_df.mlt_file]
    
    # arrs = [org_arr]
    # arrs.extend(mlt_arrs)
    # arrs.append(inp_arr)
    
    # names = ["org"]
    # names.extend([mf.split('.')[0].split('_')[-1] for mf in kh1_df.mlt_file])
    # names.append("MF6 input")
    
    # fig,axes = plt.subplots(1,kh1_df.shape[0]+2,figsize=(5*kh1_df.shape[0]+2,5))
    # for i,ax in enumerate(axes.flatten()):
    #     arr = np.log10(arrs[i])
    #     arr[ib[0]==0] = np.nan
    #     cb = ax.imshow(arr)
    #     plt.colorbar(cb,ax=ax)
    #     ax.set_title(names[i],loc="left")
    # plt.tight_layout()

    # RUN TEST IES ------------------------------------------------------
    # print("Running PESTPP-IES test...")
    # pst = pyemu.Pst(final_pst)
    # pst_ies = os.path.join(TEMP_DIR,"test_ies.pst")
    # pst.pestpp_options["ies_num_reals"] = 30  # starting with a real small ensemble!
    # pst.pestpp_options['ies_parameter_ensemble'] = 'prior_pe.jcb'
    # pst.pestpp_options["ies_bad_phi_sigma"] = 2.0 #middle ground value
    # pst.control_data.noptmax = -2
    
    # pst.write(pst_ies, version=2)
    # pyemu.os_utils.run(f"pestpp-ies test_ies.pst", cwd=TEMP_DIR)

    # # load and check phi
    # pst = pyemu.Pst(pst_ies)
    # pst.phi

    # EVAL --------------------------------------------------------------------------

    # pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_file))
    # pst.phi
    # pst.phi_components

    # nnz_phi_components = {k:pst.phi_components[k] for k in pst.nnz_obs_groups} 
    # nnz_phi_components

    # phicomp = pd.Series(nnz_phi_components)
    # plt.pie(phicomp, labels=phicomp.index.values);

    # print('Target phi:',pst.nnz_obs)
    # print('Current phi:', pst.phi)

    # figs = pst.plot(kind="1to1");
    # pst.res.loc[pst.nnz_obs_names,:]
    # plt.show()

if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations