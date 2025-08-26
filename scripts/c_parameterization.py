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
    pyemu.os_utils.run("mf6", cwd=PEST_DIR)
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
    pst = define_mult_array(
        pf, TEMP_DIR,
        tag=f'{MODEL_NAME}.npf_k_layer',
        sr=sr,
        ib=ib,
        grid_gs=grid_gs,
        fine_gs=fine_gs,
        lb=0.001, ub=100, ulb=1e-4, uub=1e3,
        add_coarse=True,
        pp_space=10,
        lays=np.arange(NLAY-2).tolist()
        )
    # confining layer
    pst = define_mult_array(
        pf, TEMP_DIR,
        tag=f'{MODEL_NAME}.npf_k_layer',
        sr=sr,
        ib=ib,
        grid_gs=grid_gs,
        fine_gs=fine_gs,
        lb=0.01, ub=100, ulb=1e-6, uub=1e2,
        add_coarse=True,
        pp_space=10,
        lays=[NLAY-1]
        )
    #bottom layer
    pst = define_mult_array(
        pf, TEMP_DIR,
        tag=f'{MODEL_NAME}.npf_k_layer',
        sr=sr,
        ib=ib,
        grid_gs=grid_gs,
        fine_gs=fine_gs,
        lb=0.01, ub=100, ulb=1e2, uub=1e5,
        add_coarse=True,
        lays=[NLAY])

    # RECHARGE ------------------------------------------------------
    define_mult_array(pf, TEMP_DIR,
            tag=f'{MODEL_NAME}.rcha_recharge',
            sr=sr,
            ib=ib,
            grid_gs=grid_gs,
            lb=0.1, ub=10, ulb=0, uub=1e-1,
            add_coarse=True)

    wel(pf, TEMP_DIR,
        name='mbr',
        tag=f'{MODEL_NAME}.wel_mbr_stress_period_data',
        grid_gs=grid_gs,
        q_bounds=[0.01, 100],
        q_ultbounds=[0.01, 10])

    wel(pf, TEMP_DIR,
        name='influx',
        tag=f'{MODEL_NAME}.wel_influx_stress_period_data.txt',
        grid_gs=grid_gs,
        q_bounds=[0.01, 100],
        q_ultbounds=[0.01, 10])

    wel(pf, TEMP_DIR,
        name='outflux',
        tag=f'{MODEL_NAME}.wel_outflux_stress_period_data.txt',
        grid_gs=grid_gs,
        q_bounds=[0.01, 100],
        q_ultbounds=[0.01, 10])

    drn(pf, TEMP_DIR,
        name='drn_riv',
        tag=f'{MODEL_NAME}.drn_riv_stress_period_data.txt',
        grid_gs=grid_gs,
        cond_bounds=[0.01, 100],
        cond_ultbounds=[0.01, 1000],
        head_bounds=[-2, 2],
        head_ultbounds=[5, 15])

    chd(pf, TEMP_DIR,
        name='chd_pw',
        tag=f'{MODEL_NAME}.chd_pw_stress_period_data.txt',
        grid_gs=None,
        head_bounds=[-2, 2],
        head_ultbounds=[5, 15])

    chd(pf, TEMP_DIR,
        name='chd_conf',
        tag=f'{MODEL_NAME}.chd_conf_stress_period_data.txt',
        grid_gs=None,
        head_bounds=[-2, 2],
        head_ultbounds=[11, 15])

    # check
    # [f for f in os.listdir(TEMP_DIR) if f.endswith(".tpl")]
    pst = pf.build_pst()
    
    # OBSERVATIONS ------------------------------------------------------
    for f in ["cum.csv"]:
        df = pd.read_csv(os.path.join(TEMP_DIR, f))
        pf.add_observations(
            f,
            index_cols=list(df.columns.values)[0],
            use_cols=list(df.columns.values)[1:],
            prefix='flux',
            obsgp='flux')


    spring_f = os.path.join(TEMP_DIR, 'obs_results.csv')
    spring_obs = pd.read_csv(spring_f)
    pf.add_observations(
        'obs_results.csv',
        index_cols=[list(spring_obs.columns.values)[0]],
        use_cols=list(spring_obs.columns.values)[1:], # skip the index column
        prefix='heads',
        obsgp='heads',
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
        f"extract_spring_obs(gwf=None, model_name='{MODEL_NAME}', samples_path=r'{sample_rel}')", is_pre_cmd=False)

    pst = pf.build_pst()
    # pst_file = f'{MODEL_NAME}.pst'
    # pst.write(os.path.join(TEMP_DIR, pst_file),version=2)


    # WRITE PEST -------------------------------------------------------
    print("Writing PEST template file...")

    pst_file = f'{MODEL_NAME}.pst'
    pst.write(os.path.join(TEMP_DIR, pst_file), version=2)

    # RUN PESTPP-IES --------------------------------------------------
    print("Running PESTPP-IES 2nd time...")

    pyemu.os_utils.run("pestpp-ies {0}".format(pst_file), cwd=TEMP_DIR)

    # UPDATE OBSERVATIONS ------------------------------------------------------
    print("Updating observation weights...")
    pst.observation_data.loc[:, 'weight'] = 0

    tspringdf = pd.read_csv(os.path.join(TRUTH_DIR, 'obs_results.csv'), index_col=0)
    tcumdf = pd.read_csv(os.path.join(TRUTH_DIR, 'cum.csv'), index_col=0)

    for col in tspringdf.columns:
        if col in ['spring-head', 'pk2-head', 'pk4-head', 'spring-pk4-diff','spring-pk2-diff','pk4-pk2-diff']:
            new_col = col.replace('-', '_')
            pst.observation_data.loc[f'oname:heads_otype:lst_usecol:{new_col}_kper:0','obsval'] = tspringdf[col].iloc[0]
            pst.observation_data.loc[f'oname:heads_otype:lst_usecol:{new_col}_kper:0','weight'] = tspringdf[col].iloc[2]

    for col in tcumdf.columns:
        if col in ['drn']:
            new_col = col.replace('-', '_')  # replace spaces with underscores
            pst.observation_data.loc[f'oname:flux_otype:lst_usecol:{new_col}_totim:1','obsval'] = tcumdf[col].iloc[0]
            pst.observation_data.loc[f'oname:flux_otype:lst_usecol:{new_col}_totim:1','weight'] = tcumdf[col].iloc[2]

    
    ## ADD FORECASTS ------------------------------------------------------
    forecasts = [
        'oname:heads_otype:lst_usecol:spring_flux_kper:0',
        'oname:heads_otype:lst_usecol:spring_pk4_diff_kper:0',
    ]
    pst.pestpp_options['forecasts'] = forecasts
    
    
    # PRIOR PARAMETER COVARIANCE --------------------------------------------------
    print("Adding parameter covariance...")
    pe_f = os.path.join(TEMP_DIR, 'prior_pe.jcb')

    pe = pf.draw(num_reals=NREALS, use_specsim=True)
    pe.enforce() # enforces parameter bounds
    pe.to_binary(pe_f) #writes the parameter ensemble to binary file
    
    final_pst = os.path.join(TEMP_DIR, pst_file)
    pst.write(final_pst, version=2)
    
    # TEST RUN PESTPP-IES --------------------------------------------------
    test_pst = os.path.join(TEMP_DIR,"test.pst")
    # grab the first realization from the ensemble
    pst.parameter_data.loc[:,"parval1"] = pe.loc[pe.index[20],pst.par_names].values
    pst.control_data.noptmax = 0
    pst.write(test_pst, version=2)
    
    pyemu.os_utils.run("pestpp-glm test.pst", cwd=TEMP_DIR)
    
    # PLOT TEST REALIZATION --------------------------------------------------
    # df = pd.read_csv(os.path.join(TEMP_DIR,"mult2model_info.csv"))
    # kh1_df = df.loc[df.model_file.str.contains("npf_k_layer4"),:]

    # org_arr = np.loadtxt(os.path.join(TEMP_DIR,kh1_df.org_file.iloc[0]))
    # inp_arr = np.loadtxt(os.path.join(TEMP_DIR,kh1_df.model_file.iloc[0]))
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
    print("Running PESTPP-IES test...")
    pst = pyemu.Pst(final_pst)
    pst_ies = os.path.join(TEMP_DIR,"test_ies.pst")
    pst.pestpp_options["ies_num_reals"] = 30  # starting with a real small ensemble!
    pst.pestpp_options['ies_parameter_ensemble'] = 'prior_pe.jcb'
    pst.pestpp_options["ies_bad_phi_sigma"] = 2.0 #middle ground value
    pst.control_data.noptmax = -2
    
    pst.write(pst_ies, version=2)
    pyemu.os_utils.run(f"pestpp-ies test_ies.pst", cwd=TEMP_DIR)

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