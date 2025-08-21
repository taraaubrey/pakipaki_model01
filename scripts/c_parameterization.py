import os
import stat
import shutil
import pyemu
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
    v_space = pyemu.geostats.ExpVario(contribution=1.0, #sill
                                        a=1000, # range of correlation; length units of the model. In our case 'meters'
                                        anisotropy=1.0, #name says it all
                                        bearing=0.0 #angle in degrees East of North corresponding to anisotropy ellipse
                                        )

    # geostatistical structure for spatially varying parameters
    grid_gs = pyemu.geostats.GeoStruct(variograms=v_space, transform='log')

    # plot the gs if you like:
    # _ = grid_gs.plot()

    # get the IDOMAIN array
    ib = gwf.dis.idomain.array

    # setup pilot points

    # # set up pst file
    # K ----------------------------------------------
    pst = define_mult_array(pf, TEMP_DIR,
            tag=f'{MODEL_NAME}.npf_k_layer',
            sr=sr,
            ib=ib,
            grid_gs=grid_gs,
            lb=0.01, ub=100, ulb=1e-4, uub=1e3,
            add_coarse=True,
            lays=[0, 1, 2, 3, 4, 5])
    pst = define_mult_array(pf, TEMP_DIR,
            tag=f'{MODEL_NAME}.npf_k_layer',
            sr=sr,
            ib=ib,
            grid_gs=grid_gs,
            lb=0.01, ub=100, ulb=1e-6, uub=1e2,
            add_coarse=True,
            lays=[6])
    pst = define_mult_array(pf, TEMP_DIR,
            tag=f'{MODEL_NAME}.npf_k_layer',
            sr=sr,
            ib=ib,
            grid_gs=grid_gs,
            lb=0.01, ub=100, ulb=1e2, uub=1e5,
            add_coarse=True,
            lays=[7])

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
    # [f for f in os.listdir(template_ws) if f.endswith(".tpl")]
    pst = pf.build_pst()
    
    # OBSERVATIONS ------------------------------------------------------
    for f in ["cum.csv"]:
        df = pd.read_csv(os.path.join(TEMP_DIR, f))
        pf.add_observations(
            f,
            index_cols=["totim"],
            # use_cols=list(df.columns.values)[1:],
            prefix=f.split('.')[0],
            obsgp=f.split(".")[0])


    spring_f = os.path.join(TEMP_DIR, 'obs_results.csv')
    spring_obs = pd.read_csv(spring_f)
    pf.add_observations(
        'obs_results.csv',
        # index_cols=[list(spring_obs.columns.values)[0]],
        use_cols=list(spring_obs.columns.values)[1:], # skip the index column
        prefix='springobs',
        obsgp='springobs',
    )

    # add the heads and budget observations
    # files = [f for f in os.listdir(TEMP_DIR) if f.startswith(f"{MODEL_NAME}_hdslay")]
    # for f in files:
    #     pf.add_observations(
    #         f,
    #         prefix=f.split(".")[0],
    #         obsgp=f.split(".")[0])

    # FORWARD RUN SCRIPT --------------------------------------------------
    # pst = pf.build_pst()
    pf.mod_sys_cmds.append("mf6") #do this only once
    # pf.mod_sys_cmds.append(f"mp7 {MODEL_NAME}.mpsim") #do this only once
    # pst = pf.build_pst()

    sample_rel = os.path.relpath(SAMPLES, TEMP_DIR)
    # post-processing to get observations
    pf.add_py_function(
        f"{SCRIPTS_DIR}helpers.py",
        f"extract_heads_and_budget(model_name='{MODEL_NAME}')", is_pre_cmd=False)
    pf.add_py_function(
        f"{SCRIPTS_DIR}helpers.py",
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

    pyemu.os_utils.run("pestpp-ies.exe {0}".format(pst_file), cwd=TEMP_DIR)

    # UPDATE OBSERVATIONS ------------------------------------------------------

    # obs = pst.observation_data
    pst.observation_data.loc[:, 'weight'] = 0

    tspringdf = pd.read_csv(os.path.join(TRUTH_DIR, 'obs_results.csv'))
    tcumdf = pd.read_csv(os.path.join(TRUTH_DIR, 'cum.csv'), index_col=0)

    for col in tspringdf.columns[2:]:
        pst.observation_data.loc[f'oname:springobs_otype:lst_usecol:{col}','obsval'] = tspringdf[col].iloc[0]
        pst.observation_data.loc[f'oname:springobs_otype:lst_usecol:{col}','weight'] = 1/(0.3*tspringdf[col].iloc[0])

    for col in tcumdf.columns:
        pst.observation_data.loc[f'oname:cum_otype:lst_usecol:{col}','obsval'] = tcumdf[col].iloc[0]
        pst.observation_data.loc[f'oname:cum_otype:lst_usecol:{col}','weight'] = tcumdf[col].iloc[2]

    print("Updating observation weights...")

    pst.write(os.path.join(TEMP_DIR, pst_file),version=2)

    # RUN PESTPP-IES --------------------------------------------------

    pyemu.os_utils.run("pestpp-ies.exe {0}".format(pst_file), cwd=TEMP_DIR)
    
    
    # EVAL --------------------------------------------------------------------------

    # pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_file))
    # pst.phi
    # pst.phi_components

    # nnz_phi_components = {k:pst.phi_components[k] for k in pst.nnz_obs_groups} # that's a dictionary comprehension there y'all
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