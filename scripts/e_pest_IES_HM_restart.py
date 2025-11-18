import os
import pyemu
import pandas as pd

from setup import *

def main():

    # the master directory
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_ies')

    pst_name = f"local_run32.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))

    for filename, argname in zip([
        'local_run32.16.par.csv',
        'local_run32.16.obs.csv'
        'local_run32.obs+noise.csv',
        ],[
            "ies_parameter_ensemble",
            "ies_restart_observation_ensemble",
            "ies_observation_ensemble"]):

        renamed_filename = "restart_"+filename
        # copy the original restart file from the prior master dir to the renamed filename in the template dir
        shutil.copy2(
            os.path.join(m_d, filename),
            os.path.join(TEMP_DIR, renamed_filename))
        #modify/set the pestpp option
        pst.pestpp_options[argname] = renamed_filename

    # pest options
    pst.control_data.noptmax = 3

    # write a new pst file
    new_pst = f"{MODEL_NAME}_restart.pst"
    pst.write(os.path.join(TEMP_DIR, new_pst))
    
    num_workers = os.cpu_count()

    # the master directory
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_ies_restart')

    pyemu.os_utils.start_workers(
        worker_dir=TEMP_DIR, # the folder which contains the "template" PEST dataset
        exe_rel_path=f'pestpp-ies', #the PEST software version we want to run
        pst_rel_path=new_pst, # the control file to use with PEST
        num_workers=num_workers, #how many agents to deploy
        worker_root= os.path.join(TEMP_DIR, '..'),
        master_dir=m_d, #the manager directory
        )


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations