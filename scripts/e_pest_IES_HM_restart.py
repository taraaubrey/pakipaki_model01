import os
import pyemu
import pandas as pd

from setup import *

def main():

    TEMP_DIR = r"C:\Users\tfo46\e_Python\e_projects\pakipaki_01\models\local_run30\pest\local_run30_restart_temp"

    pst_name = f"local_run30.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))

    # pest options
    pst.pestpp_options['ies_parameter_ensemble'] = 'local_run30.19.par.csv'
    pst.pestpp_options['ies_observation_ensemble'] = 'local_run30.obs+noise.csv'
    pst.pestpp_options["ies_restart_observation_ensemble"] = 'local_run30.19.obs.csv'
    pst.pestpp_options['overdue_giveup_fac'] = 10
    pst.pestpp_options['overdue_giveup_minutes'] = 15
    # pst.pestpp_options["ies_n_iter_reinflate"] = 5, 7
    # pst.pestpp_options["ies_reinflate_factor"] = 1.1

    pst.control_data.noptmax = 3

    # write a new pst file
    new_pst = f"local_run30_restart.pst"
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