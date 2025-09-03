import os
import pyemu

from setup import *

def main():

    pst_name = f"{MODEL_NAME}_pp.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))
    pst.control_data.noptmax = -1

    # write a new pst file
    pst.write(os.path.join(TEMP_DIR, pst_name))
    
    num_workers = os.cpu_count()

    # the master directory
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_glm_1')

    pyemu.os_utils.start_workers(
        worker_dir=TEMP_DIR, # the folder which contains the "template" PEST dataset
        exe_rel_path=f'pestpp-glm', #the PEST software version we want to run
        pst_rel_path=pst_name, # the control file to use with PEST
        num_workers=num_workers, #how many agents to deploy
        worker_root= os.path.join(TEMP_DIR, '..'),
        master_dir=m_d, #the manager directory
        )


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations