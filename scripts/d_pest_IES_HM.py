import os
import pyemu
import pandas as pd

from setup import *

def main():
    
    # update pest file
    pst_name = f"{MODEL_NAME}.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))
    pst.pestpp_options.update(PEST_PP_OPTIONS)
    pst.control_data.noptmax = NOPTMAX
    pst.svd_data.maxsing = MAXSING
    pst.write(os.path.join(TEMP_DIR, pst_name))
    
    # start ies
    num_workers = os.cpu_count()
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_ies')
    pyemu.os_utils.start_workers(
        worker_dir=TEMP_DIR, # the folder which contains the "template" PEST dataset
        exe_rel_path=f'pestpp-ies', #the PEST software version we want to run
        pst_rel_path=pst_name, # the control file to use with PEST
        num_workers=num_workers, #how many agents to deploy
        worker_root= os.path.join(TEMP_DIR, '..'),
        master_dir=m_d, #the manager directory
        )
    
    pst.plot()


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations