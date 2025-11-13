import os
import pyemu

from setup import *

def main():

    pst_name = f"{MODEL_NAME}.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))
    pst.pestpp_options["tie_by_group"] = True
    new_pst_name = f"{MODEL_NAME}_sen.pst"
    pst.write(os.path.join(TEMP_DIR, new_pst_name))
    
    # run sen
    num_workers = os.cpu_count()
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_sen')
    pyemu.os_utils.start_workers(
        worker_dir=TEMP_DIR, # the folder which contains the "template" PEST dataset
        exe_rel_path=f'pestpp-sen', #the PEST software version we want to run
        pst_rel_path=new_pst_name, # the control file to use with PEST
        num_workers=num_workers, #how many agents to deploy
        worker_root= os.path.join(TEMP_DIR, '..'),
        master_dir=m_d, #the manager directory
        )


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations