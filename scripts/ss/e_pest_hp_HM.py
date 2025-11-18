import os
import pyemu

from setup import *

def main():

    pst_name = f"{MODEL_NAME}.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))

    # pest options
    pst.pestpp_options["uptestmin"] = 20
    pst.control_data.noptmax = 3

    # write a new pst file
    new_pst_name = f"{MODEL_NAME}_hp.pst"
    pst.write(os.path.join(TEMP_DIR, new_pst_name))
    
    num_workers = os.cpu_count()

    # the master directory
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_hp')

    pyemu.os_utils.start_workers(
        worker_dir=TEMP_DIR, # the folder which contains the "template" PEST dataset
        exe_rel_path=f'pest_hp', #the PEST software version we want to run
        pst_rel_path=new_pst_name, # the control file to use with PEST
        num_workers=num_workers, #how many agents to deploy
        worker_root= os.path.join(TEMP_DIR, '..'),
        master_dir=m_d, #the manager directory
        )


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations