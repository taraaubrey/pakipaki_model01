import os
import pyemu

from setup import *


def main():

    pst_name = f"{MODEL_NAME}.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))
    pst.control_data.noptmax = -2
    # rewrite the contorl file!
    pst.write(os.path.join(TEMP_DIR, pst_name))

    num_workers = os.cpu_count()

    # the master directory
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_local')

    pyemu.os_utils.start_workers(TEMP_DIR, # the folder which contains the "template" PEST dataset
                                'pestpp-glm', #the PEST software version we want to run
                                pst_name, # the control file to use with PEST
                                num_workers=16, #how many agents to deploy
                                worker_root= os.path.join(TEMP_DIR, '..'),
                                master_dir=m_d, #the manager directory
                                )


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations