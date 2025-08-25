import os
import pyemu

from setup import *

def main():

    pst_name = f"{MODEL_NAME}.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))
    pst.control_data.noptmax = -1

    # gte rid of high dimensional parameterization
    par = pst.parameter_data
    gr_pars = par.loc[par.pargp.apply(lambda x: "gr" in x and "sfr" not in x),"parnme"]
    par.loc[gr_pars,"partrans"] = "fixed"


    # write a new pst file
    pst.write(os.path.join(TEMP_DIR, pst_name))
    

    num_workers = os.cpu_count()

    # the master directory
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_mc')

    pyemu.os_utils.start_workers(
        worker_dir=TEMP_DIR, # the folder which contains the "template" PEST dataset
        exe_rel_path=f'pestpp-swp', #the PEST software version we want to run
        pst_rel_path=pst_name, # the control file to use with PEST
        num_workers=num_workers, #how many agents to deploy
        worker_root= os.path.join(TEMP_DIR, '..'),
        master_dir=m_d, #the manager directory
        )


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations