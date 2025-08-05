import os
import pyemu

from setup import *

def main():

    pst_name = f"{MODEL_NAME}.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))
    pst.control_data.noptmax = -1
    pst.write(os.path.join(TEMP_DIR, pst_name))

    prior_cov = pyemu.Cov.from_parameter_data(pst)
    parensemble = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=prior_cov, num_reals=250)
    # ensure that the samples respect parameter bounds in the pst control file
    parensemble.enforce()
    parensemble.to_csv(os.path.join(TEMP_DIR,"sweep_in.csv"))
    
    num_workers = os.cpu_count()

    # the master directory
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_mc')

    pyemu.os_utils.start_workers(
        worker_dir=TEMP_DIR, # the folder which contains the "template" PEST dataset
        exe_rel_path=f'{exe_path}/pestpp-swp.exe', #the PEST software version we want to run
        pst_rel_path=pst_name, # the control file to use with PEST
        num_workers=num_workers, #how many agents to deploy
        worker_root= os.path.join(TEMP_DIR, '..'),
        master_dir=m_d, #the manager directory
        )


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations