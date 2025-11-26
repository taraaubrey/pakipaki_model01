import os
import pyemu
import pandas as pd
import numpy as np

from setup import *

def main():
    
    # update pest file
    pst_name = f"{MODEL_NAME}.pst"
    pst_path = os.path.join(TEMP_DIR, pst_name)
    pst = pyemu.Pst(pst_path)
    pst.pestpp_options.update(PEST_PP_OPTIONS)
    pst.control_data.noptmax = NOPTMAX

    if MAXSING is not None:
        pst.svd_data.maxsing = MAXSING
    pst.write(pst_path, version=2)

    print(f"Saving to binary format...")
    
    obscov_path = os.path.join(TEMP_DIR, 'obscov.jcb')
    v = pyemu.Cov.from_binary(obscov_path)
    ## update logic here
    all_variances = v.x[:, 0]
    all_names = v.names
    # get index where non-zero obs are
    pred_obs = pst.observation_data[pst.observation_data].coadd py()
    nz_variances = [all_variances[all_names.index(name)] for name in pred_obs.obsnme.tolist()]
    for name in pred_obs.obsnme.tolist():
        idx = all_names.index(name)
        v.x[idx,0] = 600  # set new variance value

    obscov = pyemu.Cov(x=v.x, names=all_names, isdiagonal=True)

    obscov.to_binary(obscov_path)
    pst.pestpp_options['observation_covariance'] = 'obscov.jcb'
    pst.write(os.path.join(TEMP_DIR, pst_name))

    # update the parameter file and obs ensembles as needed
    
    # start ies
    num_workers = os.cpu_count()
    m_d=os.path.join(os.path.join(TEMP_DIR, '..'), 'master_0')
    pyemu.os_utils.start_workers(
        worker_dir=TEMP_DIR, # the folder which contains the "template" PEST dataset
        exe_rel_path=f'pestpp-ies', #the PEST software version we want to run
        pst_rel_path=pst_name, # the control file to use with PEST
        num_workers=num_workers, #how many agents to deploy
        worker_root= os.path.join(TEMP_DIR, '..'),
        master_dir=m_d, #the manager directory
        )
    
    # pst.plot()


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations