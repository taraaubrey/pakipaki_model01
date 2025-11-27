import os
import pyemu
import pandas as pd
import numpy as np

from setup import *

def main():
    
    # update pest file
    pst_name = f"{MODEL_NAME}.pst"
    pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))
    pst.pestpp_options.update(PEST_PP_OPTIONS)
    pst.control_data.noptmax = NOPTMAX

    if MAXSING is not None:
        pst.svd_data.maxsing = MAXSING
    pst.write(os.path.join(TEMP_DIR, pst_name))

    cycles = np.arange(0,10,1)
    standard_deviation = [100, 80, 60, 40, 20, 10, 5, 2, 1, 0.5]

    for i, cycle in enumerate(cycles):
        if cycle == 0:
            pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))
            pst.observation_data.loc[pst.forecast_names,'standard_deviation'] = standard_deviation[i]
            pst.observation_data.loc[pst.forecast_names,'weight'] = 1/pst.observation_data.loc[pst.forecast_names,'standard_deviation']

            # update obs cov
            obscov_path = os.path.join(TEMP_DIR, 'obscov.jcb')
            v = pyemu.Cov.from_binary(obscov_path)
            all_variances = v.x[:, 0]
            all_names = v.row_names

            # get index where non-zero obs are
            # get index where non-zero obs are
            pred_obs = pst.observation_data[pst.forecast_names]
            pred_variances = [all_variances[all_names.index(name)] for name in pred_obs.obsnme.tolist()]
            for name in pred_obs.obsnme.tolist():
                idx = all_names.index(name)
                v.x[idx,0] = standard_deviation[i]**2  # set new variance value


            idx_names = pst.observation_data.loc[pst.forecast_names]
            idx_variances = [all_variances[all_names.index(name)] for name in idx_names]
            obscov = pyemu.Cov(x=np.array(nz_variances).reshape(-1, 1), names=idx_names, isdiagonal=True)
            

    
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