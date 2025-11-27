import os
import pyemu
import numpy as np
import shutil

from setup import *

def main():
    
    # update global run settings
    og_pst = f"{MODEL_NAME}.pst"
    og_obscov = f"obscov.jcb"

    pst = pyemu.Pst(os.path.join(TEMP_DIR, og_pst))
    pst.pestpp_options['observation_covariance'] = og_obscov
    pst.pestpp_options.update(PEST_PP_OPTIONS)
    pst.control_data.noptmax = 0
    if MAXSING is not None:
        pst.svd_data.maxsing = MAXSING

    
    cycles = np.arange(0,10,1)
    standard_deviation = [100, 80, 60, 40, 20, 10, 5, 2, 1, 0.5]

    for i, cycle in enumerate(cycles):
        RUN_NAME = f"{MODEL_NAME}_c{cycle}"
        pst_name = f"{RUN_NAME}.pst"
        pst_path = os.path.join(TEMP_DIR, pst_name)
        m_d=os.path.join(os.path.join(TEMP_DIR, '..'), f'master_c{cycle}')

        if cycle > 0:
            pst = pyemu.Pst(os.path.join(TEMP_DIR, last_pst_name))
            new_obscov_name = f'obscov_c{cycle}.jcb'
            
            v = pyemu.Cov.from_binary(os.path.join(TEMP_DIR, last_obscov_name))
            all_names = v.names
            all_variances = v.x[:, 0]

            # get index where non-zero obs are
            idx_names = nz_obs.index.tolist() + pst.forecast_names
            nz_variances = [all_variances[all_names.index(name)] for name in idx_names]
            var_x = np.array(nz_variances).reshape(-1, 1)

            for name in pst.forecast_names:
                idx = idx_names.index(name)
                var_x[idx] = standard_deviation[i]**2
            
            obscov = pyemu.Cov(x=var_x, names=idx_names, isdiagonal=True)
            obscov.to_binary(os.path.join(TEMP_DIR, new_obscov_name))
            pst.pestpp_options['observation_covariance'] = new_obscov_name

            # update ensemble (same as restarting)
            for filename, argname in zip(
                [
                    f'{last_pst_name}.{NOPTMAX}.par.csv',
                    f'{last_pst_name}.{NOPTMAX}.obs.csv'
                    f'{last_pst_name}.obs+noise.csv',
                ],
                [
                    "ies_parameter_ensemble",
                    "ies_restart_observation_ensemble",
                    "ies_observation_ensemble"
                ]):
                
                renamed_filename = f"c{cycle}_"+filename
                # copy the original restart file from the prior master dir to the renamed filename in the template dir
                shutil.copy2(
                    os.path.join(last_m_d, filename),
                    os.path.join(TEMP_DIR, renamed_filename))
                #modify/set the pestpp option
                pst.pestpp_options[argname] = renamed_filename
            
            num_workers = os.cpu_count()

        pst.write(pst_path, version=2)
        
        # start ies
        num_workers = os.cpu_count()
        pyemu.os_utils.start_workers(
            worker_dir=TEMP_DIR, # the folder which contains the "template" PEST dataset
            exe_rel_path=f'pestpp-ies', #the PEST software version we want to run
            pst_rel_path=pst_name, # the control file to use with PEST
            num_workers=num_workers, #how many agents to deploy
            worker_root= os.path.join(TEMP_DIR, '..'),
            master_dir=m_d, #the manager directory
            )
        
        # update last
        last_pst_name = pst_name
        last_obscov_name = pst.pestpp_options['observation_covariance']
        last_m_d = m_d


if __name__ == "__main__":
    main()  # run the main function to build the model and extract observations