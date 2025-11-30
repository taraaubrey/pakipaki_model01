import os
import pyemu
import numpy as np
import shutil
import argparse
import sys

from setup import *

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Perform Pareto hypothesis testing with PEST++ IES.'
    )

    parser.add_argument('--post-iter', type=int, default=3,
                       help='Posterior iteration to use (default: 3)')
    parser.add_argument('--m_d', type=str, default='master_ies',
                       help='name of pest master directory (default: master_ies)')
    return parser.parse_args()



def main(posterior_iter, og_md):
    
    # update global run settings
    og_obscov = os.path.join(og_md, f'{MODEL_NAME}.{posterior_iter}.shrunk_res.cov')
    og_pst = f"{MODEL_NAME}.pst"
    pst = pyemu.Pst(os.path.join(og_md, og_pst))
    pst.pestpp_options['observation_covariance'] = og_obscov

    # update ensemble (same as restarting)
    for filename, argname in zip(
        [
            f'{MODEL_NAME}.{posterior_iter}.par.csv',
            f'{MODEL_NAME}.{posterior_iter}.obs.csv',
            f'{MODEL_NAME}.obs+noise.csv',
        ],
        [
            "ies_parameter_ensemble",
            "ies_restart_observation_ensemble",
            "ies_observation_ensemble"
        ]):
        
        renamed_filename = f"ht_"+filename
        # copy the original restart file from the prior master dir to the renamed filename in the template dir
        shutil.copy2(
            os.path.join(og_md, filename),
            os.path.join(TEMP_DIR, renamed_filename))
        #modify/set the pestpp option
        pst.pestpp_options[argname] = renamed_filename

    pst.control_data.noptmax = 1

    cycles = np.arange(0,CYCLES,1).tolist()
    pred_variances = np.linspace(MIN, MAX, CYCLES)

    for i, cycle in enumerate(cycles):
        RUN_NAME = f"{MODEL_NAME}_c{cycle}"
        pst_name = f"{RUN_NAME}.pst"
        pst_path = os.path.join(TEMP_DIR, pst_name)
        m_d=os.path.join(os.path.join(TEMP_DIR, '..'), f'master_c{cycle}')

        if cycle > 0:
            pst = pyemu.Pst(os.path.join(TEMP_DIR, og_pst))
            new_obscov_name = f'obscov_c{cycle}.jcb'
            
            v = pyemu.Cov.from_binary(os.path.join(TEMP_DIR, og_obscov))
            all_names = v.names
            all_variances = v.x[:, 0]

            # get index where non-zero obs are
            idx_names = nz_obs.index.tolist() + pst.forecast_names
            nz_variances = [all_variances[all_names.index(name)] for name in idx_names]
            var_x = np.array(nz_variances).reshape(-1, 1)

            for name in pst.forecast_names:
                idx = idx_names.index(name)
                var_x[idx] = pred_variances[i]**2
            
            obscov = pyemu.Cov(x=var_x, names=idx_names, isdiagonal=True)
            obscov.to_binary(os.path.join(TEMP_DIR, new_obscov_name))
            pst.pestpp_options['observation_covariance'] = new_obscov_name
        
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


if __name__ == "__main__":
    args = parse_args()

    TEMP_DIR = r'C:\Users\tfo46\e_Python\e_projects\pakipaki_01\models\63168_run43\run43\pest\master_ies'

    posterior_iter = args.post_iter
    md_name = args.m_d
    og_md = os.path.join(os.path.join(TEMP_DIR, '..'), md_name)

    main(posterior_iter, og_md)  # run the main function to build the model and extract observations