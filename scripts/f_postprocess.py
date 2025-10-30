import os
import pandas as pd
from scipy import stats
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec

from f_helpers import *
import pyemu
import flopy

run_name = 'local_run3'
o_d = f'models/{run_name}/{run_name}'
temp_d = f'models/{run_name}/{run_name}'
m_d_pr = f'models/{run_name}/pest/master_ies_prior'
m_d = f'models/{run_name}/pest/master_ies'
output = f'models/{run_name}/figures'

f_act = os.path.join(m_d, f"{run_name}.phi.actual.csv")
f_meas = os.path.join(m_d, f"{run_name}.phi.meas.csv")
# observations
pr_oe_fn = os.path.join(m_d, f"{run_name}.0.obs.csv")
noise_fn = os.path.join(m_d, f"{run_name}.obs+noise.csv")
obs_data_fn = os.path.join(m_d, f"{run_name}.obs_data.csv")
# parameters
par_data_fn = os.path.join(m_d, f"{run_name}.par_data.csv")

# open pst file
pst = pyemu.Pst(os.path.join(m_d, f"{run_name}.pst"))

# input observation data
obs_data = pd.read_csv(obs_data_fn)
par_data = pd.read_csv(par_data_fn)
headers = obs_data['obsnme'].to_list()

# summarize parameters
parameter_summary = {}
parameter_summary['count'] = par_data.groupby(['pargp'])['parubnd'].describe()['count']
parameter_summary['parlbnd'] = par_data.groupby(['pargp'])['parlbnd'].describe()['mean']
parameter_summary['parubnd'] = par_data.groupby(['pargp'])['parubnd'].describe()['mean']
parameter_summary['transform'] = par_data.groupby(['pargp'])['partrans'].describe()['top']
parameter_summary['style'] = par_data.groupby(['pargp'])['parchglim'].describe()['top']
pd.DataFrame(parameter_summary).to_csv(os.path.join(output, 'output.parameter_summary.csv'))

obs_summary = obs_data[['obgnme', 'weight', 'standard_deviation']].groupby('obgnme').describe()
obs_summary.to_csv(os.path.join(output, 'output.allobs_summary.csv'))

nnzobs_summary = obs_data[obs_data['weight'] > 0][['obgnme', 'weight', 'standard_deviation']].groupby('obgnme').describe()
nnzobs_summary.to_csv(os.path.join(output, 'output.nnzobs_summary.csv'))


# # observation ensemble
# pr_oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(m_d, f"{run_name}.0.obs.csv"))
# pt_oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(m_d, f"{run_name}.{pst.control_data.noptmax}.obs.csv"))
# noise = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(m_d, f"{run_name}.obs+noise.csv"))


# plot_phi_comparison(f_act,f_meas, output)

# pst.plot('1to1')
# pst.plot('obs_v_sim')
# pst.plot('obs')

# # plot_phi_hist(pr_oe, pt_oe)
# plot_obs_budget_hist(pr_oe, pt_oe, headers, obs_data, output)

# # plot_phi_hist(pr_oe, pt_oe)
# plot_obs_tsheads_linesplot(pr_oe, pt_oe, headers, obs_data, output)


print('done')