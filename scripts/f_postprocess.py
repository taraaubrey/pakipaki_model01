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

run_name = 'local_run0'
o_d = f'models/{run_name}/{run_name}'
temp_d = f'models/{run_name}/{run_name}'
m_d_pr = f'models/{run_name}/pest/master_ies_prior'
m_d = f'models/{run_name}/pest/master_ies'
output = f'models/{run_name}/figures'

f_act = os.path.join(m_d, f"{run_name}.phi.actual.csv")
f_meas = os.path.join(m_d, f"{run_name}.phi.meas.csv")

pst = pyemu.Pst(os.path.join(m_d, f"{run_name}.pst"))

# observation ensemble
pr_oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(m_d, f"{run_name}.0.obs.csv"))
pt_oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(m_d, f"{run_name}.{pst.control_data.noptmax}.obs.csv"))
noise = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(m_d, f"{run_name}.obs+noise.csv"))


plot_phi_comparison(f_act,f_meas, output)

pst.plot('1to1')

plot_phi_hist(pr_oe, pt_oe)

plot_obs_hist(pr_oe, pt_oe)

print('done')


