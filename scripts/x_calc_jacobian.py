import os
import pandas as pd
import numpy as np
import pyemu

Y_sim_fn = r"C:\Users\tfo46\e_Python\e_projects\pakipaki_01\models\local_run24\pest\master_ies\local_run24.3.obs.csv"
B_fn = r"C:\Users\tfo46\e_Python\e_projects\pakipaki_01\models\local_run24\pest\master_ies\local_run24.3.par.csv"
pst_fn = r"C:\Users\tfo46\e_Python\e_projects\pakipaki_01\models\local_run24\pest\master_ies\local_run24.pst"

Y_df = pd.read_csv(Y_sim_fn, index_col=0)
B_df = pd.read_csv(B_fn, index_col=0)
Ne = B.shape[0]

Y_sim = Y_df.values # obs ensemble
B = B_df.values # parameter ensemble
Ne = B.shape[0] # ensemble size

# means
Y_mean = Y_sim.mean(axis=0)  # mean of each observation ensemble
B_mean = B.mean(axis=0)  # mean of each parameter ensemble

# deviations matrix
Y_dev = np.subtract(Y_sim,Y_mean) / np.sqrt(Ne - 1)  # deviations from mean, scaled
B_dev = np.subtract(B,B_mean) / np.sqrt(Ne - 1)  # deviations from mean, scaled

# calculate Jacobian
# obs ( Nobs x Ne )  , param ( Ne x Npar )  => Jacobian ( Nobs x Npar )
B_inv = np.linalg.pinv(B_dev)  # pseudo-inverse of B_dev (Moore-Penrose)
J = Y_dev.T @ B_inv.T  # Jacobian matrix
# J_dir = r'C:\Users\tfo46\e_Python\e_projects\pakipaki_01\models\local_run24\figures'
# np.savetxt(os.path.join(J_dir, 'Jacobian.txt'), J, delimiter=",")

J_df = pd.DataFrame(J, columns=B_df.columns, index=Y_df.columns)


# sort the matrix 
J_sort = np.sort(J, axis=0)
J_sort = np.sort(J_sort, axis=1)

# SVD
"""
very small singular values indicate directions in parameter space that have little effect on obs space. Non-identifiable or weakly identifiable parameters.
"""
U, s, Vt = np.linalg.svd(J)


# ------------------------------------------------------------
# FOSM : Schur complement
# ------------------------------------------------------------
pst = pyemu.Pst(pst_fn)
# using pyemu to create Jco object
jco =pyemu.Jco(J, row_names=Y_df.columns, col_names=B_df.columns)

sc = pyemu.Schur(jco=jco, pst=pst)


# ------------------------------------------------------------
"""
POSTERIOR PARAMETER UNCERTAINTY
Shows the percentage decrease in parameter uncertainty expected through calibration.
"""
# ------------------------------------------------------------

par_sum = sc.get_parameter_summary().sort_values("percent_reduction",ascending=False)
par_sum.iloc[0:50,:]['percent_reduction'].plot(kind='bar')
plt.title('Percent Reduction')

par_sum.iloc[0:50,:][['prior_var','post_var']].plot(kind='bar');

# parameters not informed by calibration
uninformed_pars = par_sum[par_sum['percent_reduction'] < 1]

# which parameter groups are not being informed?
uninf_par_names = uninformed_pars.index.values
pst.parameter_data.loc[uninf_par_names, 'pargp'].unique()


# ------------------------------------------------------------
# forecast uncertainty
# ------------------------------------------------------------
forecast_df = sc.get_forecast_summary()

# get the forecast summary then make a bar chart of the percent_reduction column
fig = plt.figure()
ax = plt.subplot(111)
ax = forecast_df.percent_reduction.plot(kind='bar',ax=ax,grid=True)
ax.set_ylabel("percent uncertainy\nreduction from calibration")
ax.set_xlabel("forecast")
plt.show()


# ------------------------------------------------------------
# parameter contributions to forecast uncertainty
# ------------------------------------------------------------

par_contrib = sc.get_par_group_contribution()

base = par_contrib.loc["base",:]
par_contrib = 100.0 * (base - par_contrib) / base
par_contrib.sort_index().head()



# ------------------------------------------------------------
# ------------------------------------------------------------
# DATA WORTH
# ------------------------------------------------------------
# ------------------------------------------------------------
obs_nz = sc.pst.observation_data.loc[pst.nnz_obs_names]
nn_obs_dict={}
for obsgp in pst.nnz_obs_groups:
    values = obs_nz.loc[obs_nz.obgnme==obsgp, 'obsnme'].tolist()
    nn_obs_dict[obsgp]=values

df_worth = sc.get_removed_obs_importance(nn_obs_dict)