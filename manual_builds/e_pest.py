import os
import pyemu
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import flopy

from a_setup import *


def plot_ensemble_arr(pe, tmp_d, numreals, modelname, par_dat):
    sim = flopy.mf6.MFSimulation.load(sim_ws=tmp_d, verbosity_level=0) #modflow.Modflow.load(fs.MODEL_NAM,model_ws=working_dir,load_only=[])
    gwf= sim.get_model()

    sr = pyemu.helpers.SpatialReference.from_namfile(
            os.path.join(tmp_d, f"{modelname}.nam"),
            delr=gwf.dis.delr.array, delc=gwf.dis.delc.array)

    pp_file=os.path.join(tmp_d,par_dat)
    df_pp = pyemu.pp_utils.pp_tpl_to_dataframe(os.path.join(tmp_d,f"{par_dat}.dat.tpl"))
    #same name order
    df_pp.sort_values(by='parnme', inplace=True)
    pe.sort_index(axis=1, inplace=True)


    fig = plt.figure(figsize=(12, 10))
    # generate random values
    for real in range(numreals):
        df_pp.loc[:,"parval1"] = pe.iloc[real,:].values
        # save a pilot points file

        # pyemu.pp_utils.write_pp_file(pp_file, df_pp)
        # interpolate the pilot point values to the grid
        ident_arr = pyemu.geostats.fac2real(f'{pp_file}.dat' , factors_file=pp_file+".fac",out_file=None)


        ax = fig.add_subplot(int(numreals/5)+1, 5, real+1, aspect='equal')
        mm = flopy.plot.PlotMapView(model=gwf, ax=ax, layer=0)
        ca = mm.plot_array(np.log10(ident_arr), masked_values=[1e30], vmin=0, vmax=0.1,)

        plt.scatter(df_pp.x, df_pp.y, marker='x', c='k', alpha=0.5)
        
        mm.plot_grid(alpha=0.5)
        mm.plot_inactive()
        ax.set_title(real+1)
        ax.set_yticks([])
        ax.set_xticks([])

    cb = plt.colorbar(ca, shrink=0.5)
    fig.tight_layout()
    return

pst_name = f"{MODEL_NAME}.pst"
pst = pyemu.Pst(os.path.join(TEMP_DIR, pst_name))

prior_cov = pyemu.Cov.from_parameter_data(pst)
x = prior_cov.as_2d.copy()
x[x==0] = np.nan
plt.imshow(x)
plt.colorbar();
parensemble = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst,cov=prior_cov, num_reals=250)
# ensure that the samples respect parameter bounds in the pst control file
parensemble.enforce()


par = pst.parameter_data
parnmes = par.loc[par.pargp=='npfklayer1pp'].parnme.values
pe_k = parensemble.loc[:,parnmes].copy()
TEMP_DIR = r"C:\Users\tfo46\e_Python\a_rbm\rapid-gwm-build\examples\manual_builds\models\local1\pest\local1_template"
# use the hbd convenienc function to plot several realisations
plot_ensemble_arr(pe_k, TEMP_DIR, 10, 'local1', 'npfklayer1pp_inst0pp')


v = pyemu.geostats.ExpVario(contribution=0.5, a=500, anisotropy=1.0, bearing=0.0)
gs = pyemu.utils.geostats.GeoStruct(variograms=[v])
pp_tpl = os.path.join(TEMP_DIR,"npfklayer1pp_inst0pp.dat.tpl")
cov = pyemu.helpers.geostatistical_prior_builder(pst=pst, struct_dict={gs:pp_tpl})
# display
plt.imshow(cov.to_pearson().x,interpolation="nearest")
plt.colorbar()
cov.to_dataframe().head()

parensemble = pyemu.ParameterEnsemble.from_gaussian_draw(pst=pst, cov=cov, num_reals=250,)
# ensure that the samples respect parameter bounds in the pst control file
parensemble.enforce()

pe_k = parensemble.loc[:,parnmes].copy()
# use the hbd convenienc function to plot several realisations
plot_ensemble_arr(pe_k, TEMP_DIR, 10, 'local1', 'npfklayer1pp_inst0pp')

print('here')