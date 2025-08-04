import matplotlib.pyplot as plt
import os
import shutil
import gridit as gi
import numpy as np
import flopy as fp
import pandas as pd
from scipy.ndimage import uniform_filter

import helpers

from utils import *
from a_setup import *


# create directories if they do not exist
for d in [SPATIAL_DIR, FIG_DIR]:
    if not os.path.exists(d):
        os.makedirs(d)
# if dir exists delete it
if os.path.exists(MODEL_DIR):
    shutil.rmtree(MODEL_DIR)
os.makedirs(MODEL_DIR)  # create model directory

grid = gi.Grid.from_vector(DOMAIN, RES)

# save grid
grid_gpd = grid.cell_geodataframe()
# grid_gpd.to_file(os.path.join(SPATIAL_DIR, f'{MODEL_NAME}_grid.shp'), driver='ESRI Shapefile')

# top & bottom
top = grid.array_from_raster(TOP)
bottom = grid.array_from_raster(BOTTOM)  # get the bottom elevation

# open domain
arr = grid.array_from_vector(DOMAIN)
# arr = np.where(top.data > 15, 0, arr)

idomain = np.stack([arr] * NLAY, axis=0)
# confining layer active map
conf_area = grid.array_from_vector(CONF_AREA_ACTIVE)
idomain[-2] = np.where(conf_area.data == 0, 0, idomain[-2])  # set the confining layer active where the confining area is active
idomain[-1] = np.where(conf_area.data == 0, 0, idomain[-1])

#limeston_inactive
# limestone = grid.array_from_vector(LIMESTONE_INACTIVE)
# arr_limestone = np.where(limestone == 1, 0, arr)
# # add arr_limestone to idomain
# idomain = np.vstack([idomain, arr_limestone[np.newaxis, :, :], arr_limestone[np.newaxis, :, :]])

nrow = grid.shape[0]
ncol = grid.shape[1]
delr = np.ones(ncol) * RES
delc = np.ones(nrow) * RES


# dis
shallow_elevation = 0 # top of clay /confining layer
top = top.data  # use the top elevation for the first layer, only where idomain is active
bottom_shallow = np.where(bottom.data < shallow_elevation, shallow_elevation, bottom.data)  # set bottom elevation to top elevation where it is higher
thickness_shallow = top.data - bottom_shallow  # calculate the thickness of the layers
# thickness_deeper = np.ones_like(thickness_shallow) * 10  # set the thickness of the deeper layers

conf_elev = -10
grav_elev = -20
conf_bottom = np.where(bottom.data < conf_elev, conf_elev, bottom.data)  # set bottom elevation to conf_elev where it is higher
grav_bottom = np.where(bottom.data < grav_elev, grav_elev, bottom.data)  # set bottom elevation to conf_elev where it is higher

min_b = 1 # min thickness

thicknesses = []
botm = []
b0 = top.copy()
for i in range(NLAY):
    if i < 6:
        b = np.where(thickness_shallow/6 < min_b, min_b, thickness_shallow/6)
        ibotm = b0 - b  # calculate the bottom elevation for each layer
    elif i == 6:
        b = np.where(b0 - conf_bottom < min_b, min_b, b0 - conf_bottom)  # set the thickness of the confining layer
        ibotm = b0 - b  # calculate the bottom elevation for each layer
    else:   
        b = np.where(b0 - grav_bottom < min_b, min_b, b0 - grav_bottom)  # set the thickness of the confining layer
        ibotm = b0 - b  # calculate the bottom elevation for each layer
    thicknesses.append(b)
    botm.append(ibotm)
    b0 = ibotm.copy()  # update the bottom elevation for the next layer

botm = np.array(botm)
b_arr = np.array(thicknesses)
# plot_array_layers(b_arr, figsize=(15, 5), cmap='viridis', titles=None)

# update idomain to 0 where layer 7 and 8 have min thickness
idomain[-2] = np.where(b_arr[-2] <= min_b, 0, idomain[-2])  # layer 7
idomain[-1] = np.where(b_arr[-1] <= min_b, 0, idomain[-1])  # layer 8

# drains
drain_arr = grid.array_from_vector(DRAINS)
# grid_gpd['drn'] = drain_arr.data.flatten().tolist()
# grid_gpd.to_file(os.path.join(SPATIAL_DIR, f'{MODEL_NAME}_drn.shp'), driver='ESRI Shapefile')

spring_arr = grid.array_from_vector(SPRING)  # get the spring array
drain_arr = np.where(spring_arr.data == 1, 1, drain_arr.data)  # convert to binary array
drain_elev = grid.array_from_raster(TOP, resampling='min') * drain_arr * idomain[0]  # use the top elevation for drains, only where idomain is active
drain_input = extract_value_with_indices(drain_elev, layer=0, val_col='elev', mask_value=0)  # extract non-NaN values from the drain elevation array

# add top drains
drain_top = top * np.where(drain_elev > 0, 0, 1) * idomain[0]
drn_top_input = extract_value_with_indices(drain_top, layer=0, val_col='elev', mask_value=0)  # extract non-NaN values


# mbr
mbr_arr = grid.array_from_vector(MBR)
mbr_indices = []
for i in range(NLAY-2):  # remove mbr from bottom 2 layers
    in_idomarr, idom_i = get_interior_indices(idomain[i], layer=i)  # idomainget interior indices of the mbr area
    mbr_active = np.logical_and(mbr_arr.data, in_idomarr)  # mbr area that is active in the model domain
    mbr_indices.extend(get_indices(mbr_active, layer=i))
mbr_df = pd.DataFrame({'index': mbr_indices})

# chd - Poukawa boundary
top_min = grid.array_from_raster(TOP, resampling='min')
chd_pw_arr = grid.array_from_vector(POUKAWA_BOUNDARY)
chd_pw_active = np.logical_and(chd_pw_arr.data, in_idomarr)  # mbr area that is active in the model domain
chd_pw_active = chd_pw_active * top_min.data
chd_pw_indices = get_indices(chd_pw_active, layer=0, value=True)

chd_pw_df = pd.DataFrame({'index': [i[0] for i in chd_pw_indices]})


# influx boundaries
influx_arr = grid.array_from_vector(INFLUX_BOUNDARY)
influx_indices = []
for i in range(NLAY-2): # remove influx from bottom 2 layers
    in_idomarr, idom_i = get_interior_indices(idomain[i], layer=i)  # idomainget interior indices of the mbr area
    influx_active = np.logical_and(influx_arr.data, in_idomarr)  # mbr area that is active in the model domain
    influx_indices.extend(get_indices(influx_active, layer=i))
influx_df = pd.DataFrame({'index': influx_indices})

# outflux boundaries
outflux_arr = grid.array_from_vector(OUTFLUX_BOUNDARY)
outflux_indices = []
for i in range(NLAY-2):
    in_idomarr, idom_i = get_interior_indices(idomain[i], layer=i)  # idomainget interior indices of the mbr area
    outflux_active = np.logical_and(outflux_arr.data, in_idomarr)  # mbr area that is active in the model domain
    outflux_indices.extend(get_indices(outflux_active, layer=i))
outflux_df = pd.DataFrame({'index': outflux_indices})

# chd - confined boundary inflow
chd_conf_arr = idomain[-1]
chd_conf_indices = get_indices(chd_conf_arr, layer=NLAY-1)
chd_conf_df = pd.DataFrame({'index': chd_conf_indices})
# chd_conf_in = influx_df.copy()
# chd_conf_out = outflux_df.copy()

# k
all_k = []
for i in range(NLAY):
    if i == 2:
        k = np.ones((nrow, ncol)) * 0.001  # horizontal hydraulic conductivity in m/day
    elif i < 6:
        k = np.ones((nrow, ncol)) * 10
    elif i == 6:
        k = np.ones((nrow, ncol)) * 0.001
    else:
        k = np.ones((nrow, ncol)) * 1000  # horizontal hydraulic conductivity in m/day
    all_k.append(k)
k_hor = np.array(all_k)  # horizontal hydraulic conductivity

# --------------------------------------------------------------------------
# other model parameters
init_h = np.stack([top] * NLAY)  # initial head, based on average drain elevation

#drn
drain_input['cond'] = RES**2 * 1
drn_top_input['cond'] = 1

# chd
chd_pw_df['head'] = [i[1] for i in chd_pw_indices]  # head for Poukawa boundary
chd_conf_df['head'] = 13  # head for confined boundary inflow
# chd_conf_in['head'] = 13  # head for confined boundary
# chd_conf_out['head'] = 12  # head for confined boundary

# mbr
mbr_df['flux'] = 2 # m3/d
mbr_df = mbr_df[~mbr_df['index'].isin(chd_pw_df['index'])]  # remove chd indices from mbr

# influx/outfluc
influx_df['flux'] = 1  # m3/d
outflux_df['flux'] = -1  # m3/d

# remove where mbr is active
influx_df = influx_df[~influx_df['index'].isin(mbr_df['index'])]  # remove mbr indices from influx
outflux_df = outflux_df[~outflux_df['index'].isin(mbr_df['index'])]  # remove mbr indices from outflux

# recharge
recharge = np.ones_like(idomain.data[0]) * 0.0001 * idomain.data[0]

icell_type = np.zeros((NLAY, nrow, ncol), dtype=int)  # cell type for each layer

# SAVE MODEL PARAMETERS ---------------------------------------------
k_fn = os.path.join(MODEL_DIR, f'{MODEL_NAME}.npf_k_layer')
rch_fn = os.path.join(MODEL_DIR, f'{MODEL_NAME}.rcha_recharge.txt')
drn_fn = os.path.join(MODEL_DIR, f'{MODEL_NAME}.drn_riv_stress_period_data.txt')
chdpw_fn = os.path.join(MODEL_DIR, f'{MODEL_NAME}.chd_pw_stress_period_data.txt')
chdconf_fn = os.path.join(MODEL_DIR, f'{MODEL_NAME}.chd_conf_stress_period_data.txt')
mbr_fn = os.path.join(MODEL_DIR, f'{MODEL_NAME}.wel_mbr_stress_period_data.txt')
influx_fn = os.path.join(MODEL_DIR, f'{MODEL_NAME}.wel_influx_stress_period_data.txt')
outflux_fn = os.path.join(MODEL_DIR, f'{MODEL_NAME}.wel_outflux_stress_period_data.txt')

for i in range(NLAY):
    ilay = i + 1  # layer number starts from 1
    np.savetxt(k_fn + f'{ilay}.txt', k_hor[i])
    # np.savetxt(os.path.join(MODEL_DIR, f'{MODEL_NAME}.npf_icelltype_layer{ilay}.txt'), icell_type[i])  # save cell type for each layer

np.savetxt(rch_fn, recharge)  # save bottom elevation for each layer
savedf2txt(drain_input, filename=drn_fn)
savedf2txt(chd_pw_df, filename=chdpw_fn)
savedf2txt(chd_conf_df, filename=chdconf_fn)
savedf2txt(mbr_df, filename=mbr_fn)
savedf2txt(influx_df, filename=influx_fn)
savedf2txt(outflux_df, filename=outflux_fn)

# 2 BUILD A MODEL -------------------------------------------------------

sim = fp.mf6.MFSimulation(sim_name=MODEL_NAME, # name of simulation
                          version='mf6', # version of MODFLOW
                          exe_name=r'C:\Users\tfo46\e_Python\a_rbm\rapid-gwm-build\examples\bin\mf6.exe', # absolute path to MODFLOW executable
                          sim_ws=MODEL_DIR, # path to workspace where all files are stored
                         )

tdis = fp.mf6.ModflowTdis(simulation=sim, # add to the simulation called sim (defined in prevous code cell)
                          time_units='DAYS', 
                          nper=1, # number of stress periods
                          perioddata=[(1, 1, 1)], # period length, number of steps, timestep multiplier
                         )

ims = fp.mf6.ModflowIms(simulation=sim, 
                        complexity='COMPLEX',
                        print_option='ALL'
                       )

gwf = fp.mf6.ModflowGwf(simulation=sim, 
                        modelname=MODEL_NAME, # model name
                        model_nam_file=f"{MODEL_NAME}.nam", # name of nam file
                        save_flows=True, # make sure all flows are stored in binary output file
                       )

dis = fp.mf6.ModflowGwfdis(model=gwf, # add to groundwater flow model called gwf
                           length_units='METERS', 
                           nlay=NLAY, 
                           nrow=nrow, 
                           ncol=ncol,
                           delr=delr, 
                           delc=delc, 
                           top=top, 
                           botm=botm,
                           idomain=idomain, # 3D array of active cells 
                           xorigin=grid.bounds[0],
                           yorigin=grid.bounds[1],
                          )

npf = fp.mf6.ModflowGwfnpf(model=gwf, #node property flow package
                           save_specific_discharge=True, # save the specific discharge for every cell
                        #    icelltype=icell_type, # cell type for each layer
                           icelltype=0, # 0 means constant saturated thickness
                        #    k = k_hor
                           k=k_hor, # horizontal k value
                          )

ic = fp.mf6.ModflowGwfic(model=gwf, 
                         strt=init_h, # initial head, only used for iterative solution in steady model (arbitrary)
                        )

rch = fp.mf6.ModflowGwfrcha(model=gwf,
                            # recharge=recharge, # recharge array for each cell
                            recharge={0: {'filename': os.path.basename(rch_fn)}}, # recharge for each cell
                            pname='rch' # package name
                           )

drn_riv = fp.mf6.ModflowGwfdrn(model=gwf, # add drain package to model gwf (created in previous code cell)
                            #    stress_period_data={0: drain_input.values.tolist()},
                           stress_period_data={0: {'filename': os.path.basename(drn_fn)}},
                            pname='drn_r', # package name
                            )

chd_pw = fp.mf6.ModflowGwfchd(model=gwf, # add chd package to model gwf (created in previous code cell)
                              stress_period_data={0: chd_pw_df[['index', 'head']].values.tolist()},
                            # stress_period_data={0: {'filename': os.path.basename(chdpw_fn)}},
                            pname='chd_pw', # package name
                            save_flows=True, # save flows for this package 
                           )

chd_conf = fp.mf6.ModflowGwfchd(model=gwf, # add chd package to model gwf (created in previous code cell)
                            # stress_period_data={0: chd_conf_df.values.tolist()},
                            stress_period_data={0: {'filename': os.path.basename(chdconf_fn)}},
                            pname='chd_conf', # package name
                            save_flows=True, # save flows for this package 
                           )

wel = fp.mf6.ModflowGwfwel(model=gwf,
                        #    stress_period_data={0: mbr_df[['index', 'flux']].values.tolist()},
                           stress_period_data={0: {'filename': os.path.basename(mbr_fn)}},
                           pname='mbr' # package name
                          )
influx = fp.mf6.ModflowGwfwel(model=gwf,
                            #   stress_period_data={0: influx_df[['index', 'flux']].values.tolist()},
                           stress_period_data={0: {'filename': os.path.basename(influx_fn)}},
                           pname='influx' # package name
                          )
outflux = fp.mf6.ModflowGwfwel(model=gwf,
                            # stress_period_data={0: outflux_df[['index', 'flux']].values.tolist()},
                           stress_period_data={0: {'filename': os.path.basename(outflux_fn)}},
                           pname='outflux' # package name
                          )


oc = fp.mf6.ModflowGwfoc(model=gwf, # add output control to model gwf (created in previous code cell)
                         budget_filerecord=f"{MODEL_NAME}.cbc", # file name where all budget output is stored
                         head_filerecord=f"{MODEL_NAME}.hds", # file name where all head output is stored
                         saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
                        )

# --------------------------------------------------------

print('Writing model files...')
# sim.set_all_data_external()
sim.write_simulation()  # write all model files to disk
print('Running model...')
success, _ = sim.run_simulation()  # run the model


# TEST OBS ------------------------------------------------------------

helpers.extract_heads_and_budget(ws=MODEL_DIR)  # extract heads and budget from model output
helpers.extract_spring_obs(gwf=gwf, ws=MODEL_DIR)  # extract spring observations from model output

# OUPUT PLOTS --------------------------------------------------------
# visual check
pmv = fp.plot.PlotMapView(model=gwf, layer=0) # create view of layer 0
pmv.plot_array(top *idomain[0], masked_values=[1e30], alpha=0.5, cmap='viridis') # plot top elevation

pmv.plot_bc(name='chd_pw', color='purple') # add 'chd' cells
pmv.plot_bc(name='mbr', color='orange') # add 'wells' cells
pmv.plot_bc(name='drn_r', color='blue') # add 'wells' cells
pmv.plot_bc(name='chd_conf', color='orange') # add 'chd' cells
# pmv.plot_bc(name='chd_conf_out', color='gold') # add 'chd' cells
pmv.plot_bc(name='influx', color='green') # add 'influx' cells
pmv.plot_bc(name='outflux', color='red') # add 'outflux' cells
# pmv.plot_inactive(color='lightgray', alpha=0.5) # plot inactive cells
# pmv.plot_grid(colors='silver', lw=0.01); # add grid

# save to figures
plt.savefig(os.path.join(FIG_DIR, f'{MODEL_NAME}_domain.png'), dpi=300, bbox_inches='tight') # save figure

# --------------------------------------------------------------

# plot heads
idom_plt = np.where(idomain[0] == 1, np.nan, idomain[0])  # create a mask for the active domain
head = gwf.output.head().get_data()
for i in range(NLAY):
    
    pmv = fp.plot.PlotMapView(model=gwf, layer=i)

    pmv.plot_bc(name='chd_pw', color='purple') # add 'chd' cells
    pmv.plot_bc(name='mbr', color='orange') # add 'wells' cells
    pmv.plot_bc(name='drn_r', color='blue') # add 'wells' cells
    pmv.plot_bc(name='chd_conf', color='orange') # add 'chd' cells
    # pmv.plot_bc(name='chd_conf_out', color='gold') # add 'chd' cells
    pmv.plot_bc(name='influx', color='green') # add 'influx' cells
    pmv.plot_bc(name='outflux', color='red') # add 'outflux' cells
    
    # pmv.plot_inactive(color='lightgray', alpha=0.5) # plot inactive cells
    # pmv.plot_array(idom_plt, cmap='gray', alpha=0.2) # plot head array for layer 0
    # pmv.plot_grid(colors='silver', lw=0.1); # add grid
    cs = pmv.contour_array(head[i], linewidths=1, colors='k') # contour plot of heads
    plt.clabel(cs, fmt='%1.1f'); # add contour labels with one decimal place

    # save to figures
    plt.savefig(os.path.join(FIG_DIR, f'{MODEL_NAME}_heads{i}.png'), dpi=300, bbox_inches='tight') # save figure
    plt.close()

# -------------------------------------------------------------
# create a cross section grid & plot heads

cross_col = 15 #which row to show cross section
crossview = fp.plot.PlotCrossSection(model=gwf, line={'column': cross_col})
strtArray = crossview.plot_array(head, masked_values=[1e30], alpha=0.8) # plot the active domain
crossview.plot_grid(colors='black', lw=1); # add grid
crossview.plot_inactive(color='lightgray') # plot inactive cells
cb = plt.colorbar(strtArray, shrink=0.5) # add color bar
# strtArray = crossview.plot_array(head, masked_values=[1e30], alpha = 0.5) # plot the array of heads in cross section
plt.savefig(os.path.join(FIG_DIR, f'{MODEL_NAME}_xsection_col{cross_col}.png'), dpi=300, bbox_inches='tight') # save figure
plt.close()  # close the figure to avoid memory issues




