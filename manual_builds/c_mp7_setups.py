import flopy as fp
import geopandas as gpd
import matplotlib.pyplot as plt

from a_setup import *
from mp7_utils import *

# open the model workspace
sim = fp.mf6.MFSimulation.load(sim_ws=MODEL_DIR)
gwf = sim.get_model(MODEL_NAME)

# open spring 
sample_locations = gpd.read_file(SAMPLES)
sample_locations.set_index('obsnme', inplace=True)

particles = cylinder_partloc(gwf, sample_locations, sim_ws=MODEL_DIR, radius=10, height=0.5, npartob=100, write_partcoord=True)

# Specify the paricle data
particledata = fp.modpath.ParticleData(
    partlocs=particles[['k', 'i', 'j']].values.tolist(),
    structured=True,
    localx=particles['ri'].values.tolist(),
    localy=particles['rj'].values.tolist(),
    localz=particles['rk'].values.tolist(),
    drape=0,
    )
# Group partcles in a group
pg = fp.modpath.ParticleGroup(particledata=particledata)

# Create a modpath model and call it mp
mp = fp.modpath.Modpath7(modelname=MODEL_NAME, # name of the model
                         model_ws=MODEL_DIR, # path to workspace where all files are stored
                         flowmodel=gwf, # groundwater flow model to get the flow from
                         exe_name=r'C:\Users\tfo46\e_Python\a_rbm\rapid-gwm-build\examples\bin\mp7.exe', # absolute path to MODPATH7 executable'
                        )

# Add the Basic package to the model called mp and specify the porosity
mpbas = fp.modpath.Modpath7Bas(model=mp, porosity=0.3)

# Add the MODPATH simulation package
mpsim = fp.modpath.Modpath7Sim(model=mp, # add to model called mp
                               particlegroups=pg, # particle group pg defined above
                               simulationtype='endpoint',
                               stoptimeoption='extend',
                               trackingdirection='backward',
                              )

mp.write_input()
mp.run_model(silent=False)

# fname = os.path.join(MODEL_DIR, f"{MODEL_NAME}.mppth")
# fname_end = os.path.join(MODEL_DIR, f"{MODEL_NAME}.mpend")
# plf = fp.utils.PathlineFile(fname)
# endpoints = fp.utils.EndpointFile(fname_end)

# xorigin = gwf.dis.xorigin.get_data()
# yorigin = gwf.dis.yorigin.get_data()

# drn_cells = gwf.drn.stress_period_data.get_data()[0]['cellid'].tolist()
# mbr_cells =gwf.wel[0].stress_period_data.get_data()[0]['cellid'].tolist()
# pw_cells = gwf.chd[0].stress_period_data.get_data()[0]['cellid'].tolist()
# conf_cells = gwf.chd[1].stress_period_data.get_data()[0]['cellid'].tolist()

# endpoints = fp.utils.EndpointFile(fname_end).get_alldata()

# pclasses = []
# for j in range(len(particles)):
#     end_pt = endpoints[j]
#     end_cell = gwf.modelgrid.intersect(
#             x=end_pt['x'],
#             y=end_pt['y'],
#             z=end_pt['z'],
#             local=True,
#             forgive=True)
        
#     if end_cell in drn_cells:
#         pclass = 'drain'
#     elif end_cell in mbr_cells:
#         pclass = 'mbr'
#     elif end_cell in pw_cells:
#         pclass = 'poukawa boundary'
#     elif end_cell in conf_cells:
#         pclass = 'confining_area'
#     elif end_cell[0] == 0:
#         pclass = 'rainfall_recharge'
#     else:
#         pclass = 'unknown'

#     particles.at[j, 'class'] = pclass
#     particles.at[j, 'x_end'] = end_pt['x'] + xorigin
#     particles.at[j, 'y_end'] = end_pt['y'] + yorigin
#     particles.at[j, 'z_end'] = end_pt['z']
#     particles.at[j, 'i_end'] = end_cell[1]
#     particles.at[j, 'j_end'] = end_cell[2]
#     particles.at[j, 'k_end'] = end_cell[0]
#     particles.at[j, 'cellface'] = end_pt['cellface']
#     particles.at[j, 'status'] = end_pt['status']

# # calculate percent of particles in each class per site
# for site in particles['site'].unique():
#     site_particles = particles[particles['site'] == site]
#     total_particles = len(site_particles)
#     if total_particles > 0:
#         class_counts = site_particles['class'].value_counts(normalize=True) * 100
#         print(f"Site: {site}")
#         for pclass, percent in class_counts.items():
#             print(f"  {pclass}: {percent:.2f}%")
#     else:
#         print(f"Site: {site} has no particles.")


# # PLOT -----------------------------------------------------

# pmv = fp.plot.PlotMapView(model=gwf, layer=0)

# pmv.plot_bc(name='chd_pw', color='purple') # add 'chd' cells
# pmv.plot_bc(name='mbr', color='orange') # add 'wells' cells
# pmv.plot_bc(name='drn_r', color='blue') # add 'wells' cells
# pmv.plot_bc(name='chd_conf', color='orange') # add 'chd' cells
# pmv.plot_bc(name='influx', color='green') # add 'influx' cells
# pmv.plot_bc(name='outflux', color='red') # add 'outflux' cells

# # pmv.contour_array(head[i], levels=np.arange(15, 30, 0.5), linewidths=1, colors='k')
# class_colors = {
#     'rainfall_recharge': 'green',
#     'drain': 'blue',
#     'mbr': 'orange',
#     'poukawa boundary': 'purple',
#     'confining_area': 'black',
#     'unknown': 'gray'
# }
# pclasses = []
# particles['class'] = None  # default class for particles
# for j in range(len(particles)):
#     pline = plf.get_data(partid=j)
#     if len(pline) > 0:
#         lx = pline['x'][-1] + xorigin
#         ly = pline['y'][-1] + yorigin
#         lz = pline['z'][-1]
#         last_cell = gwf.modelgrid.intersect(
#             x=pline['x'][-1],
#             y=pline['y'][-1],
#             z=pline['z'][-1],
#             local=True,
#             forgive=True)
        
#         if last_cell in drn_cells:
#             pclass = 'drain'
#         elif last_cell in mbr_cells:
#             pclass = 'mbr'
#         elif last_cell in pw_cells:
#             pclass = 'poukawa boundary'
#         elif last_cell in conf_cells:
#             pclass = 'confining_area'
#         elif last_cell[0] == 0:
#             pclass = 'rainfall_recharge'
#         else:
#             pclass = 'unknown'

#         pclasses.append(pclass)
#         particles.at[j, 'class'] = pclass

#     plt.plot(pline['x'] + xorigin, pline['y'] + yorigin, class_colors[pclass], lw=0.5, alpha=0.5)  # plot the particle path for layer i

# # save to figures
# plt.title(f'Particle paths')
# plt.savefig(os.path.join(FIG_DIR, f'{MODEL_NAME}_particles.png'), dpi=300, bbox_inches='tight') # save figure
# plt.close()  # close the figure to avoid memory issues


# print('here')