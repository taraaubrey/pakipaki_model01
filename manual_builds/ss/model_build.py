import os
import numpy as np
import pandas as pd
import gridit as gi

from utils import *

def build_model(FIG_DIR):
    # model domain
    RES = 100
    NLAY = 3
    NLAY_THICKNESS = 10  # thickness of each layer in meters

    #% paths
    DOMAIN = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\dom_domains\model_domains\model_domain_20250526.shp"
    TOP = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data/b_derived/dem_elevation_derivatives/dem_clipped.tif"
    BOTTOM = r"C:\Users\tfo46\e_Python\a_rbm\rapid-gwm-build\examples\pakipaki\models\derived_data\basement_z.tif"
    DRAINS = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\a_source\sw_surface_water\HBRC_drains\HBRC_Drains.shp"
    MBR = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\20250618_mbr.shp"
    
    # create directories if they do not exist
    if not os.path.exists(FIG_DIR):
        os.makedirs(FIG_DIR)

    # model data
    model_data = {}

    # open domain
    grid = gi.Grid.from_vector(DOMAIN, RES)
    arr = grid.array_from_vector(DOMAIN)
    arr = arr[np.newaxis, :, :]  # shape (1, nrow, ncol)
    idomain = np.broadcast_to(arr.data, (NLAY, grid.shape[0], grid.shape[1]))
    nrow = grid.shape[0]
    ncol = grid.shape[1]
    delr = np.ones(ncol) * RES
    delc = np.ones(nrow) * RES

    # top
    top = grid.array_from_raster(TOP)
    bottom = grid.array_from_raster(BOTTOM)  # get the bottom elevation
    thickness = top - bottom  # calculate the thickness of the layers

    # bottoms
    min_b = 2
    max_b = 50

    thickness0 = np.where(top - bottom < min_b, min_b, 10)  # set values larger than max_b to max_b
    bot0 = top.data - thickness0
    thickness1 = np.where((bot0 - bottom) < min_b, min_b, 15)  # set values larger than max_b to max_b
    bot1 = bot0 - thickness1
    thickness2 = np.where((bot1 - bottom) < min_b, min_b, 20)  # set values larger than max_b to max_b
    bot2 = bot1 - thickness2

    botm = np.array([bot0, bot1, bot2])

    # drains
    drain_arr = grid.array_from_vector(DRAINS)
    drain_elev = grid.array_from_raster(TOP, resampling='min') * drain_arr.data * idomain[0]  # use the top elevation for drains, only where idomain is active
    drain_cond = 0.9
    drain_input = extract_value_with_indices(drain_elev, layer=0, val_col='elev', mask_value=0)  # extract non-NaN values from the drain elevation array
    drain_input['cond'] = drain_cond

    # add top drains
    drain_top = top * np.where(drain_elev > 0, 0, 1) * idomain[0]
    drn_top_input = extract_value_with_indices(drain_top, layer=0, val_col='elev', mask_value=0)  # extract non-NaN values
    drn_top_input['cond'] = 1

    # recharge
    recharge = 0.0001

    # k
    k_0 = np.ones((nrow, ncol)) * 1000  # horizontal hydraulic conductivity in m/day
    k_1 = np.ones((nrow, ncol)) * 10  # vertical hydraulic conductivity in m/day
    k_2 = np.ones((nrow, ncol)) * 10000  # horizontal hydraulic conductivity in m/day
    k_hor = np.array([k_0, k_1, k_2])  # horizontal hydraulic conductivity

    # mbr
    mbr_arr = grid.array_from_vector(MBR)
    out_mbrarr, mbr_i = get_exterior_indices(mbr_arr, layer=0)  # get exterior indices of the mbr area
    in_idomarr, idom_i = get_interior_indices(idomain[0], layer=0)  # idomainget interior indices of the mbr area
    # return where the mbr is active (1) and where the idomain is active (1)
    mbr_active = np.logical_and(out_mbrarr, in_idomarr)  # mbr area that is active in the model domain
    mbr_indices = get_indices(mbr_active, layer=0)
    mbr_df = pd.DataFrame({'index': mbr_indices})
    mbr_df['flux'] = 10 # m3/d

    # ghb
    # get where in_idomarr is active (1) and where the mbr area is not active (1)
    ghb_active = np.logical_and(in_idomarr, ~mbr_active)  # mbr area that is active in the model domain
    ghb_indices = []
    for i in range(NLAY):
        ghb_i = get_indices(ghb_active, layer=i)
        ghb_indices.extend(ghb_i)
    ghb_df = pd.DataFrame({'index': ghb_indices})
    ghb_df['bhead'] = 40
    ghb_df['cond'] = 0.8  # conductance for GHB

    # correct drns where not ghb
    drain_input = drain_input[~drain_input['index'].isin(ghb_df['index'])]  # remove ghb indices from drain input

    # other model parameters
    init_h = drain_input['elev'].mean()  # initial head, based on average drain elevation

    return model_data