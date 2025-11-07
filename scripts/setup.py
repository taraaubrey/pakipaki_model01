MODEL_NAME = 'local_run28'  # name of the model
"""
- Removed recharge from parameterization.
- Added shapefile for pilot points.
- Updated IES settings for reinflation.
- 200 realizations for IES.
"""
# model domain
RES = 25
NLAY = 1
NLAY_THICKNESS = 10  # thickness of each layer in meters

# TDIS
NPER = 4  # number of stress periods
NSTEPS = 5

#% paths
DOMAIN = r"data/model2_domain.shp"
TOP = r"data/model2_dem.tif"
BOTTOM = r"data/basement_z.tif"
DRAINS = r"data/drains.shp"
SPRING_DRAIN = r"data/spring_drain.shp"
DRAIN_ZONES = r"data/drain_zones.shp"
MBR = r"data/model2_mbr.shp"
LIMESTONE_INACTIVE = r"data/model2_limestone_inactive_bottom.shp"
CONF_AREA_ACTIVE = r"data/confining_area.shp"

#past
PW_PAST = r"data/poukawa_past.shp"
DRAINS_PAST = r"data/drains_past.shp"
WETLANDA_PAST = r"data/wetlandB_past.shp"
WETLAND_INFLUENCE = r"data\wetland_influence_area.shp"

SPRING = r"data/spring.shp"

POUKAWA_BOUNDARY = r"data/model2_chd.shp"
INFLUX_BOUNDARY = r"data/model2_influx.shp"
OUTFLUX_BOUNDARY = r"data/model2_outflux.shp"

# TS FILES
START_DATE = '2023-12-6 00:00'  # start date for the simulation
END_DATE = '2024-1-27 23:59'  # end date for the simulation
AWANUI_TS = r"data/river-level-for-awanui-s.csv" #mm
AWANUI_Q = r"data/river-flow-for-awanui-st.csv" #m3/d
POUKAWA_TS = r"data/river-level-for-poukawa.csv" #mRL
POUKAWA_Q = r"data/river-flow-for-poukawa-s.csv" #m3/d
SPRING_TS = r"data/spring_LEVEL_D_clean.csv" #mRL
PK4_TS = r"data/pk4_LEVEL_D_clean.csv" #mRL
PK2_TS = r"data/pk2_LEVEL_D_clean.csv" #mRL
CONF_TS = r"data/1501_Harrisons.csv" # mabgl

# Directories / paths
BIN_DIR = f'bin/' # relative from model_dir
SCRIPTS_DIR = r'scripts'
FIG_DIR = f'models/{MODEL_NAME}/figures'  # directory for figures
MODEL_DIR = f'models/{MODEL_NAME}/{MODEL_NAME}' # model workspace to be used
SPATIAL_DIR = f'models/{MODEL_NAME}/spatial'  # directory for spatial ./data
PEST_DIR = f'models/{MODEL_NAME}/pest/{MODEL_NAME}'  # directory for pest files
TEMP_DIR = f'models/{MODEL_NAME}/pest/{MODEL_NAME}_template'  # directory for temporary files
HPTEMP_DIR = f'models/{MODEL_NAME}/pest/{MODEL_NAME}_hptemplate'  # directory for temporary files
TRUTH_DIR = r'truth'

# Particle locations
SAMPLES = r"data/sample_locations.shp"

PP_SHP = r"data/pilot_points.shp"
GR_SHP = r"data/pp_grid.shp" # polygon of grid zone


# model parameters
SS_PRIOR = {
    'initial': 1e-1,
    'lb': 1e-4,
    'ub': 1e5,
    'ulb': 1e-6,
    'uub': 1,
}

KH_PRIOR = {
    'initial': 100,
    'lb': 1e-6,
    'ub': 1e8,
    'ulb': 1e-6,
    'uub': 1e4,
}

RCH = {
    'initial_rf': 2.5e-4, # rainfall recharge 0.00025
    'initial_mbr': 1.6e-3, # mbr recharge
    'lb': 7.2e-1,
    'ub': 1.8,
    'ulb': 1.3e-3,
    'uub': 2.4e-3,
}

GHB_SW = {
    'initial_cond_aw': 100,
    'initial_cond_pw': 100,
    'initial_cond_spr': 100,
    'cond_lb': 1e-6,
    'cond_ub': 1e8,
    'cond_ulb': 1e-4,
    'cond_uub': 1e6,
    'head_lb': -0.5,
    'head_ub': 0.5,
    'head_ulb': 5,
    'head_uub': 12,
}

GHB_CONF = {
    'initial_cond': 1,
    'initial_head_offset': 2,
    'initial_head_min': 14,
    'cond_lb': 1e-6,
    'cond_ub': 1e8,
    'cond_ulb': 1e-4,
    'cond_uub': 1e6,
    'head_lb': -0.5,
    'head_ub': 0.5,
    'head_ulb': 10,
    'head_uub': 16,
}

AWANUI_water_offset = 0.5  # m (above the minimum elevation as initial head)
CONF_past_min = 12.95 # m ; represents the lowest water level recorded in confining area in past

# truth
GHB_Q = -6.9  # m3/d
GHB_Qstd = 1.8  # m3/d +/- 6 m3/d (assume 4 std is full range) 12/4

HEAD_offset = 1 # from top of model domain (water level can't be above this)
HEAD_std = 0.025  # m  top heads +/- 0.5 cm (assume 4 std is full range) 1/4
PK4_std = 0.01  # m +/- 2 cm (assume 4 std is full range) 0.04/4

PHI_OBS = {
        # important budget and flux obs
        'arr-h': 60, # head obs
        'budget': 50, # flux obs
        # temporal observations
        'arr-confq': 60, # flux obs
        'arr-awq': 10, # flux obs
        'ts-heads': 100, # head obs
        'recession': 100, # head obs
}

# Pest options
NREALS = 200  # number of realizations for parameterization
NREALS_PRIOR = 50  # number of realizations for parameterization
REINFLATE_ITERS = None  # iterations at which to reinflate

if REINFLATE_ITERS:
    NOPTMAX = sum(REINFLATE_ITERS)  # number of optimization iterations
else:
    NOPTMAX = 3 

PEST_PP_OPTIONS = {
    'ies_num_reals': NREALS,
    'ies_parameter_ensemble': 'prior_pe.jcb',
    'ies_no_noise': False,
    'ies_verbose_level': 1,
    'overdue_giveup_fac': 10,
    'overdue_giveup_minutes': 15,
    'ies_bad_phi_sigma': 2,
    'par_sigma_range': 6,
    # 'ies_save_rescov': True, doesn't scale beyond 20,000 nnzobs
    # 'ies_reinflate_factor': 1,
    # 'ies_n_iter_reinflate': 3,
}
