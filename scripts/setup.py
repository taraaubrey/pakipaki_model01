MODEL_NAME = 'run44'  # name of the model
"""
- Removed recharge from parameterization.
- Added shapefile for pilot points.
- Updated IES settings for reinflation.
- 200 realizations for IES.
"""
# model domain
RES = 10
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

SHALLOW_K = r"data\20251127_shallow_mean_top_neg4_mRL_filled.tif"
CONF_K = r"data\20251121_confined_mean_neg6_neg14_mRL_clipped.tif"
CONF_THICKNESS = 10  # m

#past
PW_PAST = r"data/poukawa_past.shp"
DRAINS_PAST = r"data/drains_past.shp"
WETLANDA_PAST = r"data/wetlandB_past.shp"
WETLAND_INFLUENCE = r"data\wetland_influence_area.shp"

SPRING = r"data/spring.shp"

POUKAWA_BOUNDARY = r"data/model2_chd.shp"
AWANUI_BOUNDARY = r"data/awanui_boundary.shp"
AWANUI_PAST = r"data/awanui_past.shp"
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
POST_DIR = f'models/{MODEL_NAME}/postprocessing'  # directory for post-processed files
MODEL_DIR = f'models/{MODEL_NAME}/{MODEL_NAME}' # model workspace to be used
SPATIAL_DIR = f'models/{MODEL_NAME}/spatial'  # directory for spatial ./data
PEST_DIR = f'models/{MODEL_NAME}/pest/{MODEL_NAME}'  # directory for pest files
TEMP_DIR = f'models/{MODEL_NAME}/pest/{MODEL_NAME}_template'  # directory for temporary files
HPTEMP_DIR = f'models/{MODEL_NAME}/pest/{MODEL_NAME}_hptemplate'  # directory for temporary files
TRUTH_DIR = f'models/{MODEL_NAME}/truth'
# Particle locations
SAMPLES = r"data/sample_locations.shp"

PP_SHP = r"data/pilot_points.shp"
GR_SHP = r"data/pp_grid.shp" # polygon of grid zone

# PARAMS
RIV_W = 2
RIV_L = 2
PAST_f = 2e-1 # (factor on conductance to make much lower)
SP_W = 2

SP_f = 1e-1 # (factor on conductance to make much lower)
AW_f = 1e0
PW_f = 1e0
RCH_f = 0.2 # (anything lower than e-4; little impact on recharge)
KH_f = 10 # multiplier for kh

AW_WL = 7.22  # mRL
PW_WL = 7.66  # mRL
SP_WL = 7.59  # mRL
PAST_OFFSET = 1.0 # m

AWANUI_water_offset = 0.5  # m (above the minimum elevation as initial head)
CONF_past_min = 12.95 # m ; represents the lowest water level recorded in confining area in past
HEAD_offset = 1 # from top of model domain (water level can't be above this)
HEAD_std = 0.025  # m
PK4_std = 0.01  # m +/- 2 cm (assume 4 std is full range) 0.04/4

# truth ++++++++++++++++++++++++++++++++++++++
GHB_Q = -2.8  # m3/d
GHB_Qstd = 3 # std
GHB_SPRING_Q = -1.4  # m3/d
GHB_SPRING_Qstd = 0.46 # std
# budget
SW_TOTAL = -308.2  # m3/d (must be less than)
SW_STD = 81.6 # m3/d
RCH1_TOTAL = 150 # m3/d
RCH2_TOTAL = 60 # m3/d
RCH1_STD = 50  # m3/d
RCH2_STD = 30  # m3/d

# model parameters
SS_PRIOR = {
    'iconvert': True,
    'initial': 1e-4,
    'lb': 1e-2,
    'ub': 1e2,
    'ulb': 1e-6,
    'uub': 1,
}

SY_PRIOR = {
    'initial': 0.3,
    'lb': 1e-1,
    'ub': 1e1,
    'ulb': 1e-2,
    'uub': 1,
}

KH_PRIOR = {
    'lb': 1e-2,
    'ub': 1e2,
    'ulb': 1e-8,
    'uub': 1e8,
}

RCH = {
    'initial_rf': 2.5e-4, # total recharge
    'initial_mbr': 1.57e-3, # mbr recharge
    'lb': 1e-1,
    'ub': 1e1,
    'ulb': 1e-6,
    'uub': 1e6,
}

RCH_PARMS = {
    'lb': 1e-3,
    'ub': 1e3,
    'ulb': 1e-6,
    'uub': 1e6,
}

GHB_SW = {
    # 'initial_cond_aw': 100,
    # 'initial_cond_pw': 100,
    # 'initial_cond_spr': 1000,
    'cond_lb': 1e-2,
    'cond_ub': 1e2,
    'cond_ulb': 1e-6,
    'cond_uub': 1e6,
}

PHI_OBS = {
        'top-heads': 3,
        'budget-rch': 1,
        'budget-sw': 1,
        'spring-q': 1,
        'aw-q': 1,
        'pk4-aw-diff': 1,
        'pk4-spring-diff': 3,
        'pk4-diff': 1,
        'pk4-heads': 3,
}

# parameterization settings (obs)
TIME_SUBSAMPLE = 10  # time subsampling for pilot points
SPACE_SUBSAMPLE = 5  # spatial subsampling for pilot points

# Pest options
NREALS = 250  # number of realizations for parameterization
NREALS_PRIOR = 50  # number of realizations for parameterization
REINFLATE_ITERS = [3] # iterations at which to reinflate
NOPTMAX = 0
MAXSING = None

if REINFLATE_ITERS:
    NOPTMAX = sum(REINFLATE_ITERS) + len(REINFLATE_ITERS) - 1 # number of optimization iterations
else:
    REINFLATE_ITERS = 0

PEST_PP_OPTIONS = {
    'ies_num_reals': NREALS,
    'ies_parameter_ensemble': 'prior_pe.jcb',
    'ies_no_noise': False,
    'ies_verbose_level': 2,
    'overdue_giveup_fac': 10,
    'overdue_giveup_minutes': 15,
    'ies_bad_phi_sigma': 2.5,
    'par_sigma_range': 6,
    'ies_reinflate_factor': 1,
    'ies_n_iter_reinflate': REINFLATE_ITERS,
}


#
CYCLES = 10  # number of cycles
MIN = 0.1  # minimum std for obs in reinflation
MAX = 100.0  # maximum std for obs in reinflation