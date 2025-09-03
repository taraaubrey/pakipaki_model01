MODEL_NAME = 'local2'  # name of the model
# model domain
RES = 50
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
WETLANDA_PAST = r"data/wetlandA_past.shp"

SPRING = r"data/spring.shp"

POUKAWA_BOUNDARY = r"data/model2_chd.shp"
INFLUX_BOUNDARY = r"data/model2_influx.shp"
OUTFLUX_BOUNDARY = r"data/model2_outflux.shp"

# TS FILES
START_DATE = '2023-12-7'  # start date for the simulation
END_DATE = '2024-1-28'  # end date for the simulation
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
TRUTH_DIR = r'truth'

# Particle locations
SAMPLES = r"data/sample_locations.shp"

# Pest options
NREALS_PRIOR = 100  # number of realizations for parameterization
NREALS = 500
NOPTMAX = 3  # number of optimization iterations