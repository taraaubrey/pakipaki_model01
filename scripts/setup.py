MODEL_NAME = 'local2'  # name of the model
# model domain
RES = 10
NLAY = 8
NLAY_THICKNESS = 10  # thickness of each layer in meters

#% paths
DOMAIN = r"../data/model2_domain.shp"
TOP = r"../data/model2_dem.tif"
BOTTOM = r"../data/basement_z.tif"
DRAINS = r"../data/model2_drains.shp"
MBR = r"../data/model2_mbr.shp"
LIMESTONE_INACTIVE = r"../data/model2_limestone_inactive_bottom.shp"
CONF_AREA_ACTIVE = r"../data/confining_area.shp"

SPRING = r"../data/spring.shp"

POUKAWA_BOUNDARY = r"../data/model2_chd.shp"
INFLUX_BOUNDARY = r"../data/model2_influx.shp"
OUTFLUX_BOUNDARY = r"../data/model2_outflux.shp"


# Directories / paths
BIN_DIR = f'../bin/linux' # relative from model_dir
SCRIPTS_DIR = r'./scripts'
FIG_DIR = f'../models/{MODEL_NAME}/figures'  # directory for figures
MODEL_DIR = f'../models/{MODEL_NAME}/{MODEL_NAME}' # model workspace to be used
SPATIAL_DIR = f'../models/{MODEL_NAME}/spatial'  # directory for spatial ./data
PEST_DIR = f'../models/{MODEL_NAME}/pest/{MODEL_NAME}'  # directory for pest files
TEMP_DIR = f'../models/{MODEL_NAME}/pest/{MODEL_NAME}_template'  # directory for temporary files
TRUTH_DIR = r'../truth'

# Particle locations
SAMPLES = r"../data/sample_locations.shp"