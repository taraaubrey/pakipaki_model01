MODEL_NAME = 'local2'  # name of the model
# model domain
RES = 10
NLAY = 8
NLAY_THICKNESS = 10  # thickness of each layer in meters

#% paths
DOMAIN = r"examples\manual_builds\data\model2_domain.shp"
TOP = r"examples\manual_builds\data\dem_clipped.tif"
BOTTOM = r"examples\pakipaki\models\derived_data\basement_z.tif"
DRAINS = r"examples\manual_builds\data\model2_drains.shp"
MBR = r"examples\manual_builds\data\model2_mbr.shp"
LIMESTONE_INACTIVE = r"examples\manual_builds\data\model2_limestone_inactive_bottom.shp"
CONF_AREA_ACTIVE = r"examples\manual_builds\data\confining_area.shp"

SPRING = r"examples\manual_builds\data\spring.shp"

POUKAWA_BOUNDARY = r"examples\manual_builds\data\model2_chd.shp"
INFLUX_BOUNDARY = r"examples\manual_builds\data\model2_influx.shp"
OUTFLUX_BOUNDARY = r"examples\manual_builds\data\model2_outflux.shp"

# Directories / paths
BIN_DIR = f'examples/bin'
SCRIPTS_DIR = r'examples/manual_builds/'
FIG_DIR = f'examples/manual_builds/models/{MODEL_NAME}/figures'  # directory for figures
MODEL_DIR = f'examples/manual_builds/models/{MODEL_NAME}/{MODEL_NAME}' # model workspace to be used
SPATIAL_DIR = f'examples/manual_builds/models/{MODEL_NAME}/spatial'  # directory for spatial data
PEST_DIR = f'examples/manual_builds/models/{MODEL_NAME}/pest/{MODEL_NAME}'  # directory for pest files
TEMP_DIR = f'examples/manual_builds/models/{MODEL_NAME}/pest/{MODEL_NAME}_template'  # directory for temporary files
TRUTH_DIR = r'examples\manual_builds\truth'

# Particle locations
SAMPLES = r"examples\manual_builds\data\sample_locations.shp"