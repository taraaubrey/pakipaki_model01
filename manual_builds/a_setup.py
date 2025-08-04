MODEL_NAME = 'local1'  # name of the model
# model domain
RES = 10
NLAY = 8
NLAY_THICKNESS = 10  # thickness of each layer in meters

#% paths
DOMAIN = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\model2_domain.shp"
TOP = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data/b_derived/dem_elevation_derivatives/dem_clipped.tif"
BOTTOM = r"C:\Users\tfo46\e_Python\a_rbm\rapid-gwm-build\examples\pakipaki\models\derived_data\basement_z.tif"
DRAINS = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\model2_drains.shp"
MBR = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\model2_mbr.shp"
LIMESTONE_INACTIVE = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\model2_limestone_inactive_bottom.shp"
CONF_AREA_ACTIVE = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\confining_area.shp"

SPRING = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\spring.shp"


POUKAWA_BOUNDARY = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\model2_chd.shp"
INFLUX_BOUNDARY = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\model2_influx.shp"
OUTFLUX_BOUNDARY = r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\b_derived\mod_model_files\pakipaki\shp\model2_outflux.shp"

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