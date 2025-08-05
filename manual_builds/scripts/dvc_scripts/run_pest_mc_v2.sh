#!/bin/bash -l

## SLURM job script for PEST Monte Carlo analysis using pyemu
## This script runs e_pest_MC.py with all required dependencies and proper path handling

## Job configuration
#SBATCH --job-name=pest_mc_analysis
#SBATCH --cpus-per-task=8          # Number of logical CPUs (adjust based on your num_workers)
#SBATCH --time=02:00:00            # Maximum runtime (adjust based on your model complexity)
#SBATCH --mem-per-cpu=2GB          # Memory per CPU (adjust based on your model size)
#SBATCH --output=pest_mc_%j.out    # Standard output file
#SBATCH --error=pest_mc_%j.err     # Standard error file

# Load required modules
module load Python SciPy-bundle

echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Define paths
PROJECT_ROOT="/home/tfo46/pakipaki_model01"
SCRIPT_DIR="$PROJECT_ROOT/manual_builds/dvc_scripts"
PYTHON_SCRIPT="e_pest_MC.py"

# Create working directory on compute node (for better I/O performance)
WORK_DIR="${SLURM_TMPDIR:-/tmp}/pest_mc_${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"

echo "Working directory: $WORK_DIR"
echo "Project root: $PROJECT_ROOT"

# Change to working directory
cd "$WORK_DIR"

# Copy the entire project to maintain the expected directory structure
echo "Copying project files..."
cp -r "$PROJECT_ROOT" ./examples/
# Note: This creates ./examples/pakipaki_model01/ structure that matches setup.py paths

# Set up Python path to find local modules
export PYTHONPATH="$WORK_DIR/examples/pakipaki_model01/manual_builds/dvc_scripts:$PYTHONPATH"

# Add PEST++ executables to PATH
export PATH="$WORK_DIR/examples/pakipaki_model01/bin:$PATH"

# Make all executables in bin directory executable
echo "Setting executable permissions..."
chmod +x ./examples/pakipaki_model01/bin/*

# Create symbolic links to match the expected "examples/bin" path structure
mkdir -p examples/bin
ln -sf "$WORK_DIR/examples/pakipaki_model01/bin/"* examples/bin/

# Verify critical executables are present and executable
echo "Checking for required executables:"
for exe in pestpp-swp mf6 mp7; do
    if [ -x "examples/bin/${exe}.exe" ]; then
        echo "  ✓ ${exe}.exe found and executable"
    else
        echo "  ✗ ${exe}.exe not found or not executable"
    fi
done

# Check if Python dependencies are available
echo "Checking Python dependencies:"
python -c "import pyemu; print('✓ pyemu available')" || echo "✗ pyemu not available - you may need to install it"
python -c "import flopy; print('✓ flopy available')" || echo "✗ flopy not available - you may need to install it"
python -c "import pandas; print('✓ pandas available')" || echo "✗ pandas not available"
python -c "import numpy; print('✓ numpy available')" || echo "✗ numpy not available"

# Check if we need to install pyemu (if not available)
if ! python -c "import pyemu" 2>/dev/null; then
    echo "Installing pyemu..."
    pip install --user pyemu
fi

if ! python -c "import flopy" 2>/dev/null; then
    echo "Installing flopy..."
    pip install --user flopy
fi

# Change to the script directory within the copied structure
cd "$WORK_DIR/examples/pakipaki_model01/manual_builds/dvc_scripts"

# List directory contents to verify setup
echo "Contents of script directory:"
ls -la

# Verify the setup.py paths make sense from current directory
echo "Verifying path structure..."
echo "Current directory: $(pwd)"
echo "Expected bin path: $(python -c "from setup import BIN_DIR; print(BIN_DIR)")"
echo "Expected temp dir: $(python -c "from setup import TEMP_DIR; print(TEMP_DIR)")"

# Set environment variables for the job
export OMP_NUM_THREADS=1  # Prevent oversubscription
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Run the PEST Monte Carlo analysis
echo "Starting PEST Monte Carlo analysis..."
echo "Command: python $PYTHON_SCRIPT"

# Run with error handling
if python "$PYTHON_SCRIPT"; then
    echo "PEST Monte Carlo analysis completed successfully"
    exit_code=0
else
    echo "PEST Monte Carlo analysis failed with exit code $?"
    exit_code=$?
fi

# Copy results back to original location
echo "Copying results back to project directory..."

# The script creates files in the examples/ structure, copy them back
echo "Copying generated model files back..."
if [ -d "../../models" ]; then
    echo "Copying model results..."
    cp -r ../../models "$PROJECT_ROOT/manual_builds/" 2>/dev/null || true
fi

# Copy the specific files mentioned in the script
if [ -f "../../models/local2/pest/local2_template/sweep_in.csv" ]; then
    echo "Copying sweep_in.csv..."
    cp "../../models/local2/pest/local2_template/sweep_in.csv" "$SCRIPT_DIR/" 2>/dev/null || true
fi

# Copy any master_mc directory
if [ -d "../master_mc" ]; then
    echo "Copying master_mc directory..."
    cp -r ../master_mc "$PROJECT_ROOT/manual_builds/" 2>/dev/null || true
fi

# Copy any log files
echo "Copying log files..."
find . -name "*.log" -exec cp {} "$SCRIPT_DIR/" \; 2>/dev/null || true

echo "Job completed at: $(date)"
echo "Exit code: $exit_code"
exit $exit_code
