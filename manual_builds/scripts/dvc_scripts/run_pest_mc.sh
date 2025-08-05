#!/bin/bash -l

## SLURM job script for PEST Monte Carlo analysis using pyemu
## This script runs e_pest_MC.py with all required dependencies

## Job configuration
#SBATCH --job-name=pest_mc_analysis
#SBATCH --cpus-per-task=8          # Number of logical CPUs (adjust based on your num_workers)
#SBATCH --time=02:00:00            # Maximum runtime (adjust based on your model complexity)
#SBATCH --mem-per-cpu=2GB          # Memory per CPU (adjust based on your model size)
#SBATCH --output=pest_mc_%j.out    # Standard output file
#SBATCH --error=pest_mc_%j.err     # Standard error file

# Load required modules
module load SciPy-bundle/2024.05-gfbf-2024a

echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Define paths
PROJECT_ROOT="/home/tfo46/pakipaki_model01"
SCRIPT_DIR="$PROJECT_ROOT/manual_builds/dvc_scripts"
BIN_DIR="$PROJECT_ROOT/bin"
PYTHON_SCRIPT="e_pest_MC.py"


# Create working directory on compute node (for better I/O performance)
WORK_DIR="${HOME}/pest_mc_${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR"

echo "Working directory: $WORK_DIR"
echo "Project root: $PROJECT_ROOT"

# Change to working directory
cd "$WORK_DIR"

# Copy the entire project structure to maintain relative paths
echo "Copying project files..."
cp -r "$PROJECT_ROOT/manual_builds" .
cp -r "$PROJECT_ROOT/bin" .



# Set up Python path to find local modules
# create a for loop of req files and use bash

export PYTHONPATH="$WORK_DIR/manual_builds/dvc_scripts:$PYTHONPATH"

# Add PEST++ executables to PATH
export PATH="$WORK_DIR/bin:$PATH"

# Make all executables in bin directory executable
echo "Setting executable permissions..."
chmod +x ./bin/*

# Verify critical executables are present and executable
echo "Checking for required executables:"
for exe in pestpp-swp mf6 mp7; do
    if [ -x "./bin/${exe}.exe" ]; then
        echo "  ✓ ${exe}.exe found and executable"
    else
        echo "  ✗ ${exe}.exe not found or not executable"
    fi
done

# Check if Python dependencies are available
echo "Checking Python dependencies:"
python -c "import pyemu; print('✓ pyemu available')" || echo "✗ pyemu not available"
python -c "import flopy; print('✓ flopy available')" || echo "✗ flopy not available"
python -c "import pandas; print('✓ pandas available')" || echo "✗ pandas not available"
python -c "import numpy; print('✓ numpy available')" || echo "✗ numpy not available"

# Change to the script directory
cd "$WORK_DIR/manual_builds/dvc_scripts"

# List directory contents to verify setup
echo "Contents of script directory:"
ls -la

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
    echo "PEST Monte Carlo analysis failed"
    exit_code=$?
fi

# Copy results back to original location
echo "Copying results back to project directory..."

# Copy any output files back (adjust these paths based on your actual output)
if [ -d "../models" ]; then
    echo "Copying model results..."
    cp -r ../models "$PROJECT_ROOT/manual_builds/"
fi

# Copy any generated PEST files
if [ -f "sweep_in.csv" ]; then
    cp sweep_in.csv "$SCRIPT_DIR/"
fi

# Copy any log files or other outputs
cp -f *.log "$SCRIPT_DIR/" 2>/dev/null || true

echo "Job completed at: $(date)"
exit $exit_code
