#!/bin/bash -l

## SLURM job script for PEST Monte Carlo analysis using pyemu
## This script runs e_pest_MC.py with all required dependencies

## Job configuration
#SBATCH --job-name=pest_mc_analysis
#SBATCH --cpus-per-task=16          # Number of logical CPUs (adjust based on your num_workers)
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
SCRIPT_DIR="$PROJECT_ROOT/manual_builds/scripts/dvc_scripts"
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
cp -r "$PROJECT_ROOT" .

# Change to the script directory
cd "$WORK_DIR/manual_builds/scripts/dvc_scripts"

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

echo "Job completed at: $(date)"
exit $exit_code
