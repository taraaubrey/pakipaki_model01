#!/bin/bash -l

## SLURM job script for post-processing PEST++ IES results

## Job configuration
#SBATCH --job-name=GW_PP                                # Job name
#SBATCH --cpus-per-task=64                              # Number of logical CPUs (adjust based on your num_workers)
#SBATCH --time=01:00:00                                 # Maximum runtime (adjust based on your model complexity)
#SBATCH --mem-per-cpu=500MB                             # Memory per CPU (adjust based on your model size)
#SBATCH --output=runs/SLURM_output/%j/pp_%j.out    # Standard output file
#SBATCH --error=runs/SLURM_output/%j/pp_%j.err     # Standard error file

# run this script from the project root directory with `sbatch runs/run_postprocess.sh`

# Load required modules
module load SciPy-bundle/2024.05-gfbf-2024a
start_time=$(date +%s)


# Model run description
# Define paths
OLD_RUN_ID="62576"
TAG="run34_rch250"
PROJECT_ROOT="${HOME}/runs/SLURM_output/${OLD_RUN_ID}/${TAG}"
WORK_DIR="${HOME}/runs/SLURM_output/${SLURM_JOB_ID}/${OLD_RUN_ID}_${TAG}"


echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"


mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

echo "Working directory: $WORK_DIR"
rsync -a "$PROJECT_ROOT/" .

echo "#####################################################################"
echo "Post-Processing Script"
echo "Start time: $(date)"
echo "#####################################################################"
echo ""
echo "Note: Using setup.py configuration from current directory"
echo ""

# Define commands to run
CMDS=(
    "python scripts/f_analyze_phi.py"
    "python scripts/g_postprocess.py --subset-method all"
    "python scripts/g_postprocess.py --subset-method phi --phi-percentile 25 75"
    "python scripts/g_postprocess.py --subset-method head --head-threshold 15"
)

# Run each command in sequence
for cmd in "${CMDS[@]}"; do
    echo ""
    echo "#####################################################################"
    echo "Running command: $cmd"
    echo "#####################################################################"
    echo "Start time: $(date)"

    if ! eval "$cmd"; then
        echo ""
        echo "ERROR: Command failed: $cmd"
        echo "Exiting..."
        exit 1
    fi

    echo "End time: $(date)"
    echo ""
done

# Calculate and display total runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
minutes=$((runtime / 60))
seconds=$((runtime % 60))

echo "#####################################################################"
echo "All commands completed successfully"
echo "#####################################################################"
echo "Total runtime: ${minutes}m ${seconds}s ($runtime seconds)"
echo "Job completed at: $(date)"
echo "#####################################################################"

exit 0
