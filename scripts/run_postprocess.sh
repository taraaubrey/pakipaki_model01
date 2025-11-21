#!/bin/bash
#
# Post-processing script for PEST++ IES results
# Runs phi analysis followed by post-processing with three different subset methods
#
# Usage: bash scripts/run_postprocess.sh
#
# Requirements:
#   - setup.py must be configured with correct MODEL_NAME and paths
#   - PEST++ IES run must be completed
#

# Track start time
start_time=$(date +%s)

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
