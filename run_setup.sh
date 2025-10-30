#!/bin/bash
# Bash script to run complete PEST++ IES workflow
# Run from repository root

echo "======================================================================"
echo "RUNNING COMPLETE PEST++ IES WORKFLOW"
echo "======================================================================"

# Activate virtual environment
echo ""
echo "Activating virtual environment..."
source .venv/Scripts/activate

# Run model setup
echo ""
echo "======================================================================"
echo "STEP 1: Running a_model_setup.py"
echo "======================================================================"
python scripts/a_model_setup.py
if [ $? -ne 0 ]; then
    echo "ERROR: a_model_setup.py failed with exit code $?"
    read -p "Press enter to exit..."
    exit $?
fi

# Run parameterization
echo ""
echo "======================================================================"
echo "STEP 2: Running b_parameterization.py"
echo "======================================================================"
python scripts/b_parameterization.py
if [ $? -ne 0 ]; then
    echo "ERROR: b_parameterization.py failed with exit code $?"
    read -p "Press enter to exit..."
    exit $?
fi

# # Run prior
# echo ""
# echo "======================================================================"
# echo "STEP 2: Running c_pest_IES_prior.py"
# echo "======================================================================"
# python scripts/c_pest_IES_prior.py
# if [ $? -ne 0 ]; then
#     echo "ERROR: c_pest_IES_prior.py failed with exit code $?"
#     read -p "Press enter to exit..."
#     exit $?
# fi

# Run PEST++ IES history matching
echo ""
echo "======================================================================"
echo "STEP 3: Running d_pest_IES_HM.py"
echo "======================================================================"
echo "WARNING: This step may take several hours to complete!"
python scripts/d_pest_IES_HM.py
if [ $? -ne 0 ]; then
    echo "ERROR: d_pest_IES_HM.py failed with exit code $?"
    read -p "Press enter to exit..."
    exit $?
fi

# Run post-processing figures
echo ""
echo "======================================================================"
echo "STEP 4: Running g_postprocess_figures.py"
echo "======================================================================"
python scripts/g_postprocess_figures.py
if [ $? -ne 0 ]; then
    echo "ERROR: g_postprocess_figures.py failed with exit code $?"
    read -p "Press enter to exit..."
    exit $?
fi

echo ""
echo "======================================================================"
echo "COMPLETE WORKFLOW FINISHED SUCCESSFULLY!"
echo "======================================================================"
echo ""
echo "Results are available in:"
echo "  - models/{MODEL_NAME}/pest/master_ies/ (PEST++ results)"
echo "  - models/{MODEL_NAME}/figures/ (output figures)"
echo "======================================================================"

read -p "Press enter to exit..."
