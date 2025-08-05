# Running PEST Monte Carlo Analysis on SLURM

This guide explains how to run your `e_pest_MC.py` script on the SLURM cluster.

## Files Created

1. **`run_pest_mc_v2.sh`** - Main SLURM job script (recommended)
2. **`run_pest_mc.sh`** - Alternative job script 
3. **`check_requirements.py`** - Requirements checker script

## Quick Start

### 1. Check Requirements First
```bash
cd /home/tfo46/pakipaki_model01/manual_builds/dvc_scripts
python check_requirements.py
```

This will verify:
- Python packages (pyemu, flopy, pandas, numpy)
- PEST++ executables (pestpp-swp.exe, mf6.exe, mp7.exe)
- Project directory structure
- setup.py configuration

### 2. Fix Any Missing Dependencies
If the check shows missing packages, install them:
```bash
pip install --user pyemu flopy pandas numpy matplotlib
```

### 3. Submit the Job
```bash
sbatch run_pest_mc_v2.sh
```

### 4. Monitor the Job
```bash
# Check job status
squeue -u $USER

# View output (replace JOBID with actual job ID)
tail -f pest_mc_JOBID.out
tail -f pest_mc_JOBID.err
```

## What the Job Script Does

The `run_pest_mc_v2.sh` script:

1. **Sets up the environment**:
   - Loads Python and SciPy modules
   - Allocates 8 CPUs and 2GB RAM per CPU
   - Sets 2-hour time limit

2. **Copies files to compute node**:
   - Copies entire project to local storage for better I/O performance
   - Sets up directory structure to match your setup.py paths
   - Makes PEST++ executables executable

3. **Installs missing dependencies**:
   - Automatically installs pyemu and flopy if not available

4. **Runs your script**:
   - Changes to correct directory
   - Runs `e_pest_MC.py` with proper environment

5. **Copies results back**:
   - Copies generated model files back to original location
   - Copies sweep_in.csv and other outputs
   - Copies log files for debugging

## Your Script Analysis

Your `e_pest_MC.py` script:
- Generates 250 parameter realizations using Gaussian draw
- Uses `pestpp-swp` for sweep analysis
- Uses all available CPU cores (`os.cpu_count()`)
- Creates a master directory for parallel execution

## Expected Outputs

After successful completion, you should see:
- `sweep_in.csv` - Parameter ensemble file
- Model results in the master_mc directory
- PEST++ output files
- Job log files

## Troubleshooting

### Common Issues:

1. **Import errors for pyemu/flopy**:
   - Install with: `pip install --user pyemu flopy`

2. **Executable not found**:
   - Check that executables exist in `/home/tfo46/pakipaki_model01/bin/`
   - Run: `chmod +x /home/tfo46/pakipaki_model01/bin/*`

3. **Path issues**:
   - Your setup.py expects paths starting with "examples/"
   - The job script handles this by creating the proper directory structure

4. **Memory issues**:
   - Increase `--mem-per-cpu` if you get out-of-memory errors
   - Reduce `num_reals` from 250 to a smaller number for testing

5. **Time limit exceeded**:
   - Increase `--time` if your model takes longer than 2 hours

### Debugging Steps:

1. Run the requirements checker first
2. Test with a smaller number of realizations (change `num_reals=250` to `num_reals=10`)
3. Check the job output files for specific error messages
4. Verify your model runs locally first

## Customization

### Adjust Resources:
Edit the SBATCH directives in `run_pest_mc_v2.sh`:
```bash
#SBATCH --cpus-per-task=16         # More CPUs for faster execution
#SBATCH --time=04:00:00            # Longer time limit
#SBATCH --mem-per-cpu=4GB          # More memory per CPU
```

### Change Number of Realizations:
Edit `e_pest_MC.py` and change:
```python
num_reals=250  # Reduce for testing, e.g., num_reals=50
```

### Use Different PEST++ Version:
Edit `e_pest_MC.py` and change:
```python
'pestpp-swp'  # Could be 'pestpp-ies', 'pestpp-glm', etc.
```

## File Locations

- **Scripts**: `/home/tfo46/pakipaki_model01/manual_builds/dvc_scripts/`
- **Executables**: `/home/tfo46/pakipaki_model01/bin/`
- **Job outputs**: Same directory as scripts
- **Model results**: Will be created in the model directory structure

Run `python check_requirements.py` first, then submit with `sbatch run_pest_mc_v2.sh`!
