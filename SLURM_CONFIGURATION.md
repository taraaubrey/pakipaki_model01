# SLURM Configuration for PEST++ IES on University HPC

**Analysis Date**: 2025-11-18
**Based on**: local_run33 performance data

---

## 📊 Performance Requirements Summary

### Current Baseline (Local Machine - 8 workers)
- **Total runtime**: 565.85 minutes (9.43 hours)
- **Model runs**: 3,548 total
- **Single MODFLOW run**: ~6.7 seconds
- **Memory per model**: ~3 MB (MODFLOW) + ~50 MB (Python) = **~100 MB active**

### Scaling Calculation

**Models per worker per minute**:
- 3,548 models / 565.85 min / 8 workers = **0.78 models/worker/min**

**For different target runtimes**:

| Target Runtime | Required Workers | Speedup |
|----------------|------------------|---------|
| 2 hours | 38 workers | 4.7× |
| 1 hour | 76 workers | 9.4× |
| 30 min | 152 workers | 18.9× |
| 20 min | 228 workers | 28.3× |
| 10 min | 456 workers | 56.6× |
| 5 min | 912 workers | 113× |

---

## 🎯 Recommended SLURM Configurations

### Option 1: Balanced (30-Minute Runtime) - **RECOMMENDED**

```bash
#!/bin/bash
#SBATCH --job-name=pest_ies_run33
#SBATCH --cpus-per-task=152         # 152 workers for 30-min runtime
#SBATCH --time=01:00:00             # 1 hour (30 min run + 30 min buffer)
#SBATCH --mem-per-cpu=500MB         # 500 MB per worker (safe overhead)
#SBATCH --output=runs/SLURM_output/%j/pest_ies_%j.out
#SBATCH --error=runs/SLURM_output/%j/pest_ies_%j.err
#SBATCH --nodes=1                   # Single node (shared memory easier)

# Load required modules (adjust to your HPC environment)
module load python/3.10
module load gcc/11.2.0  # For compiled dependencies

# Activate virtual environment
source /path/to/your/venv/bin/activate

# Set number of PEST++ workers (should match cpus-per-task)
export PEST_NWORKERS=152

# Run PEST++ IES
cd /path/to/pakipaki_01
python scripts/e_pest_IES_HM.py

# Or run directly with pestpp-ies if configured
# cd models/local_run33/pest/master_ies
# pestpp-ies local_run33.pst /n 152
```

**Why this configuration**:
- ✅ Reasonable runtime (30 min vs 9.4 hours local)
- ✅ Moderate resource request (likely to be scheduled quickly)
- ✅ Safe memory allocation (500 MB × 152 = 76 GB total)
- ✅ Fits on most single HPC nodes (up to 192 cores typical)

---

### Option 2: Fast (10-Minute Runtime)

```bash
#!/bin/bash
#SBATCH --job-name=pest_ies_run33_fast
#SBATCH --cpus-per-task=456         # 456 workers for 10-min runtime
#SBATCH --time=00:30:00             # 30 min (10 min run + 20 min buffer)
#SBATCH --mem-per-cpu=500MB         # 500 MB per worker
#SBATCH --output=runs/SLURM_output/%j/pest_ies_%j.out
#SBATCH --error=runs/SLURM_output/%j/pest_ies_%j.err
#SBATCH --nodes=4                   # Likely need multiple nodes
#SBATCH --ntasks-per-node=114       # 456 / 4 = 114 per node

# Load modules
module load python/3.10
module load openmpi/4.1.0  # For multi-node PEST++

source /path/to/your/venv/bin/activate

export PEST_NWORKERS=456

cd /path/to/pakipaki_01
python scripts/e_pest_IES_HM.py
```

**Why this configuration**:
- ✅ Very fast runtime (10 minutes)
- ⚠️ Requires multi-node setup (456 cores)
- ⚠️ May have longer queue time
- ⚠️ Requires MPI or panther networking between nodes

---

### Option 3: Ultra-Fast (5-Minute Runtime)

```bash
#!/bin/bash
#SBATCH --job-name=pest_ies_run33_ultra
#SBATCH --cpus-per-task=912         # 912 workers for 5-min runtime
#SBATCH --time=00:20:00             # 20 min (5 min run + 15 min buffer)
#SBATCH --mem-per-cpu=500MB         # 500 MB per worker
#SBATCH --output=runs/SLURM_output/%j/pest_ies_%j.out
#SBATCH --error=runs/SLURM_output/%j/pest_ies_%j.err
#SBATCH --nodes=8                   # 912 / 114 ≈ 8 nodes
#SBATCH --ntasks-per-node=114

module load python/3.10
module load openmpi/4.1.0

source /path/to/your/venv/bin/activate

export PEST_NWORKERS=912

cd /path/to/pakipaki_01
python scripts/e_pest_IES_HM.py
```

**Why this configuration**:
- ✅ Extremely fast (5 minutes)
- ⚠️ Heavy resource request (may be hard to schedule)
- ⚠️ Diminishing returns (fixed overheads become significant)
- ⚠️ Complex multi-node coordination

---

### Option 4: Conservative (1-Hour Runtime)

```bash
#!/bin/bash
#SBATCH --job-name=pest_ies_run33_conservative
#SBATCH --cpus-per-task=76          # 76 workers for 1-hour runtime
#SBATCH --time=02:00:00             # 2 hours (1 hour run + 1 hour buffer)
#SBATCH --mem-per-cpu=500MB         # 500 MB per worker
#SBATCH --output=runs/SLURM_output/%j/pest_ies_%j.out
#SBATCH --error=runs/SLURM_output/%j/pest_ies_%j.err
#SBATCH --nodes=1

module load python/3.10

source /path/to/your/venv/bin/activate

export PEST_NWORKERS=76

cd /path/to/pakipaki_01
python scripts/e_pest_IES_HM.py
```

**Why this configuration**:
- ✅ Light resource request (quick scheduling)
- ✅ Fits easily on single node
- ✅ Still 9× faster than local
- ⚠️ Longer runtime

---

## 🔧 Adjusting Your Original Config

Your original request was:
```bash
#SBATCH --cpus-per-task=250
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2GB
```

**Analysis**:
- 250 workers → Expected runtime: **~18 minutes** (+ overhead)
- 2 GB/cpu is **EXCESSIVE** (model only needs ~100 MB)
- 2 hour time limit is good (safe buffer)

### Optimized Version:

```bash
#!/bin/bash
#SBATCH --job-name=pest_ies_run33
#SBATCH --cpus-per-task=250         # For ~18 min runtime
#SBATCH --time=01:00:00             # Reduced to 1 hour (plenty of buffer)
#SBATCH --mem-per-cpu=500MB         # Reduced from 2GB (more efficient)
#SBATCH --output=runs/SLURM_output/%j/pest_ies_%j.out
#SBATCH --error=runs/SLURM_output/%j/pest_ies_%j.err
#SBATCH --nodes=1                   # Add node specification
#SBATCH --partition=compute         # Adjust to your HPC partition name

# Load required modules
module load python/3.10

# Activate environment
source ~/venvs/pakipaki/bin/activate

# Set PEST++ worker count
export PEST_NWORKERS=250

# Navigate to project
cd /path/to/pakipaki_01

# Run PEST++ IES
python scripts/e_pest_IES_HM.py

echo "PEST++ IES completed at $(date)"
```

**Key changes**:
1. ✅ Reduced `--mem-per-cpu` from 2GB → 500MB (saves resources, faster scheduling)
2. ✅ Reduced `--time` from 2:00:00 → 1:00:00 (still safe, helps scheduling)
3. ✅ Added `--nodes=1` (explicit single-node)
4. ✅ Added environment activation
5. ✅ Added completion message

**Expected performance with 250 workers**:
- Runtime: ~18 minutes
- Total memory: 250 × 500 MB = 125 GB
- Speedup: 31.4× vs local

---

## 📋 Memory Sizing Guide

### Actual Memory Usage Per Worker

From your model:
- **MODFLOW memory**: 2.98 MB (from mfsim.lst)
- **Python overhead**: ~50-100 MB (numpy, pandas, flopy)
- **File I/O buffers**: ~20 MB
- **Total active**: ~100-150 MB per worker

### Recommended Memory Allocations

| Scenario | mem-per-cpu | Total (250 workers) | Notes |
|----------|-------------|---------------------|-------|
| **Tight** | 200 MB | 50 GB | Risky, may fail on large stress periods |
| **Safe** ⭐ | 500 MB | 125 GB | Recommended (3-5× actual usage) |
| **Conservative** | 1 GB | 250 GB | Very safe but wasteful |
| **Your original** | 2 GB | 500 GB | Excessive (20× actual usage) |

**Recommendation**: Use **500 MB per CPU** for balance between safety and efficiency.

---

## ⏱️ Runtime Sizing Guide

### Expected Runtimes by Worker Count

| Workers | Expected Runtime | Time Limit | Speedup vs Local |
|---------|------------------|------------|------------------|
| 38 | 2 hours | 3:00:00 | 4.7× |
| 76 | 1 hour | 2:00:00 | 9.4× |
| 152 | 30 min | 1:00:00 | 18.9× |
| **250** | **~18 min** | **1:00:00** | **31.4×** |
| 456 | 10 min | 0:30:00 | 56.6× |
| 912 | 5 min | 0:20:00 | 113× |

**Buffer recommendations**:
- Add 2-3× the expected runtime for SLURM time limit
- Accounts for: queue wait variability, I/O overhead, iteration transitions

---

## 🔍 Checking Your HPC Capabilities

Before submitting, check your cluster's resources:

```bash
# Check available partitions and their limits
sinfo -o "%P %.10l %.6D %.6c %.8z %.10m"

# Check node specifications
scontrol show node | grep -E "NodeName|CPUTot|RealMemory"

# Check your account limits
sacctmgr show assoc where user=$USER format=user,account,partition,maxcpus,maxnodes

# Test job allocation (without running)
salloc --cpus-per-task=250 --mem-per-cpu=500MB --time=01:00:00 --test-only
```

**Key questions**:
1. What's the max CPUs per node? (Determines if you need multi-node)
2. What's the max memory per node? (Validates your mem-per-cpu × cpus-per-task)
3. What partition should you use? (compute, long, high-mem, etc.)
4. Are there any account limits? (May restrict total cores)

---

## 🚀 Complete SLURM Submission Script

### Recommended Production Script

```bash
#!/bin/bash
#SBATCH --job-name=pest_ies_local_run33
#SBATCH --account=YOUR_ACCOUNT_NAME      # Replace with your account
#SBATCH --partition=compute              # Replace with your partition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=250
#SBATCH --mem-per-cpu=500MB
#SBATCH --time=01:00:00
#SBATCH --output=runs/SLURM_output/%j/pest_ies_%j.out
#SBATCH --error=runs/SLURM_output/%j/pest_ies_%j.err
#SBATCH --mail-type=END,FAIL            # Email when job ends/fails
#SBATCH --mail-user=your.email@university.edu

# Print job info
echo "========================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start Time: $(date)"
echo "========================================="

# Create output directory if it doesn't exist
mkdir -p runs/SLURM_output/$SLURM_JOB_ID

# Load modules
module purge
module load python/3.10
module load gcc/11.2.0

# Show loaded modules
module list

# Activate virtual environment
source ~/venvs/pakipaki/bin/activate

# Verify Python packages
python -c "import flopy, pyemu, numpy, pandas; print('All packages loaded successfully')"

# Set environment variables
export PEST_NWORKERS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=1  # Prevent nested parallelism
export MKL_NUM_THREADS=1  # Intel MKL threading control

# Navigate to project directory
cd $SLURM_SUBMIT_DIR  # Directory where sbatch was called
# Or use absolute path:
# cd /home/username/pakipaki_01

# Run PEST++ IES
echo "Starting PEST++ IES with $PEST_NWORKERS workers..."
python scripts/e_pest_IES_HM.py

# Check exit status
if [ $? -eq 0 ]; then
    echo "========================================="
    echo "PEST++ IES completed successfully"
    echo "End Time: $(date)"
    echo "========================================="
else
    echo "========================================="
    echo "PEST++ IES failed with exit code $?"
    echo "End Time: $(date)"
    echo "========================================="
    exit 1
fi

# Optional: Copy results to archive location
# cp -r models/local_run33/pest/master_ies /archive/location/

# Print resource usage
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Elapsed,MaxRSS,MaxVMSize,CPUTime,State
```

### Usage

```bash
# Submit the job
sbatch slurm_pest_ies.sh

# Check job status
squeue -u $USER

# Monitor output in real-time
tail -f runs/SLURM_output/JOBID/pest_ies_JOBID.out

# Cancel if needed
scancel JOBID
```

---

## 📊 Quick Reference Table

| Use Case | cpus-per-task | mem-per-cpu | time | Expected Runtime |
|----------|---------------|-------------|------|------------------|
| **Quick test** | 38 | 500MB | 03:00:00 | ~2 hours |
| **Moderate** | 76 | 500MB | 02:00:00 | ~1 hour |
| **Balanced** ⭐ | 152 | 500MB | 01:00:00 | ~30 min |
| **Your config** | 250 | 500MB | 01:00:00 | ~18 min |
| **Fast** | 456 | 500MB | 00:30:00 | ~10 min |
| **Ultra** | 912 | 500MB | 00:20:00 | ~5 min |

---

## 🐛 Troubleshooting

### Job Fails with "Out of Memory"

**Solution**: Increase `--mem-per-cpu` to 1GB:
```bash
#SBATCH --mem-per-cpu=1GB
```

### Job Stays in Queue Too Long

**Solution**: Reduce `--cpus-per-task` to request fewer cores:
```bash
#SBATCH --cpus-per-task=76  # More likely to be scheduled quickly
```

### Workers Not Starting

**Check**: PEST++ may need explicit worker configuration in `e_pest_IES_HM.py`:
```python
# In e_pest_IES_HM.py, look for:
pyemu.os_utils.start_workers(
    worker_root,
    exe_rel_path='pestpp-ies',
    pst_rel_path=pst_name,
    num_workers=os.environ.get('PEST_NWORKERS', os.cpu_count()),  # Use env var
    master_dir=master_dir,
)
```

### Multi-Node Jobs Fail

**Solution**: Use panther master/agent mode or MPI-enabled PEST++:
```bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=114

# In script, use:
mpirun -np 456 pestpp-ies local_run33.pst
```

---

## ✅ Final Recommendation

For your university HPC with 250 CPUs requested:

```bash
#SBATCH --job-name=pest_ies_run33
#SBATCH --cpus-per-task=250
#SBATCH --time=01:00:00              # Reduced from 02:00:00
#SBATCH --mem-per-cpu=500MB          # Reduced from 2GB
#SBATCH --output=runs/SLURM_output/%j/pest_ies_%j.out
#SBATCH --error=runs/SLURM_output/%j/pest_ies_%j.err
```

**Expected Results**:
- Runtime: ~18 minutes (vs 9.4 hours local)
- Memory: 125 GB total (vs 500 GB with 2GB/cpu)
- Speedup: 31.4×
- **Total wall time savings**: 9.1 hours per run!

This configuration is **well-balanced** for HPC scheduling and will likely be allocated quickly.
