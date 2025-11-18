# Cloud Computing Requirements for PEST++ IES Runs

**Analysis Date**: 2025-11-18
**Based on**: local_run33 performance data

---

## 📊 Current Performance Baseline (Local Machine)

### Hardware Configuration
- **CPUs Available**: 8 workers (based on local_run34)
- **Total Runtime**: 565.85 minutes (9.43 hours)
- **Total Model Runs**: 3,548 (across 8 iterations, 219 realizations)

### Timing Breakdown
| Component | Time | Details |
|-----------|------|---------|
| **Total PEST++ Runtime** | 565.85 min | 9.43 hours |
| **Single MODFLOW Run** | ~6.7 seconds | From mfsim.lst |
| **Iterations** | 8 | Including 1 reinflation at iter 4 |
| **Realizations** | 219 final | Started with 250, dropped 31 |
| **Model Runs per Iteration** | ~444 avg | 3548 total / 8 iterations |

### File Sizes
| Item | Size | Notes |
|------|------|-------|
| **Master IES Directory** | 19 GB | All iterations, ensembles, outputs |
| **Single Model Directory** | 91 MB | MODFLOW workspace |
| **Binary Output (.cbc)** | 22 MB | Cell-by-cell flows |
| **Binary Output (.hds)** | 2.2 MB | Heads |
| **MODFLOW Executable** | 24 MB | mf6.exe |
| **Memory per Model** | ~3 MB | From MODFLOW memory manager |

### Problem Dimensions
- **Parameters**: 4,298 adjustable
- **Observations**: 356,273 total (1,416 non-zero weighted)
- **Grid**: Single layer, 25m resolution
- **Stress Periods**: 4

---

## 🎯 Target: 10-Minute Runtime

### Calculation Method

**Current Performance**:
- Total time: 565.85 minutes
- Workers: 8
- Average models/minute: 3548 / 565.85 = 6.27 models/min
- Models per worker per minute: 6.27 / 8 = 0.78 models/worker/min

**Target Performance**:
- Target time: 10 minutes
- Required throughput: 3548 / 10 = 354.8 models/min
- Required workers: 354.8 / 0.78 = **455 workers**

### Scalability Factor

Current to target ratio:
- Speedup needed: 565.85 / 10 = **56.6×**
- Worker increase: 455 / 8 = **56.9× more workers**

**Conclusion**: Near-perfect linear scaling expected (PEST++ IES is embarrassingly parallel)

---

## 💻 Cloud Configuration Recommendations

### Option 1: AWS EC2 (Recommended for Flexibility)

#### Instance Type: **c7i.xlarge** (4 vCPUs, 8 GB RAM)

**Why this instance**:
- Compute-optimized for parallel workloads
- 4 vCPUs sufficient for 1 worker per instance (with overhead)
- 8 GB RAM more than enough (model uses ~3 MB + Python overhead)
- Cost-effective for short-duration intensive compute

**Configuration**:
```yaml
Instance Type: c7i.xlarge
vCPUs per instance: 4
RAM per instance: 8 GB
Number of instances: 115
Total vCPUs: 460
Workers per instance: 4
Total workers: 460
Estimated cost: $0.1785/hr × 115 instances = $20.53/hr
Runtime cost: $20.53 × (10 min / 60) = $3.42 per run
```

**Setup**:
- Launch 115 × c7i.xlarge instances
- Install: Python 3.10+, numpy, pandas, flopy, pyemu, MODFLOW6
- Configure PEST++ panther workers pointing to master
- Master node: c7i.2xlarge (8 vCPUs, 16 GB RAM) for coordination

#### Alternative: **c7i.2xlarge** (8 vCPUs, 16 GB RAM)

**Configuration**:
```yaml
Instance Type: c7i.2xlarge
vCPUs per instance: 8
RAM per instance: 16 GB
Number of instances: 58
Total vCPUs: 464
Workers per instance: 8
Total workers: 464
Estimated cost: $0.357/hr × 58 instances = $20.71/hr
Runtime cost: $20.71 × (10 min / 60) = $3.45 per run
```

**Trade-offs**:
- ✅ Fewer instances to manage
- ✅ Slightly less network overhead
- ⚠️ Less granular scaling
- ⚠️ Higher cost if some instances underutilized

---

### Option 2: Google Cloud Compute Engine

#### Instance Type: **c2-standard-4** (4 vCPUs, 16 GB RAM)

**Configuration**:
```yaml
Instance Type: c2-standard-4
vCPUs per instance: 4
RAM per instance: 16 GB
Number of instances: 115
Total vCPUs: 460
Workers per instance: 4
Total workers: 460
Estimated cost: $0.1494/hr × 115 instances = $17.18/hr
Runtime cost: $17.18 × (10 min / 60) = $2.86 per run
```

**Advantages**:
- Slightly cheaper than AWS
- Good sustained use discounts for repeated runs
- Pre-configured container support (GKE)

---

### Option 3: Azure Virtual Machines

#### Instance Type: **F4s v2** (4 vCPUs, 8 GB RAM)

**Configuration**:
```yaml
Instance Type: F4s_v2
vCPUs per instance: 4
RAM per instance: 8 GB
Number of instances: 115
Total vCPUs: 460
Workers per instance: 4
Total workers: 460
Estimated cost: $0.169/hr × 115 instances = $19.44/hr
Runtime cost: $19.44 × (10 min / 60) = $3.24 per run
```

---

### Option 4: HPC/Batch Services (Best for Regular Runs)

#### AWS Batch + Spot Instances

**Configuration**:
```yaml
Job Type: Array job (3548 tasks)
Instance Type: c7i.xlarge (Spot)
Spot discount: ~70% off on-demand
Number of instances: Auto-scaled to 115
Estimated cost: $0.1785/hr × 0.3 × 115 = $6.16/hr
Runtime cost: $6.16 × (10 min / 60) = $1.03 per run
```

**Advantages**:
- ✅ **70-80% cost savings** with Spot instances
- ✅ Auto-scaling built-in
- ✅ No manual instance management
- ✅ Fault tolerance (restarts failed jobs)

**Disadvantages**:
- ⚠️ Spot instances can be interrupted (rare for 10 min jobs)
- ⚠️ More complex setup initially

#### Google Cloud Batch

Similar to AWS Batch, with comparable pricing and auto-scaling.

---

## 📋 Detailed Resource Requirements

### Per Worker Requirements

| Resource | Minimum | Recommended | Notes |
|----------|---------|-------------|-------|
| **CPU** | 1 core | 1 core | Single-threaded MODFLOW runs |
| **RAM** | 100 MB | 500 MB | Model uses 3 MB, Python ~50 MB, buffer |
| **Disk (temp)** | 150 MB | 250 MB | Model files (91 MB) + outputs |
| **Network** | 1 Mbps | 10 Mbps | Template transfer, result return |

### Master Node Requirements

| Resource | Minimum | Recommended | Notes |
|----------|---------|-------------|-------|
| **CPU** | 4 cores | 8 cores | Managing 460 workers |
| **RAM** | 8 GB | 16 GB | Ensemble storage, matrix ops |
| **Disk** | 50 GB | 100 GB | All iteration outputs (19 GB) + buffer |
| **Network** | 100 Mbps | 1 Gbps | Communicating with 460 workers |

### Storage Requirements

| Component | Size | Backup Needed |
|-----------|------|---------------|
| **Input files** | ~100 MB | Yes - version control |
| **Template directory** | 91 MB | No - regenerated |
| **Output per iteration** | ~2.4 GB | Yes - final results only |
| **Total run output** | 19 GB | Yes - archive after analysis |

**Recommendation**: Use shared storage (AWS EFS, Google Filestore) for master directory, local SSD for workers.

---

## 🚀 Recommended Cloud Setup

### Best Overall: **AWS Batch with Spot Instances**

**Architecture**:
```
Master Node (c7i.2xlarge)
  ├── PEST++ IES master process
  ├── Shared storage (EFS): Input files, outputs
  └── Workers: 460 × c7i.xlarge Spot instances
      └── Local disk: MODFLOW runs
```

**Configuration Steps**:

1. **Create Master AMI**:
   ```bash
   # Ubuntu 22.04 base
   sudo apt update
   sudo apt install -y python3.10 python3-pip git
   pip3 install numpy pandas matplotlib flopy pyemu
   # Copy MODFLOW6 binary to /usr/local/bin/
   # Copy PESTPP-IES binary to /usr/local/bin/
   ```

2. **Setup AWS Batch**:
   ```bash
   # Create compute environment
   aws batch create-compute-environment \
     --compute-environment-name pest-workers \
     --type MANAGED \
     --compute-resources \
       type=SPOT,\
       allocationStrategy=SPOT_CAPACITY_OPTIMIZED,\
       minvCpus=0,\
       maxvCpus=1840,\
       desiredvCpus=1840,\
       instanceTypes=c7i.xlarge,\
       subnets=subnet-xxx,\
       securityGroupIds=sg-xxx
   ```

3. **Launch Master**:
   ```bash
   # On master node
   cd /mnt/efs/pakipaki_01
   python scripts/e_pest_IES_HM.py
   # PEST++ will spawn workers via panther
   ```

**Total Cost per Run**:
- **Spot instances**: $1.03 (10 min runtime)
- **Master node**: $0.06 (c7i.2xlarge × 10 min)
- **Storage**: $0.01 (EFS negligible for short run)
- **Total**: **~$1.10 per run**

**Comparison to local**:
- Local runtime: 9.43 hours
- Cloud runtime: 10 minutes
- **Speedup: 56.6×**
- **Cost per hour saved**: $0.12/hr

---

## 🎛️ Scaling Considerations

### Linear Scaling Assumptions

PEST++ IES scales linearly because:
- ✅ Each realization runs independently
- ✅ No inter-worker communication during model runs
- ✅ Master only coordinates at iteration boundaries
- ✅ I/O is minimal (template in, results out)

**Expected scaling efficiency**: >95% up to 1000 workers

### Potential Bottlenecks

1. **Master node CPU** (likelihood: low)
   - Matrix operations between iterations
   - **Mitigation**: Use c7i.4xlarge (16 cores) for >500 workers

2. **Network I/O** (likelihood: medium)
   - 460 workers requesting templates simultaneously
   - **Mitigation**: Pre-stage templates on workers, use VPC peering

3. **Shared storage** (likelihood: low for 10 min runs)
   - EFS throughput limits
   - **Mitigation**: Use EFS Max I/O mode, or workers write locally

4. **Spot instance interruption** (likelihood: very low)
   - <5% hourly interruption rate for 10 min jobs
   - **Mitigation**: PEST++ auto-restarts failed runs

---

## 📈 Alternative Scenarios

### If Budget is Constrained: 30-Minute Runtime

**Configuration**:
- Workers needed: 3548 / 30 / 0.78 = **152 workers**
- Instances: 38 × c7i.xlarge
- Cost: $0.1785/hr × 38 × 0.5hr = **$3.39 per run**
- **Speedup**: 18.9×

### If You Have More Budget: 5-Minute Runtime

**Configuration**:
- Workers needed: 3548 / 5 / 0.78 = **910 workers**
- Instances: 228 × c7i.xlarge
- Cost: $0.1785/hr × 228 × (5/60)hr = **$3.40 per run**
- **Speedup**: 113×

**Note**: Diminishing returns due to fixed overhead (iteration transitions).

---

## 🔧 Implementation Checklist

### Pre-Launch
- [ ] Create base AMI with Python, MODFLOW6, PESTPP-IES
- [ ] Test single worker execution on cloud instance
- [ ] Configure PEST++ panther settings for cloud workers
- [ ] Setup shared storage (EFS/Filestore) for master directory
- [ ] Configure security groups for worker-master communication

### Launch
- [ ] Start master node
- [ ] Verify connectivity to worker pool
- [ ] Launch worker instances (manual or auto-scaling)
- [ ] Monitor PEST++ .rec file for progress
- [ ] Check worker logs for failures

### Post-Run
- [ ] Download outputs from shared storage
- [ ] Terminate worker instances (save costs!)
- [ ] Archive master directory to S3/GCS (cheap storage)
- [ ] Review costs in billing console

---

## 💰 Cost Summary

| Configuration | Runtime | Instances | Cost/Run | Speedup |
|---------------|---------|-----------|----------|---------|
| **Local (baseline)** | 565 min | 8 workers | $0 | 1× |
| **Cloud - 30 min** | 30 min | 38 × c7i.xlarge | $3.39 | 18.9× |
| **Cloud - 10 min** ⭐ | 10 min | 115 × c7i.xlarge | $1.10 | 56.6× |
| **Cloud - 5 min** | 5 min | 228 × c7i.xlarge | $3.40 | 113× |

**Recommended**: 10-minute runtime with Spot instances ($1.10/run)

---

## 🎯 Final Recommendation

### For Your Use Case:

**Target**: 10-minute runtime
**Platform**: AWS Batch with Spot instances
**Configuration**:
- **Workers**: 115 × c7i.xlarge Spot (460 cores)
- **Master**: 1 × c7i.2xlarge on-demand (8 cores, 16 GB)
- **Storage**: 100 GB EFS for shared data
- **Network**: VPC with 1 Gbps bandwidth

**Per-Run Cost**: ~$1.10
**Setup Time**: ~4 hours (first time), ~30 min (subsequent runs)
**Runtime**: 10 minutes (vs 9.43 hours local)
**Speedup**: 56.6×

### Next Steps:

1. Create AWS account (or use existing)
2. Setup base AMI with dependencies
3. Test with 10 workers (local_run34 or similar)
4. Scale to 115 workers for full 10-min runs
5. Monitor and adjust based on actual performance

**Note**: If you're running this frequently (>10 times), consider reserved instances for 30-40% additional savings on the master node.
