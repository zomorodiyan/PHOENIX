# Build & Run

## Requirements

- **Compiler**: gfortran (GCC Fortran) with OpenMP support
- **Python 3**: For toolpath generation and thermal history plotting (`numpy`, `matplotlib`)
- **ParaView**: For VTK visualization (optional)

## Install dependencies

A helper script is provided that installs everything above on common Linux
distributions (apt / dnf / yum / pacman / zypper). It auto-detects the package
manager, **skips anything already installed**, and verifies the toolchain
(including a small OpenMP compile test):

```bash
cd fortran_new
bash install_deps.sh
```

The script installs:

- `gfortran` and `libgomp` (OpenMP runtime)
- `python3` and `pip`
- `numpy` and `matplotlib` (skipped if already importable, e.g. via pip)

## Build

```bash
cd fortran_new
bash compile.sh
```

This cleans previous build artifacts (`.o`, `.mod`), compiles all modules in dependency order with `-fopenmp -O3 -march=native`, and links into `cluster_main`.

## Run

```bash
bash run.sh <case_name> [omp_threads] &
```

| Argument | Description | Default |
|----------|-------------|---------|
| `case_name` | Name for this run (creates `result/<case_name>/` directory) | Required |
| `omp_threads` | Number of OpenMP threads | 4 |

The script automatically updates `case_name` in `input_param.txt` and sets `OMP_NUM_THREADS`.

**Examples:**

```bash
bash run.sh baseline 4 &          # 4 threads
bash run.sh highpower 8 &         # 8 threads
```

## Monitor

```bash
# Watch output in real-time
tail -f result/mycase/mycase_output.txt

# Check progress (look for "progress%" line)
grep "progress" result/mycase/mycase_output.txt | tail -1
```

## Stop

```bash
# Stop all running simulations
kill $(pgrep -f cluster_main)

# Stop a specific run
ps aux | grep cluster_main    # Find PID
kill <PID>
```

## Clean

```bash
bash clean.sh    # Removes .o, .mod, cluster_main (preserves results)
```

## Workflow

1. `bash install_deps.sh` (first time only — installs gfortran, python3, numpy, matplotlib, OpenMP)
2. Edit `inputfile/input_param.txt` (geometry, materials, numerics)
3. Generate or select a toolpath in `ToolFiles/`
4. `bash compile.sh`
5. `bash run.sh mycase 4 &`
6. Open VTK files in ParaView: `result/mycase/mycase_vtkmov*.vtk`

## ML Experiment Automation (Parallel + Isolated Cases)

PHOENIX includes an automation runner for grid sweeps that is race-condition-safe and ML-data-ready.

- Generates per-case toolpaths (`.crs`) from swept `scan_speed` and `hatch_spacing`
- Materializes isolated case workspaces under `fortran_new/automation/runs/<run_id>/`
- Runs cases in parallel batches (`max_concurrent_runs`)
- Collects KPI table + VTK artifacts index under `fortran_new/automation/results/<run_id>/`

### 1) Configure sweep

Edit example config:

```bash
fortran_new/automation/configs/ml_grid_example.json
```

Key fields:

- `sweep`: Cartesian grid of parameters (include `scan_speed`, `hatch_spacing`)
- `toolpath.template`: fixed geometry defaults for toolpath generation
- `execution.max_concurrent_runs`: number of parallel PHOENIX jobs
- `ml_outputs.vtk_patterns`: mesh files to collect (default: defect, maxtemp, vtkmov)

### 2) Run dry-run (materialize only)

```bash
cd fortran_new
python3 automation/run_experiments.py \
  --config automation/configs/ml_grid_example.json \
  --run-id test_dryrun \
  --dry-run
```

### 3) Run experiments

```bash
cd fortran_new
python3 automation/run_experiments.py \
  --config automation/configs/ml_grid_example.json \
  --run-id train_001
```

### 4) Resume interrupted run

```bash
cd fortran_new
python3 automation/run_experiments.py \
  --config automation/configs/ml_grid_example.json \
  --run-id train_001 \
  --resume
```

### Output files

Generated in `fortran_new/automation/results/<run_id>/`:

- `cases_manifest.csv` - per-case status, pid, retries, timing, params
- `kpi_summary.csv` - scalar KPI table per successful case
- `failures.csv` - failed cases with reason/log path
- `ml_index.csv` - collected VTK file paths, bytes, SHA256 checksums
- `run_metadata.json` - config snapshot + execution metadata
