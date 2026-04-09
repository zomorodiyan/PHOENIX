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

1. Edit `inputfile/input_param.txt` (geometry, materials, numerics)
2. Generate or select a toolpath in `ToolFiles/`
3. `bash compile.sh`
4. `bash run.sh mycase 4 &`
5. Open VTK files in ParaView: `result/mycase/mycase_vtkmov*.vtk`
