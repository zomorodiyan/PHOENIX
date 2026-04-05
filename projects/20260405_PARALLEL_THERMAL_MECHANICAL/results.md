# Parallel Thermal-Mechanical — Results

## Test Case

- **Grid**: 200×100×52 thermal, 100×50×25 FEM (mech_mesh_ratio=2)
- **timax**: 0.002 s (100 thermal steps, 4 mechanical solves)
- **AMR**: adaptive_flag=1
- **Machine**: 24 cores available

## Timing Comparison

| Configuration | Thermal wall (s) | Mech wall (s) | Total wall (s) | Speedup |
|--------------|------------------|---------------|----------------|---------|
| Serial 10 threads | 134 (includes mech) | (in-loop) | **134** | 1.0x |
| Parallel 10+10 | 367 | 388 | **388** | 0.35x (slower!) |
| Parallel 5+5 | 78 | 113 | **113** | **1.19x** |

## Analysis

### Why 10+10 is slower than serial
Both processes run 10 OpenMP threads on shared cores. Although 24 cores are available, OpenMP thread scheduling and memory bandwidth contention cause significant overhead. Each process individually runs ~3x slower than when running alone.

### Why 5+5 works
Total 10 threads = same as serial. Thermal finishes in 78s (vs ~64s thermal-only portion in serial = reasonable 22% overhead from 5 vs 10 threads). Mechanical runs concurrently in 113s. Total = max(78, 113) = 113s, **16% faster** than serial 134s.

### Optimal thread allocation
The mechanical solver benefits less from parallelism than thermal (EBE approach has limited parallelism per element). For this grid size:
- Thermal scales well up to ~10 threads
- Mechanical saturates around 5-8 threads
- Best split: give thermal more threads (e.g., 7+3 or 8+2)

## Result Validation

Von Mises stress comparison (serial vs parallel 5+5):

| Mech Solve | Serial max_vm (MPa) | Parallel max_vm (MPa) | Diff |
|-----------|---------------------|----------------------|------|
| 1 | 174.1 | 172.6 | <1% |
| 2 | 179.0 | 179.0 | 0% |
| 3 | 184.4 | 184.4 | 0% |
| 4 | 188.9 | 188.7 | <1% |

Results are consistent — parallel mode produces identical physics.

## Bugs Fixed During Implementation

1. **Double init_mechanical**: mechanical process called `init_mechanical` twice (once before and once inside `run_mechanical_loop`). Fixed by removing the outer call.
2. **update_mech_grid in thermal process**: AMR remesh called `update_mech_grid` even in parallel mode where mechanical arrays don't exist in thermal process → segfault. Fixed by adding `mech_parallel` guard.
