# Adaptive Mesh — Execution Log

## Phase 1: Remove local_enthalpy, microstructure, crack_risk

### [2026-03-30 00:53:04] Starting Phase 1 execution

### [2026-03-30 00:53:10] Tasks 1.1-1.4: Replaced delt_eff with delt
- `mod_discret.f90`: removed `use local_enthalpy`, changed `delt_eff(i,j,k)` → `delt`
- `mod_sour.f90`: removed `use local_enthalpy`, changed `delt_eff(i,j,k)` → `delt`
- `mod_flux.f90`: removed `use local_enthalpy`, changed `delt_eff(i,j,k)` → `delt`

### [2026-03-30 00:54:00] Task 1.5: Rewrote main.f90
- Removed `use local_enthalpy`, `use microstructure_mod`, `use crack_risk_mod`
- Removed `is_local` variable and all local/global branching
- Set `ilo/ihi/jlo/jhi/klo/khi` to full domain once before time loop
- Removed `call allocate_skipped`, `call get_enthalpy_region`, `call update_localfield`
- Removed `call compute_delt_eff`, `call update_skipped` and timing blocks
- Removed micro/crack allocate/update/report calls
- Removed local/global wall-clock accumulation
- Simplified convergence (always global criteria)
- Always call `heat_fluxes` and `write_meltpool_history`
- Kept only global update path for end-of-timestep

### [2026-03-30 00:55:00] Task 1.6: Cleaned mod_param.f90
- Removed `local_half_x`, `local_half_y`, `local_depth_z`, `localnum`
- Removed `micro_flag`, `crack_flag`
- Removed `namelist / local_solver /`
- Removed from `namelist / output_control /`: `micro_flag`, `crack_flag`
- Removed `READ(10, NML=local_solver)`

### [2026-03-30 00:55:10] Task 1.7: Cleaned input_param.txt
- Removed `&local_solver` line
- Removed `micro_flag=0, crack_flag=0` from `&output_control`

### [2026-03-30 00:55:20] Task 1.9: Cleaned mod_timing.f90
- Removed `t_skipped_mgmt`, `t_local_step`, `t_global_step`, `n_local_step`, `n_global_step`
- Updated `nmod` from 17 to 16
- Removed local/global step reporting section

### [2026-03-30 00:55:40] Task 1.10: Updated thermal history monitoring points
- Updated 10 monitoring points in `mod_print.f90` for `my_toolpath.crs` (diamond/island scan centered at 1.5mm, 1.5mm)
- Updated Python plot labels

### [2026-03-30 00:56:00] Task 1.11: Deleted files and cleaned compile.sh
- Deleted: `mod_local_enthalpy.f90`, `mod_microstructure.f90`, `mod_crack_risk.f90`
- Removed from `compile.sh`

### [2026-03-30 00:56:30] Task 1.12 Step A: First compile attempt
- Build succeeded
- Grep found `is_local` in `mod_bound.f90` — `bound_enthalpy` still had `is_local` parameter

### [2026-03-30 00:57:00] Fix: Updated mod_bound.f90
- Removed `is_local` parameter from `bound_enthalpy` signature
- Removed entire `if (is_local)` branch (local boundary logic), kept only global path
- Updated call in `main.f90` to remove `.false.` argument

### [2026-03-30 00:57:30] Task 1.12 Step A: Second compile — CLEAN BUILD
- `bash compile.sh` succeeded
- Grep: only 1 comment in `mod_species.f90` mentioning `delt_eff` — acceptable

### [2026-03-30 00:58:11] Task 1.12 Step B: Started reference case
- `input_param.txt`: `case_name='ref_no_local'`, `toolpath_file='./ToolFiles/my_toolpath.crs'`
- Command: `bash run.sh ref_no_local 4 &` (PID 32938)
- Simulation running in background...

---

## Phase 2: Implement Adaptive Mesh Module

### [2026-03-30 00:59:00] Task 2.1: Added adaptive mesh parameters to mod_param.f90
- Added `adaptive_flag`, `amr_local_half_x`, `amr_local_half_y`, `amr_dx_fine`, `remesh_interval`
- Added `namelist / adaptive_mesh /`
- Added `READ(10, NML=adaptive_mesh)` in `read_data`

### [2026-03-30 00:59:10] Task 2.2: Added `&adaptive_mesh` to input_param.txt
- `adaptive_flag=0` (default off)

### [2026-03-30 01:00:00] Task 2.3: Created mod_adaptive_mesh.f90
- Module `adaptive_mesh_mod` with subroutines:
  - `amr_init()` — validates nzx==1/nzy==1, allocates mapping arrays, calls initial regrid
  - `amr_check_remesh(step_idx)` — tracks beam, expands for melt pool, 20% check, clamps to domain
  - `amr_regenerate_grid()` — saves old grid, generates new 1D grids, recomputes geometry, builds maps, interpolates fields
  - `amr_generate_1d()` — core 1D grid generation with left coarse (geometric), refined (uniform), right coarse (geometric)
  - `amr_recompute_geometry()` — recalculates dxpwinv, dypsinv, fracx, fracy, volumes, areas (with OpenMP)
  - `amr_build_map_1d()` — monotone sweep to find bracketing indices + weights
  - `amr_interpolate_all_fields()` — interpolates 18 solution/property fields
  - `amr_interp_field()` — bilinear X-Y interpolation with precomputed maps (OpenMP parallel)
  - `amr_validate_grid()` — checks domain length, positive volumes, positive spacing inverses

### [2026-03-30 01:01:00] Task 2.4: Added mod_adaptive_mesh.f90 to compile.sh
- After `mod_dimen.f90`

### [2026-03-30 01:02:00] Task 2.5: Integrated into main.f90
- Added `use adaptive_mesh_mod`
- `amr_init()` after allocations (before `init_thermal_history`)
- AMR check/regenerate/validate block in time loop after laser_beam
- `update_thermal_history_indices()` and `defect_update_map()` after remesh

### [2026-03-30 01:03:00] Compile test — CLEAN BUILD

### [2026-03-30 01:03:30] Task 2.7: Added update_thermal_history_indices() to mod_print.f90
- Promoted `px`, `py` to module-level `thist_px`, `thist_py`
- New subroutine re-finds nearest cell indices in X-Y after remesh

### [2026-03-30 01:04:00] Task 2.8: Updated mod_defect.f90 for AMR support
- When `adaptive_flag==1`: allocates uniform defect mesh (`def_ni x def_nj x nk`) at `amr_dx_fine` resolution
- `defect_update_map()` builds index mapping from AMR mesh to uniform mesh (called after each remesh)
- `defect_interp_temp()` interpolates temp from AMR to uniform mesh (bilinear, OpenMP)
- `update_max_temp()` uses interpolated `temp_def` on uniform mesh when AMR is on
- `compute_defect_determ()` and `write_defect_vtk()` use `def_x/def_y` coordinates when AMR is on
- When `adaptive_flag==0`: all behavior unchanged

### [2026-03-30 01:05:00] Final compile — CLEAN BUILD
- All Phase 2 code compiles successfully

---

## Phase 3: Update Documentation

### [2026-03-30 01:05:30] Tasks 3.1-3.3: Updated docs
- `docs/architecture/modules.md`: replaced local_enthalpy with adaptive_mesh, removed micro/crack sections
- `docs/architecture/solver.md`: updated call flow for AMR, removed local/global scheduling
- `docs/architecture/overview.md`: updated program flow and dependency graph
- `docs/input-reference.md`: removed `&local_solver`, `micro_flag`, `crack_flag`; added `&adaptive_mesh`
- `docs/output.md`: removed micro/crack output files
- `docs/modules/physics.md`: replaced local_enthalpy with adaptive_mesh
- `docs/modules/postprocess.md`: removed micro/crack sections
- `docs/index.md`: updated feature list
- `docs/species/overview.md`: removed `delt_eff` reference
- `docs/modules/solver.md`: updated bound_enthalpy and discretization notes

---

## Phase 4: Validation

### [2026-03-30 01:13:00] Reduced timax to 0.005 for faster validation
- Killed old ref run (was 90ms, would take hours)
- Set timax=0.005 (5ms, ~2.5 tracks, 250 timesteps)

### [2026-03-30 01:14:48] Started reference case (timax=0.005)

### [2026-03-30 01:30:54] Reference case completed
- Wall time: 612.1s, CPU time: 3236.1s, Memory: 1335.2 MB

### [2026-03-30 01:31:30] Started AMR case — first attempt FAILED
- NaN in enthalpy residual — `amr_init()` called `amr_regenerate_grid()` before `beam_pos` was set
- Fix: removed initial regrid from `amr_init()`

### [2026-03-30 01:34:00] AMR case — second attempt FAILED
- `Non-positive volume at i,j,k=240,2,2` after first remesh
- Root cause: expansion ratio 1.3 with 162 cells → `1.3^162 ≈ 10^18.5` → `dx0 ≈ 0`
- Fix: replaced fixed ratio with `amr_find_ratio()` bisection solver

### [2026-03-30 01:42:00] AMR case — third attempt: RUNNING
- No errors, no NaN, remesh working correctly

### [2026-03-30 01:55:52] AMR case completed
- Wall time: 601.9s, CPU time: 2204.6s, Memory: 1399.6 MB

### [2026-03-30 01:59:00] Re-ran reference case (clean, no leftover data)

### [2026-03-30 02:12:00] Reference case completed (clean run)

### [2026-03-30 02:15:00] First comparison (400x400, uniform=AMR=10um)
- AMR grid was uniform everywhere — 400 cells in 4mm = 10um = amr_dx_fine
- AMR had no effect on grid spacing
- Negative melt pool lengths detected during cooling

### [2026-03-30 02:20:00] Bug fix: pool_size negative length
- **Root cause**: during cooling, beam center temp < tsolid → length interpolation extrapolated to negative
- **Fix**: rewrote `mod_dimen.f90` pool_size with proper cooling-phase detection and clamping
- Also fixed width detection: replaced `cycle outer_jmax` pattern with sequential j-scan, added value clamping

### [2026-03-30 02:25:00] Re-ran with 200x200 cells (uniform 20um vs AMR 10um refined)
- ref_200: 200x200x52 uniform, 129.5s wall time, 340 MB
- amr_200: 200x200x52 AMR, 197.3s wall time, 429 MB

### [2026-03-30 02:35:00] AMR grid verified in VTK
- Coarse cells: 252um → 186um → 137um → ... (geometric expansion from laser)
- Refined cells: exactly 10um near laser position
- Grid non-uniformity clearly visible

### [2026-03-30 02:40:00] Final comparison (200x200)
- Thermal history: max error 0.35 K — PASS
- Defect: within 0.72% — PASS
- Negative melt pool length: 0 in both cases — FIXED
- Width still has large relative errors at transitions — known limitation of pool_size interpolation
- See `results.md` for full data

### [2026-03-30 02:45:00] Restored input_param.txt to defaults
- 400x400x52, timax=0.09, B26.crs, adaptive_flag=0

### [2026-03-30 03:00:00] Bug fix: pool_size zero-drop during turnaround
- **Symptom**: meltpool l/d/w drops to 0 at t=3.81e-4 while tpeak=1834K >> tsolid=1563K
- **Root cause**: during laser-off turnaround, beam position `(istart,jstart)` moves to next track start, but melt pool is still solidifying at previous track. `pool_size` searches from beam position → finds nothing
- **Fix**: when `temp(istart,jstart,nk) <= tsolid`, find `tpeak` location `(ic,jc)` in domain and use it as search center instead of beam position
- All length/depth/width searches now use `(ic,jc)` instead of `(istart,jstart)`

