# Adaptive Mesh Refinement (X-Y Movable)

## Objective

Implement a movable adaptive structured mesh in the X-Y plane that follows the laser/melt pool, with fixed Z-direction mesh. The refined region dynamically tracks the melt pool and expands if the melt pool grows beyond the initial refinement zone. This replaces the removed local enthalpy solver with a true mesh refinement approach.

## Execution Tracking

- **`log.md`**: Step-by-step execution log with system timestamps (`[YYYY-MM-DD HH:MM:SS]`). Records commands, compilation output, errors, and decisions.
- **`results.md`**: Final validation data, comparison tables, and performance metrics.

---

## Phase 1: Remove Local Enthalpy Solver

**Goal**: Completely remove `mod_local_enthalpy.f90` and all references. After removal, all solves are global with uniform `delt`.

### Task 1.1: Remove `mod_local_enthalpy.f90` from build

- **File**: `compile.sh`
- **Change**: Remove line `mod_local_enthalpy.f90 \` from the gfortran compilation list (line 24)

### Task 1.2: Replace `delt_eff` with `delt` in `mod_discret.f90`

- **File**: `mod_discret.f90`
- **Line 12**: Remove `use local_enthalpy, only: delt_eff`
- **Line 361**: Change `apnot(i,j,k)=den(i,j,k)/delt_eff(i,j,k)*volume(i,j,k)` to `apnot(i,j,k)=den(i,j,k)/delt*volume(i,j,k)`

### Task 1.3: Replace `delt_eff` with `delt` in `mod_sour.f90`

- **File**: `mod_sour.f90`
- **Line 13**: Remove `use local_enthalpy, only: delt_eff`
- **Line 213**: Change `volht=volume(i,j,k)*hlatnt*den(i,j,k)/delt_eff(i,j,k)` to `volht=volume(i,j,k)*hlatnt*den(i,j,k)/delt`

### Task 1.4: Replace `delt_eff` with `delt` in `mod_flux.f90`

- **File**: `mod_flux.f90`
- **Line 10**: Remove `use local_enthalpy, only: delt_eff`
- **Line 77**: Change `accul=accul+volume(i,j,k)*den(i,j,k)*dh1/delt_eff(i,j,k)` to `accul=accul+volume(i,j,k)*den(i,j,k)*dh1/delt`

### Task 1.5: Remove all local_enthalpy usage from `main.f90`

- **File**: `main.f90`
- **Line 33**: Remove `use local_enthalpy`
- **Line 35**: Remove `use microstructure_mod`
- **Line 36**: Remove `use crack_risk_mod`
- **Line 44**: Remove `integer ilo, ihi, jlo, jhi, klo, khi`
- **Line 45**: Remove `logical is_local`
- **Line 58**: Remove `call allocate_skipped(ni, nj, nk)`
- **Line 60**: Remove `if (micro_flag == 1) call allocate_microstructure(ni, nj, nk)`
- **Line 61**: Remove `if (crack_flag == 1) call allocate_crack_risk(ni, nj, nk)`
- **Line 91**: Remove `call get_enthalpy_region(...)` — replace with setting full-domain indices:
  ```fortran
  ilo = 2; ihi = nim1; jlo = 2; jhi = njm1; klo = 2; khi = nkm1
  ```
  (Note: `ilo/ihi/jlo/jhi/klo/khi` are still needed as arguments to solver subroutines. Keep them as local variables but always set to full domain.)
- **Line 92**: Remove `call update_localfield(ilo, ihi, jlo, jhi, klo, khi)`
- **Lines 96-99**: Remove `call compute_delt_eff()` and its timing block
- **Lines 101-105**: Remove the local/global logging block (or simplify to just log the step)
- **Line 274**: Remove the `if (.not. is_local)` conditional around `heat_fluxes` — always call it
- **Lines 279-284**: Remove the `else` branch that zeroes flux diagnostics
- **Line 290-294**: Simplify convergence: remove `is_local` branches, always use global convergence criteria:
  ```fortran
  if(resorh.lt.conv_res_heat .and. ratio.le.ratio_upper .and. ratio.ge.ratio_lower) exit iter_loop
  ```
- **Line 312**: Remove `call update_skipped(...)` and its timing block
- **Lines 319-320**: Remove `if (micro_flag == 1) call update_microstructure(delt)` and `if (crack_flag == 1) call update_crack_risk(delt)`
- **Lines 340-346**: Remove local/global wall-clock accumulation block
- **Lines 348-414**: Simplify end-of-timestep update — remove the `if (is_local) ... else ...` branching. Keep only the global update path (lines 393-413)
- **Line 418**: Remove `if (.not. is_local)` guard on `write_meltpool_history` — always call it
- **Line 427**: Remove `if (micro_flag == 1) call report_microstructure()`
- **Line 428**: Remove `if (crack_flag == 1) call compute_crack_report()`

### Task 1.6: Clean up `mod_param.f90`

- **File**: `mod_param.f90`
- **Line 28**: Remove `real(wp) local_half_x, local_half_y, local_depth_z`
- **Line 31**: Remove `localnum` from the integer declaration
- **Line 46**: Remove `namelist / local_solver / localnum, local_half_x, local_half_y, local_depth_z`
- **Line 48-49**: Remove `integer micro_flag` and `integer crack_flag`
- **Line 50**: Remove `micro_flag, crack_flag` from `namelist / output_control /`
- **Line 59-60**: Remove default assignments `micro_flag = 0` and `crack_flag = 0`
- **Line 96-97**: Remove `READ (10, NML=local_solver)` and its comment

### Task 1.7: Clean up `input_param.txt`

- **File**: `inputfile/input_param.txt`
- **Line 28**: Remove `&local_solver localnum=4, local_half_x=1.0e-3, local_half_y=2.0e-4, local_depth_z=2.0e-4 /`
- **Line 30**: Remove `micro_flag=0, crack_flag=0` from `&output_control` namelist

### Task 1.8: Remove `localfield` from VTK output (optional — repurpose later)

- **File**: `mod_print.f90` line 244 — the `localfield` VTK output can remain since the field still exists in `mod_field_data.f90`. It will just be all zeros. Will be repurposed in Phase 2 to visualize the refined mesh region.

### Task 1.9: Remove timing variables for skipped management

- **File**: `mod_timing.f90` — remove `t_skipped_mgmt` variable and its reporting
- **File**: `main.f90` — remove timing accumulation for `t_skipped_mgmt`
- **File**: `mod_timing.f90` — also remove `t_local_step`, `t_global_step`, `n_local_step`, `n_global_step` and their reporting if they exist

### Task 1.10: Update thermal history monitoring points for `my_toolpath.crs`

The toolpath is a **circular/diamond island scan** centered at ~(1.5mm, 1.5mm), 26 tracks, hatch spacing 0.1mm, longest tracks at y≈1.15-1.75mm spanning ~2.3mm in x. The current 10 monitoring points (hardcoded in `mod_print.f90`) were designed for `B26.crs` (track at y=0.5mm) and need updating.

- **File**: `mod_print.f90`, subroutine `init_thermal_history`
- Update `px`, `py`, `pz` arrays and the Python plot labels:

| Point | Coordinates (mm) | Description |
|-------|-------------------|-------------|
| P1 | (0.50, 1.45, 0.695) | Central track (track 13), near start, surface |
| P2 | (1.50, 1.45, 0.695) | Central track, mid-length, surface |
| P3 | (2.50, 1.45, 0.695) | Central track, near end, surface |
| P4 | (1.50, 1.45, 0.660) | Central track centre, 40 um depth |
| P5 | (1.50, 1.45, 0.600) | Central track centre, 100 um depth |
| P6 | (1.50, 1.50, 0.695) | Midpoint between track 13 and 14, surface |
| P7 | (1.50, 1.65, 0.695) | 2 hatches from centre (track 15), surface |
| P8 | (1.50, 2.06, 0.695) | Near edge of scanned region (track 19), surface |
| P9 | (1.50, 3.50, 0.695) | Far from scanned region, substrate surface |
| P10 | (1.50, 1.45, 0.200) | Deep substrate below centre |

Also update `finalize_thermal_history` Python labels array to match.

### Task 1.10b: Fix `pool_size` in `mod_dimen.f90`

**Problem**: During laser-off turnaround, `pool_size` searches for the melt pool centered on `(istart, jstart)` — the **current beam position**. But during turnaround the beam has already moved to the next track's start position while the melt pool is still solidifying at the **previous track**. Result: length/depth/width drop to 0 even though `tpeak >> tsolid`.

**Root cause**: All search loops (`do i=istart,...`, width scan from `j=jstart`) use the beam position as origin, not the melt pool position.

**Fix**: When `temp(istart, jstart, nk) <= tsolid` (beam center is below solidus), find the actual `tpeak` location `(ipeak, jpeak)` and use it as the search center.

- **File**: `mod_dimen.f90`
- Add local variables `ic, jc` (search center indices)
- Default `ic=istart, jc=jstart`
- When `temp(istart,jstart,nk) <= tsolid`: scan domain for `tpeak` location → set `ic, jc`
- Replace all `istart`/`jstart` references in length/depth/width search with `ic`/`jc`
- Add `alen = max(0.0_wp, alen)` and `width = max(0.0_wp, width)` clamping
- Add division-by-zero guards: only interpolate when `|T(i)-T(i±1)| > 1 K`
- Width: replaced `cycle outer_jmax` pattern with sequential j-scan to avoid skipping entire i-columns
- Width interpolation: clamp result to `[y(jl), y(jh)]` range

### Task 1.11: Delete removed module files and clean `compile.sh`

- Delete `fortran_new/mod_local_enthalpy.f90`
- Delete `fortran_new/mod_microstructure.f90`
- Delete `fortran_new/mod_crack_risk.f90`
- **File**: `compile.sh` — remove `mod_local_enthalpy.f90 \`, `mod_microstructure.f90 \`, `mod_crack_risk.f90 \` from compilation list

### Task 1.12: Verify — compile, check references, and run reference case

**Step A**: Compile and confirm no dangling references
- Run `bash compile.sh` and confirm clean build
- Run `grep -r "local_enthalpy\|delt_eff\|localnum\|is_local\|n_skipped\|update_skipped\|allocate_skipped\|compute_delt_eff\|get_enthalpy_region\|update_localfield\|local_half_x\|local_half_y\|local_depth_z\|local_solver\|microstructure_mod\|crack_risk_mod\|micro_flag\|crack_flag\|allocate_microstructure\|allocate_crack_risk\|update_microstructure\|update_crack_risk\|report_microstructure\|compute_crack_report" fortran_new/ --include="*.f90"` — should return no matches (except possibly comments)

**Step B**: Run reference simulation with `my_toolpath.crs`
- Modify `input_param.txt`: set `toolpath_file='./ToolFiles/my_toolpath.crs'` and `case_name='ref_no_local'`
- Toolpath covers X: [0.14, 2.86] mm, Y: [0.24, 2.77] mm, Z: 0.6975 mm — fits within the 4mm x 4mm x 0.7mm domain
- Run: `bash run.sh ref_no_local 4 &`
- Record from output:
  - Total wall-clock time
  - Peak memory usage
  - Melt pool dimensions (length, depth, width) at representative timesteps
  - Temperature field at final timestep
- These results serve as the **baseline reference** for Phase 2 AMR comparison:
  - **Accuracy**: compare melt pool dimensions and temperature field between uniform mesh (this run) and adaptive mesh (Phase 2)
  - **Performance**: compare wall-clock time and memory usage
  - **Validation**: AMR results should match uniform mesh within acceptable tolerance

---

## Phase 2: Implement Adaptive Mesh Module

**Goal**: Create `mod_adaptive_mesh.f90` that manages a movable refined region in X-Y, regenerates the grid on demand, and interpolates solution fields onto the new mesh.

### Concept

The domain has a **fixed total number of cells** in X and Y (user-specified). The Z mesh is unchanged. At every `remesh_interval` timesteps, the mesh is regenerated:

1. A **refined region** (rectangular box) is centered on the laser beam position in X-Y
2. Inside the refined region: uniform cells of size `refined_dx` (e.g., 10 um)
3. Outside the refined region: cells are distributed using a **geometric expansion** with ratio `expand_ratio` (default 1.3)
4. The total cell count in X and Y is conserved

### Grid Distribution Algorithm (1D, applied to X and Y independently)

Given:
- Total domain length `L`, total cells `N_total`
- Refined region: `[x_lo, x_hi]` with cell size `dx_fine`
- Number of refined cells: `N_fine = (x_hi - x_lo) / dx_fine`
- Remaining cells: `N_coarse = N_total - N_fine`
- Required: at least 20% of `N_fine` cells on each side outside the refined region

**Left coarse region** (from domain start to `x_lo`):
- Length: `L_left = x_lo`
- Cells distributed with geometric expansion (ratio 1.3) growing away from the refined boundary
- Cell count: proportional to `L_left / (L_left + L_right)` of `N_coarse`

**Right coarse region** (from `x_hi` to domain end):
- Length: `L_right = L - x_hi`
- Same geometric expansion, mirrored

**Minimum cell count check**: Each coarse side must have at least `0.2 * N_fine / 2 = 0.1 * N_fine` cells. If not enough cells are available, the refined region is **not expanded** further and a warning is written to the output file.

### Task 2.1: Add adaptive mesh parameters to `mod_param.f90`

- **File**: `mod_param.f90`
- Add new module-level variables:
  ```fortran
  integer :: adaptive_flag = 0          ! 0=off, 1=on
  real(wp) :: amr_local_half_x = 1.0e-3  ! initial half-length of refined region in x
  real(wp) :: amr_local_half_y = 2.0e-4  ! initial half-length of refined region in y
  real(wp) :: amr_dx_fine = 10.0e-6       ! refined cell size (10 um)
  integer  :: remesh_interval = 20        ! regrid every N timesteps
  ```
- Add new namelist:
  ```fortran
  namelist / adaptive_mesh / adaptive_flag, amr_local_half_x, amr_local_half_y, amr_dx_fine, remesh_interval
  ```
- Read the namelist in `read_data` (after `&output_control`)

### Task 2.2: Add `&adaptive_mesh` to `input_param.txt`

- **File**: `inputfile/input_param.txt`
- Add after `&output_control`:
  ```
  &adaptive_mesh adaptive_flag=0, amr_local_half_x=1.0e-3, amr_local_half_y=2.0e-4, amr_dx_fine=10.0e-6, remesh_interval=20 /
  ```
- Default `adaptive_flag=0` so existing behavior is unchanged
- Note: when `adaptive_flag=1`, the x-zone and y-zone geometry parameters still define total cell count and domain length, but `nzx` and `nzy` must be 1 (single zone). The exponents `powrx`/`powry` are ignored in adaptive mode. Z mesh remains as defined.

### Task 2.3: Create `mod_adaptive_mesh.f90`

- **File**: `mod_adaptive_mesh.f90` (new file)
- **Module**: `adaptive_mesh_mod`
- **Dependencies**: `precision`, `geometry`, `parameters`, `sim_state`, `field_data`, `dimensions`

**Module variables**:
```fortran
real(wp) :: amr_half_x, amr_half_y        ! current refined region half-sizes
real(wp) :: amr_expand_ratio = 1.3_wp      ! geometric expansion ratio for coarse cells
integer  :: amr_step_counter = 0           ! counts steps since last remesh
logical  :: amr_needs_remesh = .false.     ! flag to trigger remesh
```

**Subroutines**:

#### `subroutine amr_init()`
- Save initial half sizes: `amr_half_x = amr_local_half_x`, `amr_half_y = amr_local_half_y`
- Validate: `adaptive_flag==1` requires `nzx==1` and `nzy==1`
- Perform initial mesh generation (call `amr_regenerate_grid`)

#### `subroutine amr_check_remesh(step_idx)`
- Increment `amr_step_counter`
- If `amr_step_counter >= remesh_interval`:
  - Update refined region to track beam: center = `(beam_pos, beam_posy)`
  - Expand if melt pool exceeds initial region: `amr_half_x = max(amr_local_half_x, alen)`, `amr_half_y = max(amr_local_half_y, width)`
  - Clamp to domain boundaries: ensure `[center - half, center + half]` stays within `[0, dimx]` and `[0, dimy]`
  - Check minimum coarse cell count (20% rule). If violated, do NOT expand further and write warning to output file (unit 9):
    ```
    WARNING: Adaptive mesh cannot expand further — insufficient coarse cells (melt pool too large)
    ```
  - Set `amr_needs_remesh = .true.`
  - Reset counter

#### `subroutine amr_regenerate_grid()`
The core algorithm. For X-direction (Y is analogous):

1. Compute refined region bounds:
   ```
   x_lo = max(0, beam_pos - amr_half_x)
   x_hi = min(dimx, beam_pos + amr_half_x)
   ```
2. Number of refined cells: `n_fine_x = nint((x_hi - x_lo) / amr_dx_fine)`
3. Remaining cells: `n_coarse_x = ncvx(1) - n_fine_x`
4. **20% check**: `n_min_coarse = nint(0.2 * n_fine_x)`. If `n_coarse_x < n_min_coarse`, do not expand `amr_half_x` further, write warning, and shrink `n_fine_x` to `ncvx(1) - n_min_coarse`
5. Split coarse cells between left and right proportionally:
   ```
   L_left = x_lo
   L_right = dimx - x_hi
   n_left = nint(n_coarse_x * L_left / (L_left + L_right))
   n_right = n_coarse_x - n_left
   ```
   (Handle edge cases: if `L_left == 0`, all coarse cells go right, etc.)
6. **Generate left coarse cells** (geometric series expanding away from refined boundary):
   - First cell (adjacent to refined region) has size close to `amr_dx_fine`
   - Each subsequent cell (moving left) grows by factor `amr_expand_ratio`
   - Total length must equal `L_left`
   - Solve for first cell size: `dx_0 * (r^n - 1) / (r - 1) = L_left` where `r = amr_expand_ratio`, `n = n_left`
   - If `r == 1`: uniform `dx_0 = L_left / n_left`
7. **Generate refined cells**: uniform `amr_dx_fine` from `x_lo` to `x_hi`
8. **Generate right coarse cells**: mirror of left (expanding away from refined boundary)
9. Assemble the full 1D velocity grid `xu(1:ni)` from these three segments
10. Compute scalar grid `x(i) = (xu(i) + xu(i+1)) / 2`

After generating new X and Y grids:
- Recompute all geometry metrics: `dxpwinv`, `dypsinv`, `fracx`, `fracy`, `volume`, `volume_u`, `volume_v`, `volume_w`, `areaij`, `areajk`, `areaik`, and all staggered variants
- Z grid and Z-related metrics remain unchanged
- **OpenMP**: Parallelize the 3D volume/area recomputation loops (`!$OMP PARALLEL DO`) — same pattern as `generate_grid` in `mod_geom.f90`

#### `subroutine amr_interpolate_fields(x_old, y_old, x_new, y_new, ni_old, nj_old)`
- Interpolate all 3D solution fields from old mesh to new mesh:
  - `enthalpy`, `temp`, `fracl`, `uVel`, `vVel`, `wVel`, `pp`, `pressure`
  - `hnot`, `tnot`, `fraclnot`, `unot`, `vnot`, `wnot`
  - `solidfield`, `max_temp`
  - `concentration` (if `species_flag==1`)
  - `den`, `vis`, `cp`, `thcon`, `gam` (property fields)
- **Interpolation method**: bilinear interpolation in X-Y (Z indices are unchanged, 1-to-1 mapping)
- For each new cell `(i_new, j_new, k)`:
  1. Find bracketing old cell indices: `x_old(i1) <= x_new(i_new) <= x_old(i2)`
  2. Compute weight: `wx = (x_new - x_old(i1)) / (x_old(i2) - x_old(i1))`
  3. Same for Y
  4. `field_new(i,j,k) = (1-wx)*(1-wy)*f(i1,j1,k) + wx*(1-wy)*f(i2,j1,k) + (1-wx)*wy*f(i1,j2,k) + wx*wy*f(i2,j2,k)`
- **Important**: since `ni`, `nj`, `nk` do not change (total cell count is fixed), array dimensions stay the same. No reallocation needed for solution fields.
- Geometry arrays also do not need reallocation (same dimensions), just update values in place.
- **OpenMP**: The interpolation loop over `(i_new, j_new, k)` must be parallelized with `!$OMP PARALLEL DO`. Each new cell's interpolation is independent. Use `PRIVATE` for loop indices and interpolation weights. Pre-compute the index mapping arrays `imap(i_new)` and `jmap(j_new)` (bracketing old indices + weights) in serial 1D sweeps first, then use them in the parallel 3D loop to avoid redundant searches.

### Task 2.4: Add `mod_adaptive_mesh.f90` to `compile.sh`

- **File**: `compile.sh`
- Add `mod_adaptive_mesh.f90 \` after `mod_dimen.f90` (line 23, where `mod_local_enthalpy.f90` used to be)
- It depends on `geometry`, `parameters`, `sim_state`, `field_data`, `dimensions`

### Task 2.5: Integrate adaptive mesh into `main.f90`

- **File**: `main.f90`
- Add `use adaptive_mesh_mod`
- After `generate_grid` and `allocate_fields`, call:
  ```fortran
  if (adaptive_flag == 1) call amr_init()
  ```
- In the time loop, before the iteration loop (after `laser_beam` and `read_coordinates`):
  ```fortran
  if (adaptive_flag == 1) then
      call amr_check_remesh(step_idx)
      if (amr_needs_remesh) then
          call amr_regenerate_grid()
          amr_needs_remesh = .false.
      endif
  endif
  ```

### Task 2.6: Validate geometry consistency

After remesh, the following must hold:
- `sum(xu(i+1)-xu(i) for i=2..nim1)` equals domain length in X (within tolerance)
- `sum(yv(j+1)-yv(j) for j=2..njm1)` equals domain length in Y
- All `volume(i,j,k) > 0`
- All `dxpwinv(i) > 0`, `dypsinv(j) > 0`
- Add a debug check subroutine `amr_validate_grid()` that verifies these and aborts with error message if violated

### Task 2.7: Update thermal history monitoring point indices after remesh

**Problem**: `init_thermal_history` (mod_print.f90) maps 10 physical coordinates to grid indices `thist_i(p), thist_j(p), thist_k(p)` once at startup. After AMR remesh, the X-Y grid positions change, so the same index points to a different physical location, causing temperature history data to jump.

**Solution**: Add `subroutine amr_update_thermal_history_indices()` in `mod_adaptive_mesh.f90` (or in `mod_print.f90` as a public subroutine):
- Store the 10 physical coordinates `px, py, pz` as **module-level arrays** in `mod_print.f90` (currently they are `parameter` locals inside `init_thermal_history` — promote them)
- After each remesh, re-run the nearest-cell search for X and Y only (Z grid is unchanged, so `thist_k` stays the same):
  ```fortran
  subroutine update_thermal_history_indices()
      integer :: p, ip, jp
      real(wp) :: dx, dy, dmin
      do p = 1, n_thist
          dmin = 1.0e30_wp
          do ip = 2, nim1
              dx = (x(ip) - thist_px(p))**2
              if (dx > dmin) cycle
              do jp = 2, njm1
                  dy = dx + (y(jp) - thist_py(p))**2
                  if (dy < dmin) then
                      dmin = dy
                      thist_i(p) = ip
                      thist_j(p) = jp
                  endif
              enddo
          enddo
      enddo
  end subroutine
  ```
- Call this from `amr_regenerate_grid()` at the end, after new grid coordinates are computed
- **Files modified**: `mod_print.f90` (promote `px`, `py` to module-level `thist_px`, `thist_py`; add public subroutine), `mod_adaptive_mesh.f90` (call after remesh)

### Task 2.8: Defect field on independent uniform mesh

**Problem**: Under AMR, the cell index `(i,j,k)` maps to a different physical location after each remesh. `max_temp(i,j,k)` accumulates per-cell history and would drift.

**Solution**: When `adaptive_flag==1`, allocate a **separate uniform mesh** at resolution `amr_dx_fine` covering the full X-Y domain, with the same Z grid. The `max_temp` and `defect_arr` arrays live on this uniform mesh.

#### Uniform defect mesh specification

```fortran
! Uniform mesh dimensions (computed once at init)
integer  :: def_ni, def_nj          ! cell counts: def_ni = nint(dimx / amr_dx_fine) + 2
real(wp), allocatable :: def_x(:)   ! cell-center x positions (1:def_ni)
real(wp), allocatable :: def_y(:)   ! cell-center y positions (1:def_nj)
! Z direction: reuse simulation mesh z(:), nk unchanged
```

- `def_ni = nint(dimx / amr_dx_fine) + 2` (include boundary nodes), similarly `def_nj`
- Example: 4mm domain, 10um cell → `def_ni = 402`, `def_nj = 402` (same as current 400-cell uniform mesh)
- Cell centers: `def_x(i) = (i - 1.5) * amr_dx_fine` for `i = 2..def_ni-1`

#### Per-timestep workflow

1. **Interpolate** `temp` from AMR mesh to uniform defect mesh (bilinear in X-Y, Z is 1:1):
   ```fortran
   subroutine amr_interp_to_defect_mesh(field_amr, field_def)
   ```
   - Pre-compute index mapping `def_imap(i_def)` → bracketing AMR i-indices + weight (changes after each remesh)
   - `!$OMP PARALLEL DO` over the 3D uniform mesh loop

2. **Update** `max_temp` on uniform mesh: `max_temp(i,j,k) = max(max_temp(i,j,k), temp_def(i,j,k))`

#### VTK output

`write_defect_vtk` already writes its own STRUCTURED_GRID VTK with explicit coordinates. Under AMR, use `def_x`, `def_y`, `z` instead of the simulation mesh `x`, `y`, `z`.

#### When `adaptive_flag==0`

No uniform defect mesh is created. Defect module operates directly on the simulation mesh as before.

#### Memory impact

Additional arrays on the uniform defect mesh (only when `adaptive_flag==1`):
- `temp_def`: 1 array × `def_ni × def_nj × nk` (temporary, can be reused)
- For 402×402×52 at 8 bytes: ~67 MB extra

#### Files modified
- `mod_defect.f90`: allocate uniform mesh and defect arrays on it; modify `update_max_temp` to use `temp_def`
- `mod_adaptive_mesh.f90`: add `amr_interp_to_defect_mesh` subroutine and index mapping update after remesh
- `main.f90`: add interpolation call before `update_max_temp` in timestep loop

---

## Phase 3: Update Documentation

### Task 3.1: Update `docs/architecture/modules.md`

- Remove `mod_local_enthalpy.f90`, `mod_microstructure.f90`, `mod_crack_risk.f90` entries
- Add `mod_adaptive_mesh.f90` entry with description

### Task 3.2: Update `docs/input-reference.md`

- Remove `&local_solver` namelist documentation
- Add `&adaptive_mesh` namelist documentation with all parameters

### Task 3.3: Update `docs/architecture/solver.md`

- Update solver call flow to reflect removal of local/global enthalpy scheduling
- Add adaptive mesh remesh step in the time loop flow

---

## Summary of Files Modified

| File | Phase | Action |
|------|-------|--------|
| `mod_local_enthalpy.f90` | 1 | **DELETE** |
| `mod_microstructure.f90` | 1 | **DELETE** |
| `mod_crack_risk.f90` | 1 | **DELETE** |
| `mod_discret.f90` | 1 | Replace `delt_eff` with `delt` |
| `mod_sour.f90` | 1 | Replace `delt_eff` with `delt` |
| `mod_flux.f90` | 1 | Replace `delt_eff` with `delt` |
| `mod_param.f90` | 1+2 | Remove `&local_solver`, `micro_flag`, `crack_flag`; add `&adaptive_mesh` |
| `main.f90` | 1+2 | Remove local_enthalpy/micro/crack logic, add adaptive mesh calls |
| `compile.sh` | 1+2 | Remove deleted modules, add new module |
| `input_param.txt` | 1+2 | Remove `&local_solver`, `micro_flag`, `crack_flag`; add `&adaptive_mesh` |
| `mod_timing.f90` | 1 | Remove skipped/local timing vars |
| `mod_adaptive_mesh.f90` | 2 | **NEW** — core adaptive mesh module |
| `mod_print.f90` | 2 | Promote `px/py` to module-level; add `update_thermal_history_indices()` |
| `mod_defect.f90` | 2 | Uniform defect mesh; `update_max_temp` uses interpolated `temp_def` |
| `docs/architecture/modules.md` | 3 | Update module reference |
| `docs/input-reference.md` | 3 | Update input parameters |
| `docs/architecture/solver.md` | 3 | Update solver call flow |

## Key Design Decisions

1. **Total cell count is conserved** — `ni`, `nj`, `nk` never change. Only the cell positions and sizes change. This means no reallocation of solution arrays.
2. **Z direction is untouched** — the refined region only moves in X-Y. Z mesh remains static as originally defined.
3. **Geometric expansion ratio 1.3** — coarse cells grow by factor 1.3 moving away from the refined boundary. This provides smooth transition.
4. **20% minimum coarse cell rule** — ensures enough resolution outside the refined region. If violated, refined region stops expanding and a warning is issued.
5. **Bilinear interpolation** — solution fields are interpolated from old mesh to new mesh in X-Y only (Z is 1:1).
6. **Remesh interval** — regridding every `remesh_interval` timesteps (default 20) avoids excessive overhead from interpolation.

## Phase 4: Validation — AMR vs Uniform Mesh Comparison

**Goal**: Confirm AMR produces results consistent with the uniform mesh baseline (Phase 1 reference case), and quantify performance gains.

### Task 4.1: Run AMR case with `my_toolpath.crs`

- Set `adaptive_flag=1` in `input_param.txt`, keep all other parameters identical to ref case
- Use same toolpath: `toolpath_file='./ToolFiles/my_toolpath.crs'`
- Set `case_name='amr_test'`
- Run: `bash run.sh amr_test 4 &`

### Task 4.2: Compare physical quantities (output.txt)

For each timestep, compare the following between `ref_no_local` and `amr_test`:
- Enthalpy residual (`resorh`)
- Melt pool dimensions: `alen`, `depth`, `width`
- Peak temperature `tpeak`
- Heat fluxes: `heatout`, `accul`, `ratio`
- Max velocities: `umax`, `vmax`, `wmax`

**Acceptance criteria**: melt pool dimensions within 5% relative error; temperature field within 2% at monitoring points.

### Task 4.3: Compare thermal history

- Diff `ref_no_local_thermal_history.txt` vs `amr_test_thermal_history.txt`
- For each of the 10 monitoring points, compute max absolute temperature difference across all timesteps
- Plot overlay of both cases (can reuse the Python plot script)

### Task 4.4: Compare melt pool history

- Diff `ref_no_local_meltpool_history.txt` vs `amr_test_meltpool_history.txt`
- For each timestep, compare:
  - Melt pool length (`alen`), depth, width
  - Melt pool volume
  - Peak temperature (`Tpeak`)
- Compute max and mean relative error for each quantity over the full simulation
- Plot overlay of both cases: length/depth/width/volume vs time (write a comparison Python script or extend the existing `plot_meltpool_history.py`)
- **Acceptance criteria**: melt pool dimensions within 5% relative error; volume within 10%

### Task 4.5: Compare defect field

- Compare `max_temp` VTK files between the two cases
- Compare `defect_report.txt`: defect fractions (LoF%, keyhole%), total defect volume
- Acceptable tolerance: defect fractions within 1% absolute difference

### Task 4.6: Compare performance

| Metric | `ref_no_local` | `amr_test` | Speedup |
|--------|---------------|------------|---------|
| Total wall-clock time (s) | | | |
| Avg time per step (ms) | | | |
| Peak memory (MB) | | | |

- Wall-clock time from `timing_report.txt`
- Peak memory from `memory_report.txt`

### Task 4.7: Write `results.md`

Create `projects/20260329_ADAPTIVE_MESH/results.md` documenting all execution results:

1. **Phase 1 results**: compilation log, grep verification, reference case output summary
2. **Phase 2 results**: compilation log, any issues encountered during implementation
3. **Phase 4 validation**:
   - Table of melt pool dimension comparison (selected timesteps)
   - Melt pool history comparison (max/mean relative error for length, depth, width, volume)
   - Thermal history max error per monitoring point
   - Defect fraction comparison
   - Performance comparison table
   - Plots (reference paths to generated `.png` files)
4. **Conclusion**: whether AMR passed validation and observed speedup

---

---

## Post-implementation Fixes

### Fix: pool_size rewrite (mod_dimen.f90)
- Replaced single-point search with **flood-fill + bounding-box** algorithm
- Flood fill from `tpeak` location finds connected molten region at surface
- Sub-cell interpolation at bounding-box edges for length, width, depth
- Handles multi-pool (reports pool containing tpeak), cooling/turnaround
- `pool_size` moved outside iter_loop (once per timestep, not per iteration)
- `pool_mask` as module-level array to avoid repeated allocate/deallocate
- OpenMP on bounding-box and depth scans

### Fix: Thermal history interpolation (mod_print.f90)
- `write_thermal_history` now does bilinear interpolation at exact physical coordinates instead of reading nearest cell value
- Eliminates temperature jumps after AMR remesh

### Fix: AMR 20% coarse cell cap (mod_adaptive_mesh.f90)
- When melt pool grows too large, compute max `half_x/y` that satisfies `n_fine <= n_total/1.2` instead of reverting to previous (possibly also too large) value

### Fix: AMR efficiency optimizations (mod_adaptive_mesh.f90)
- Skip remesh when beam moves < 5×dx_fine (avoids unnecessary remesh during turnaround)
- Removed property field interpolation (den/vis/diff recomputed each iteration)
- Module-level interpolation buffer `amr_tmp` (avoid 15× allocate/deallocate per remesh)
- Removed `amr_validate_grid` from production loop
- AMR timing tracked in `t_amr` / `n_amr_remesh`

### Fix: VTK output cleanup
- Removed `solidID` and `localfield` from vtkmov output
- Single `defect.vtk` with scalars: `defect`, `maxtemp`
- `solidfield_def` tracked on uniform defect mesh each timestep (for AMR accuracy)
- Toolpath generator: added `--domain_x/y` for plot axis limits; generated `center_rot45.crs`

---

## Notes

- Mesh refinement region can be directly visualized in Paraview from the VTK mesh itself
- When `adaptive_flag=0`, behavior is identical to the original code (uniform mesh, no remeshing)
- Coefficient arrays (`ap`, `ae`, `aw`, etc. in `mod_coeff_data.f90`) do NOT need interpolation — they are recomputed each iteration from the geometry
- All execution results (Phase 1 through Phase 4) are recorded in `results.md`
