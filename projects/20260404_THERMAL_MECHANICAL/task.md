# Thermal-Mechanical Coupling (Residual Stress)

## Objective

Integrate the EBE (Element-By-Element) FEM mechanical solver from [fortran-ebe-lpbf](https://github.com/DrZGan/fortran-ebe-lpbf) into PHOENIX for residual stress computation. One-way coupling: PHOENIX temperature field drives thermal strain → displacement → stress via quasi-static equilibrium with J2 plasticity.

## Execution Tracking

- **`log.md`**: Execution log with system timestamps.
- **`results.md`**: Single-track and multi-track results.

---

## Source Reference: fortran-ebe-lpbf

Key files to port:
- `mod_mechanical.f90` — EBE CG solver, Newton-Raphson, J2 return map, 8-color coloring
- `mod_phase.f90` — Phase tracking (POWDER/LIQUID/SOLID)
- `mod_io.f90` — VTK output (only the mechanical fields part)

NOT ported:
- `mod_thermal.f90` — PHOENIX already handles temperature
- `mod_parameters.f90` — Parameters defined in PHOENIX's input system + module headers

---

## Architecture

```
fortran_new/
├── mechanical/                    # New subfolder for mechanical solver
│   ├── mod_mechanical.f90         # EBE FEM solver (main module)
│   ├── mod_mech_material.f90      # Phase-dependent material properties
│   └── mod_mech_io.f90            # Mechanical VTK output
├── main.f90                       # Calls mechanical solver
├── mod_param.f90                  # mechanical_flag in &output_control
└── inputfile/input_param.txt      # mechanical_flag=0/1
```

---

## Phase 1: Core Module Setup

### Task 1.1: Add mechanical parameters to input system

- **File**: `mod_param.f90`
- Add variables:
  ```fortran
  integer :: mechanical_flag = 0      ! 0=off, 1=on
  integer :: mech_interval = 10       ! solve mechanical every N thermal steps
  integer :: mech_output_interval = 5  ! output mech VTK every N mechanical solves
  ```
- Add new namelist:
  ```fortran
  namelist / mechanical_params / mechanical_flag, mech_interval, mech_output_interval
  ```
- Read in `read_data` (after `&adaptive_mesh`)
- **File**: `input_param.txt`
- Add:
  ```
  &mechanical_params mechanical_flag=0, mech_interval=10, mech_output_interval=5 /
  ```

### Task 1.2: Create `mechanical/mod_mech_material.f90`

Module `mech_material` — phase-dependent mechanical properties.

```fortran
module mech_material
    use precision
    implicit none

    ! Phase constants
    integer, parameter :: MECH_POWDER = 0
    integer, parameter :: MECH_LIQUID = 1
    integer, parameter :: MECH_SOLID  = 2

    ! Material parameters (defined here, not in inputfile)
    real(wp), parameter :: E_solid   = 70.0e9_wp    ! Young's modulus, solid (Pa)
    real(wp), parameter :: E_soft    = 0.7e9_wp     ! Young's modulus, powder/liquid (Pa)
    real(wp), parameter :: nu_mech   = 0.3_wp       ! Poisson's ratio
    real(wp), parameter :: sig_yield = 250.0e6_wp   ! J2 yield stress (Pa)
    real(wp), parameter :: alpha_V   = 1.0e-5_wp    ! Volumetric thermal expansion (1/K)

    ! Solver parameters
    integer, parameter  :: cg_maxiter_mech = 20000
    real(wp), parameter :: cg_tol_mech     = 1.0e-4_wp
    integer, parameter  :: newton_maxiter  = 10
    real(wp), parameter :: newton_tol      = 1.0e-4_wp

    contains
    ! Phase determination from PHOENIX fields
    ! ...
end module mech_material
```

**Phase mapping from PHOENIX**:
- `MECH_POWDER`: `solidfield(i,j,k) <= powder_threshold` AND `temp(i,j,k) < tsolid`
- `MECH_LIQUID`: `temp(i,j,k) >= tsolid` (currently molten, in melt pool)
- `MECH_SOLID`: `solidfield(i,j,k) > powder_threshold` AND `temp(i,j,k) < tsolid` (solidified, was once melted)

Subroutine `update_mech_phase(phase, temp, solidfield)`:
- Updates `phase(Nnx, Nny, Nnz)` from PHOENIX's `temp` and `solidfield`
- Uses PHOENIX's node-centered temperature and solidfield directly

### Task 1.3: Create `mechanical/mod_mechanical.f90`

Port from `fortran-ebe-lpbf/mod_mechanical.f90`. Key changes for PHOENIX integration:

**Grid mapping**:
- PHOENIX grid: `ni × nj × nk` with `nim1 = ni-1` interior cells. Cell-centered fields.
- EBE FEM: node-centered, `Nnx × Nny × Nnz` nodes, `Nx × Ny × Nz` elements.
- Mapping: `Nx = nim1-1`, `Ny = njm1-1`, `Nz = nkm1-1`, `Nnx = nim1`, `Nny = njm1`, `Nnz = nkm1`
- PHOENIX's `x(2:nim1)` = node positions for FEM (interior cell centers = FEM nodes)
- Note: PHOENIX cell centers are FEM node positions. The FEM element spans between adjacent cell centers.

**AMR compatibility**:
- Mechanical solver runs on the **same mesh as temperature** (simulation mesh). No separate mesh — direct field mapping, no interpolation at solve time.
- With AMR, dx/dy are non-uniform → **Ke must be computed per element** using actual node spacing `dx(ie) = x(ie+1) - x(ie)`, not a single precomputed Ke. The Jacobian `J` and `det_J` vary per element.
- **When AMR remeshes**: mechanical fields (`ux, uy, uz, sig_gp, eps_gp, T_old_for_u, phase`) must be **interpolated** from old mesh to new mesh, same as the thermal fields in `amr_interpolate_all_fields`. Add mechanical field interpolation to the AMR remesh routine.
- `sig_gp` and `eps_gp` are element-centered (at Gauss points) → interpolation requires mapping old element → new element, not just node interpolation. Strategy: interpolate node-averaged stress/strain, then distribute back to GPs. Or simpler: **reset GP state after remesh** (set `sig_gp=0, eps_gp=0`) and accept a small stress discontinuity — acceptable since remesh happens every 10-20 steps and stress evolves slowly.
- `ux, uy, uz` are node-centered → use the same bilinear X-Y interpolation as other fields.

**Non-uniform dz in Z**:
- PHOENIX has non-uniform Z spacing (2 zones: 0.5mm with 10 cells + 0.2mm with 40 cells).
- The EBE solver assumes uniform spacing → need to either:
  (a) Precompute Ke per element using actual dz(k) — more complex but exact
  (b) Use average dz — approximate but simple
  (c) Require uniform Z for mechanical solver — simplest
- **Decision**: Option (a) — precompute B matrices using actual node spacing. The Jacobian `J` at each GP depends on `dx(ie), dy(je), dz(ke)` of the element. For non-uniform Z, each layer has a different `dz`, so precompute `det_J` and `B` per unique `dz` value.

**Key subroutines to port** (from fortran-ebe-lpbf `mod_mechanical.f90`):
- `init_mechanical()` → allocate state arrays, precompute Ke/B/N
- `solve_mechanical(ux,uy,uz,T_new,T_old,phase,res_out,newton_iters,cg_iters)` → Newton + CG, returns final residual and iteration counts
- `compute_residual(...)` → 8-color EBE internal force
- `ebe_matvec_mech(...)` → 8-color EBE matrix-vector product
- `solve_mech_cg(...)` → Jacobi-preconditioned CG
- `j2_return_map(...)` → von Mises radial return
- `get_stress_yield(sxx,syy,szz,von_mises,fplus,phase)` → smooth GP stresses to nodes, compute von Mises
- `update_gp_state(...)` → update stress/strain history

**OpenMP**: Keep the 8-color coloring strategy from the original code.

**Boundary conditions**:
- Bottom face (k=1): u = 0 (Dirichlet, clamped substrate base)
- All other faces: traction-free (natural BC)

### Task 1.4: Create `mechanical/mod_mech_io.f90`

Module `mech_io` — separate VTK output for mechanical results.

Output file: `<case_name>_mech_NNNNN.vtk` (VTK legacy structured grid, binary)

**Fields** (all node-centered, on the mechanical grid):

| Field | Type | Description |
|-------|------|-------------|
| `Temperature` | scalar | T at mechanical nodes |
| `ux` | scalar | x-displacement |
| `uy` | scalar | y-displacement |
| `uz` | scalar | z-displacement |
| `phase` | integer scalar | 0=powder, 1=liquid, 2=solid |
| `sxx` | scalar | Smoothed stress xx component |
| `syy` | scalar | Smoothed stress yy component |
| `szz` | scalar | Smoothed stress zz component |
| `von_mises` | scalar | Von Mises equivalent stress |
| `fplus` | scalar | Yield function (>0 means plastic) |

**Output frequency**: Every `mech_interval` timesteps (when mechanical solve happens).

---

## Phase 2: Integration into main.f90

### Task 2.1: Call mechanical solver from main.f90

```fortran
use mech_material
use mechanical_solver
use mech_io

! After allocations:
if (mechanical_flag == 1) call init_mechanical()
integer :: mech_solve_count = 0

! In time loop, after iter_loop and field updates:
if (mechanical_flag == 1 .and. mod(step_idx, mech_interval) == 0) then
    call cpu_time(t0)
    call update_mech_phase(mech_phase, temp, solidfield)
    call solve_mechanical(ux_mech, uy_mech, uz_mech, temp, temp_old_mech, mech_phase, &
        mech_res, newton_iters, cg_iters_total)
    call get_stress_yield(sxx_out, syy_out, szz_out, vm_out, fplus_out, mech_phase)
    mech_solve_count = mech_solve_count + 1
    temp_old_mech = temp

    ! Report to output.txt (unit 9)
    n_yield = count(fplus_out > 0.0_wp)
    write(9,'(A,I6,A,es10.3,A,es10.3,A,I8)') &
        '  Mech step', mech_solve_count, &
        '  res=', mech_res, '  max_vm=', maxval(vm_out), &
        '  yield_elems=', n_yield

    ! VTK output at mech_output_interval
    if (mod(mech_solve_count, mech_output_interval) == 0) then
        call write_mech_vtk(step_idx, temp, ux_mech, uy_mech, uz_mech, mech_phase, &
            sxx_out, syy_out, szz_out, vm_out, fplus_out)
    endif

    call cpu_time(t1)
    t_mech = t_mech + (t1 - t0)
    n_mech_solves = n_mech_solves + 1
endif
```

**Output to output.txt** (every mechanical solve):
```
  Mech step   1  res=1.234E-05  max_vm=2.500E+08  yield_elems=   12345
  Mech step   2  res=8.765E-06  max_vm=2.510E+08  yield_elems=   12890
```

### Task 2.2: Mechanical timing and memory reports (separate from thermal)

Mechanical solver has its own timing and memory reports — does NOT modify `mod_timing.f90` or the existing thermal reports.

**File**: `mechanical/mod_mech_io.f90` — add `write_mech_timing_report()` and `write_mech_memory_report()`

**Timing report**: `<case_name>_mech_timing_report.txt`
```
============================================
  PHOENIX Mechanical Solver Timing Report
============================================
  Solves performed:       46
  Total CPU time:       1234.567 s
  Total wall time:       123.456 s
--------------------------------------------
  Avg per solve:          26.838 s
  Avg Newton iters/solve:   4.2
  Avg CG iters/Newton:    850
--------------------------------------------
  Phase update:           12.345 s (1.0%)
  Residual assembly:     456.789 s (37.0%)
  CG solver:             678.901 s (55.0%)
  Stress smoothing:       23.456 s (1.9%)
  VTK output:             45.678 s (3.7%)
  Other:                  17.398 s (1.4%)
============================================
```

**Memory report**: `<case_name>_mech_memory_report.txt`
```
============================================
  PHOENIX Mechanical Memory Report
============================================
  Grid: Nx × Ny × Nz elements
  GP state arrays:
    sig_gp(6,8,Nx,Ny,Nz):  XXXX.X MB
    eps_gp(6,8,Nx,Ny,Nz):  XXXX.X MB
  Displacement fields:
    ux,uy,uz(Nnx,Nny,Nnz): XXXX.X MB
  Stress output:
    sxx,syy,szz,vm,fplus:  XXXX.X MB
  Total mechanical:         XXXX.X MB
============================================
```

**In main.f90**: track `t_mech` locally (not in mod_timing), pass to `write_mech_timing_report()` after time loop.

**`mod_timing.f90`**: Add `real(wp) :: t_mech = 0.0_wp` and `integer :: n_mech_solves = 0`. When `mechanical_flag == 1`, include a single `mechanical` entry in the module timing table (shows total CPU time and percentage, like `adaptive mesh` or `defect`). The detailed mechanical breakdown (Newton/CG/residual) goes in the separate `_mech_timing_report.txt`.

### Task 2.3: Update compile.sh

Add mechanical source files:
```bash
mechanical/mod_mech_material.f90 \
mechanical/mod_mechanical.f90 \
mechanical/mod_mech_io.f90 \
```

---

## Phase 3: Single-Track Validation

### Task 3.1: Run single-track with mechanical solver

- Use existing input_param.txt geometry (do NOT modify grid size)
- `adaptive_flag=0`, `mechanical_flag=1`
- `toolpath_file='./ToolFiles/single_track.crs'`
- Verify:
  - VTK output contains all fields (T, ux, uy, uz, phase, sxx, syy, szz, von_mises, fplus)
  - Displacement field is physical (expansion near melt pool, zero at bottom)
  - Stress concentrations at solidification front
  - Phase field matches PHOENIX's solidfield

### Task 3.2: Validate against fortran-ebe-lpbf

- The mechanical solver model and all parameters must be **identical** to fortran-ebe-lpbf:
  - Isotropic linear elasticity + J2 plasticity with radial return
  - E_solid=70 GPa, E_soft=0.7 GPa, nu=0.3, sig_yield=250 MPa, alpha_V=1e-5
  - 8-node hex element, 2×2×2 Gauss quadrature
  - 8-color EBE coloring, Jacobi-preconditioned CG
  - Newton-Raphson with same tolerances (newton_tol=1e-4, cg_tol=1e-4)
  - Bottom face clamped (u=0), all other faces traction-free
  - Element phase: SOLID if any node is SOLID, else POWDER
  - Per-GP phase for stress update, element-level phase for CG tangent
- Qualitative sanity check: displacement/stress patterns should be physically reasonable (no quantitative comparison since thermal models differ)

---

## Phase 4: Multi-Track Extension

### Task 4.1: Run multi-track with center_rot0.crs

- Use existing input_param.txt geometry (do NOT modify grid size)
- `adaptive_flag=0`, `mechanical_flag=1`
- `toolpath_file='./ToolFiles/center_rot0.crs'`
- Check:
  - Stress accumulation across multiple tracks
  - Track-to-track reheating effects on stress relaxation
  - Phase field correctly transitions

### Task 4.2: Mechanical history at monitoring points

Same 10 monitoring points as thermal history (`thist_px, thist_py, pz` from `mod_print.f90`). Written every mechanical solve step.

**File**: `<case_name>_mech_history.txt`
```
# Mechanical History - 10 Monitoring Points
# Columns: time(s)  T1..T10(K)  ux1..ux10(m)  uy1..uy10(m)  uz1..uz10(m)  sxx1..sxx10(Pa)  syy1..syy10(Pa)
```

Each row = one mechanical solve. 60 columns: time + 10×(T, ux, uy, uz, sxx, syy) = 1 + 60.

**Plot**: `<case_name>_mech_history.png` (generated by `finalize_mechanical()`)

Python script generates a 3×1 subplot figure:
1. **Top**: Temperature (K) vs time at 10 points (same as thermal history but at mechanical solve intervals)
2. **Middle**: Displacement magnitude `sqrt(ux²+uy²+uz²)` (μm) vs time at 10 points
3. **Bottom**: Stress sxx (MPa) vs time at 10 points, with `sig_yield` horizontal line

**Implementation** in `mod_mech_io.f90`:
- `init_mech_history()` — open file, write header with physical coordinates
- `write_mech_history(t, T, ux, uy, uz, sxx, syy)` — append one row using bilinear interpolation at exact physical coordinates (same as `write_thermal_history`)
- `finalize_mech_history()` — write Python plot script and execute

### Task 4.3: Generate deformation animation

After simulation completes, generate `<case_name>_deformation.gif` using a Python script written by the Fortran code (same pattern as `plot_thermal_history.py` and `plot_meltpool.py`).

**Script**: `<case_name>_plot_deformation.py` (written by `finalize_mechanical()` in `mod_mech_io.f90`)

**Animation spec** (based on `fdm_animation.py` from fortran-ebe-lpbf):
- Read all `<case_name>_mech_NNNNN.vtk` files in order
- Warp mesh by displacement × magnification factor (default 10x)
- Color by phase: POWDER=blue(#3366FF), LIQUID=white, SOLID=red(#FF3333)
- Overlay text: step number, deformation magnification
- Isometric camera view, consistent across frames
- Output: GIF, 100ms per frame, loop=0
- Uses PyVista (`pyvista`, `PIL`)

**`finalize_mechanical()` subroutine** (called after time loop):
- Writes the Python script to result directory
- Executes it: `python3 <script_name>.py`
- Similar to `finalize_thermal_history()` and `finalize_meltpool_history()`

```fortran
subroutine finalize_mechanical()
    ! Write Python animation script
    open(unit=lun, file=trim(file_prefix)//'plot_deformation.py', status='replace')
    write(lun,'(a)') 'import pyvista as pv'
    write(lun,'(a)') 'from PIL import Image'
    write(lun,'(a)') 'import glob, os'
    write(lun,'(a)') 'pv.OFF_SCREEN = True'
    ! ... write full script ...
    write(lun,'(a)') 'frames[0].save("'//trim(case_name)//'_deformation.gif",'
    write(lun,'(a)') '    save_all=True, append_images=frames[1:], duration=100, loop=0)'
    close(lun)
    ! Execute
    call execute_command_line('cd '//trim(result_dir)//' && python3 ...')
end subroutine
```

### Task 4.3: Performance analysis

- Timing breakdown: thermal vs mechanical
- Memory usage: thermal vs mechanical arrays
- Scalability with grid size

---

## Key Design Decisions

1. **Separate subfolder** `mechanical/` keeps FEM code organized and independent from CFD modules
2. **Same mesh as thermal**: mechanical solver runs on the simulation mesh (same as temperature). With AMR, Ke is computed per element using actual non-uniform spacing. Mechanical fields are interpolated during AMR remesh (ux/uy/uz bilinear, GP state reset to zero).
3. **Parameters in module header**: `E_solid`, `nu`, `sig_yield`, etc. defined as `parameter` constants in `mod_mech_material.f90`, not in inputfile (reduces input complexity)
4. **One-way coupling**: T → stress only. No feedback from displacement to temperature.
5. **Solve interval**: every `mech_interval=10` timesteps (mechanical equilibrium is quasi-static, doesn't need every thermal step)
6. **Separate VTK**: mechanical results in `_mech_NNNNN.vtk`, not mixed with thermal vtkmov files
7. **GP state arrays**: `sig_gp` and `eps_gp` are essential for incremental plasticity — cannot be avoided. Memory scales as 96 × N_elements × 8 bytes.

## Future: Thermal-Mechanical Overlap (deferred)

Considered running thermal and mechanical concurrently on separate thread groups (e.g., 10+10 out of 20). Analysis:
- Requires `OMP_NESTED=TRUE`, `OMP_MAX_ACTIVE_LEVELS=2`
- Each solver gets half the threads → ~35% slower individually
- Only beneficial when mechanical accounts for >35% of total time
- **Decision**: implement serial first, measure timing, then decide

## Notes

- The EBE approach avoids assembling a global stiffness matrix → memory efficient for the solver itself
- The GP state arrays are the main memory cost — for large grids, consider solving on a coarser mechanical mesh
- The 8-color coloring for OpenMP is identical to standard structured-grid FEM parallelization
- J2 plasticity with radial return is the simplest metal plasticity model — sufficient for residual stress estimation
- `fplus > 0` indicates plastic yielding has occurred — useful for identifying plastically deformed regions
