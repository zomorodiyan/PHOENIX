# Thermal-Mechanical — Execution Log

### [2026-04-04 11:41] Phase 1: Core module setup
- Created `mechanical/` subfolder
- Created `mod_mech_material.f90`: phase constants, material params, J2 return map, build_C_matrix
- Created `mod_mechanical.f90`: full EBE FEM solver ported from fortran-ebe-lpbf
- Created `mod_mech_io.f90`: VTK output, timing/memory reports, history, plot script
- Added `mechanical_flag`, `mech_interval`, `mech_output_interval` to mod_param.f90

### [2026-04-04 12:00] Phase 2: Integration
- Added `use mechanical_solver, use mech_io` to main.f90
- Mechanical solve called every mech_interval steps after iter_loop
- Reports res, max_vm, yield_elems to output.txt
- t_mech added to mod_timing.f90
- Separate mech timing/memory reports
- Updated compile.sh

### [2026-04-04 12:15] First compile — CLEAN BUILD

### [2026-04-04 12:18] Test run (200×200×52, mechanical_flag=1)
- Thermal solver runs fine
- Mechanical solver at step 10: CG never converges — grid too large (200×200×50 = 2M elements)
- Killed after 2+ minutes on a single mech solve

### [2026-04-04 12:22] Test run (50×50×52, mechanical_flag=1)
- Mech step 1 completed: res=3.6E+07, max_vm=2.1 GPa (too high)
- Mech step 2: CG still very slow
- Issues identified: tangent/residual inconsistency, Jacobi preconditioner insufficient

### [2026-04-04 12:30] Added mech_mesh_ratio parameter
- `mech_mesh_ratio=2` in `&mechanical_params` namelist
- Mechanical grid = thermal grid / ratio (e.g., 200/2 = 100 nodes per direction)
- FEM node coordinates stored in `fem_x, fem_y, fem_z`
- Temperature subsampled at every `mratio`-th cell center
- VTK output uses `fem_x/y/z` coordinates
- Compile: CLEAN BUILD

### [2026-04-04 12:35] Bug fix: per-element Ke in ebe_matvec and Jacobi preconditioner
- Root cause: ebe_matvec used precomputed Ke from reference element, but residual used per-element B with actual spacing → tangent/residual inconsistency → Newton divergence
- Fix: compute Ke per element using actual dx/dy/dz in both ebe_matvec and Jacobi preconditioner

### [2026-04-04 12:40] Bug fix: T_fem array size mismatch
- VTK/history was passing `temp(2:nim1,...)` (200×200×50) but FEM arrays are (51×51×13 with ratio=4)
- Fix: added `T_fem_last` public array in mechanical_solver, populated after each solve

### [2026-04-04 12:47] mech_debug test PASSED (50×50×12 FEM, ratio=4)
- 10 mech steps, all Newton converged (res ~0.02-0.03)
- max_vm: 30 → 127 MPa (realistic growth)
- yield_elems: 4 → 13
- 10 VTK files + history PNG generated
- Memory: 24 MB mechanical

### [2026-04-04 13:12] Single-track validation PASSED (single_track.crs, 200×200×52 thermal, 51×51×13 FEM)
- 10 mech steps, all converged
- max_vm: 30 → 127 MPa
- Mechanical: 90% of total CPU (2764s / 3077s)
- Temperature in history: 293K (correct, monitoring points away from track)
- Displacement: order 1e-7 m (reasonable)

### [2026-04-04 13:13] Multitrack attempt 1 (20 threads) — FAILED
- CG solver extremely slow with 20 threads on 51×51×13 grid (OpenMP overhead)
- 0 mech steps in 7 minutes, killed

### [2026-04-04 13:30] Fix: precompute per-Z-layer Ke
- Instead of computing 24×24 Ke at every CG iteration, precompute Ke_solid_z/Ke_soft_z for each Z-layer
- dx/dy uniform within layer, only dz varies

### [2026-04-04 13:39] Multitrack attempt 2 (4 threads) — SUCCESS
- 75 mech steps completed in ~22 min
- max_vm: 27 → 132 MPa (plateaus below yield)
- yield_elems: 4 → 54
- Mechanical: 24.6% of total CPU time
- 15 mech VTK files + history PNG generated

### [2026-04-04 14:02] Results written to results.md
