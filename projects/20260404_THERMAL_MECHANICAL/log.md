# Thermal-Mechanical — Execution Log

### [2026-04-04 11:41:26] Starting Phase 1: Core module setup
- Created `mechanical/` subfolder
- Created `mod_mech_material.f90`: phase constants, material params, J2 return map, build_C_matrix
- Created `mod_mechanical.f90`: full EBE FEM solver ported from fortran-ebe-lpbf
  - Per-element B matrix computation (supports non-uniform grid)
  - 8-color EBE matvec, Jacobi-preconditioned CG
  - Newton-Raphson with J2 plasticity
  - Temperature extraction from PHOENIX grid to FEM nodes
- Created `mod_mech_io.f90`: VTK output, timing/memory reports, history, plot script generation
- Added `mechanical_flag`, `mech_interval`, `mech_output_interval` to mod_param.f90 + input_param.txt

### [2026-04-04 12:30:00] Phase 2: Integration
- Added `use mechanical_solver, use mech_io` to main.f90
- Mechanical solve called after iter_loop, every mech_interval steps
- Reports to output.txt: res, max_vm, yield_elems
- t_mech added to mod_timing.f90 (entry in thermal timing table)
- Separate mech timing/memory reports written by mod_mech_io
- Updated compile.sh with mechanical/*.f90 files

### [2026-04-04 12:45:00] First compile — CLEAN BUILD

### [2026-04-04 12:50:00] First test (50×50×52, single_track.crs, mechanical_flag=1)
- Thermal solver runs fine
- Mechanical step 1 completed: res=3.6E+07, max_vm=2.1 GPa (too high)
- Mechanical step 2: CG solver extremely slow (>2 min without converging)

### Known issues to fix:
1. **CG convergence**: Jacobi preconditioner insufficient for this grid size. May need multigrid or SSOR preconditioner.
2. **Stress values**: max_vm=2.1 GPa >> sig_yield=250 MPa. Check Newton residual interpretation — `res` should be relative but showing 3.6E+07.
3. **Non-uniform Z grid**: The precomputed Ke_solid/Ke_soft used in ebe_matvec_mech assume uniform dx/dy/dz (line `dx_ref = x(3)-x(2)`). But compute_residual uses per-element B matrix with actual spacing. This inconsistency between tangent (CG) and residual may hurt convergence.
4. **Performance**: 50×50×50 = 125K elements. Each CG iteration = 8-color matvec over all elements. With 20K max CG iterations × 10 Newton steps = 200K matvec operations per solve. Too slow.
