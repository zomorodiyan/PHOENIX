# Solver Algorithm

## Main Program Call Flow

Detailed call graph showing which function belongs to which module (`.f90` file).

```
main.f90
│
│ ══════════════ INITIALIZATION ══════════════
│
├── read_data()                          [mod_param.f90]
├── read_toolpath()                      [mod_toolpath.f90]
├── generate_grid()                      [mod_geom.f90]
│     └── generate_1d_grid() ×3
├── allocate_fields(ni, nj, nk)          [mod_init.f90]
│     ├── allocate_field_data()          [mod_field_data.f90]
│     └── allocate_coeff_data()          [mod_coeff_data.f90]
├── allocate_source(ni, nj, nk)          [mod_sour.f90]
├── allocate_print(ni, nj, nk)           [mod_print.f90]
├── allocate_laser(ni, nj)               [mod_laser.f90]
├── allocate_defect(ni, nj, nk)          [mod_defect.f90]
├── OpenFiles()                          [mod_print.f90]
├── initialize()                         [mod_init.f90]
├── init_thermal_history()               [mod_print.f90]
├── init_meltpool_history()              [mod_print.f90]
│
├── [if species_flag == 1]
│     ├── read_species_params()          [species_solver/mod_species.f90]
│     ├── allocate_species()             [species_solver/mod_species.f90]
│     └── init_species()                 [species_solver/mod_species.f90]
├── [if mechanical_flag == 1]
│     ├── init_mechanical()              [mechanical/mod_mechanical.f90]
│     └── init_mech_history()            [mechanical/mod_mech_io.f90]
├── [if adaptive_flag == 1]
│     └── amr_init()                     [mod_adaptive_mesh.f90]
│
│ ══════════════ TIME STEPPING LOOP ══════════════
│
├── time_loop: do while (timet < timax)
│   │
│   ├── laser_beam()                     [mod_laser.f90]
│   ├── read_coordinates()               [mod_toolpath.f90]
│   ├── [if adaptive_flag == 1]
│   │     ├── amr_check_remesh(step_idx) [mod_adaptive_mesh.f90]
│   │     ├── amr_regenerate_grid()      [mod_adaptive_mesh.f90]
│   │     │     └── amr_interpolate_all_fields() [mod_adaptive_mesh.f90]
│   │     ├── amr_interp_field(concentration, conc_old) [if species_flag=1]
│   │     ├── update_mech_grid()         [if mechanical_flag=1, serial mode]
│   │     │     ├── mech_interp_node_nn()  ← nearest-neighbor for u, T_old, f_plus
│   │     │     ├── mech_interp_gp_field() ← nearest-neighbor for sig_gp
│   │     │     └── recompute eps_gp       ← B_new * u_interp
│   │     └── amr_validate_grid()        [mod_adaptive_mesh.f90]
│   │
│   │ ──────────── PREDICTION (optional) ────────────
│   │
│   ├── [if predict_flag==1 AND laser on AND tpeak > tsolid]
│   │     ├── predict_shift_integer()      [mod_predict.f90]
│   │     │     └── ishift_field() ×5      [mod_predict.f90]
│   │     └── enthalpy_to_temp()           [mod_entot.f90]
│   │
│   │ ──────────── ITERATION LOOP ────────────
│   │
│   ├── iter_loop: do while (niter < maxit)
│   │   │
│   │   │ ── ENERGY EQUATION ──
│   │   │
│   │   ├── properties()                 [mod_prop.f90]
│   │   ├── bound_enthalpy()             [mod_bound.f90]
│   │   ├── discretize_enthalpy()        [mod_discret.f90]
│   │   ├── source_enthalpy()            [mod_sour.f90]
│   │   ├── calc_enthalpy_residual()     [mod_resid.f90]
│   │   ├── enhance_converge_speed()     [mod_converge.f90]
│   │   ├── solution_enthalpy()          [mod_solve.f90]
│   │   │     └── tdma_solve_3d_2()      [mod_solve.f90]
│   │   ├── enthalpy_to_temp()           [mod_entot.f90]
│   │   ├── pool_size()                  [mod_dimen.f90]
│   │   │
│   │   │ ── MOMENTUM + PRESSURE (if melt pool exists) ──
│   │   │
│   │   ├── [if tpeak > tsolid]
│   │   │   │
│   │   │   ├── cleanuvw()               [mod_solve.f90]
│   │   │   │
│   │   │   ├── ── u-momentum ──
│   │   │   ├── bound_uv(1)              [mod_bound.f90]
│   │   │   ├── discretize_momentum(1)   [mod_discret.f90]
│   │   │   ├── source_momentum(1)       [mod_sour.f90]
│   │   │   ├── calc_momentum_residual() [mod_resid.f90]
│   │   │   ├── solution_uvw(uVel)       [mod_solve.f90]
│   │   │   │
│   │   │   ├── ── v-momentum ──
│   │   │   ├── bound_uv(2)              [mod_bound.f90]
│   │   │   ├── discretize_momentum(2)   [mod_discret.f90]
│   │   │   ├── source_momentum(2)       [mod_sour.f90]
│   │   │   ├── calc_momentum_residual() [mod_resid.f90]
│   │   │   ├── solution_uvw(vVel)       [mod_solve.f90]
│   │   │   │
│   │   │   ├── ── w-momentum ──
│   │   │   ├── bound_w()                [mod_bound.f90]
│   │   │   ├── discretize_momentum(3)   [mod_discret.f90]
│   │   │   ├── source_momentum(3)       [mod_sour.f90]
│   │   │   ├── calc_momentum_residual() [mod_resid.f90]
│   │   │   ├── solution_uvw(wVel)       [mod_solve.f90]
│   │   │   │
│   │   │   ├── ── Pressure correction ──
│   │   │   ├── bound_pp()               [mod_bound.f90]
│   │   │   ├── discretize_pp()          [mod_discret.f90]
│   │   │   ├── source_pp()              [mod_sour.f90]
│   │   │   ├── calc_pressure_residual() [mod_resid.f90]
│   │   │   ├── solution_uvw(pp)         [mod_solve.f90]
│   │   │   └── revision_p()             [mod_revise.f90]
│   │   │
│   │   │ ── CONVERGENCE CHECK ──
│   │   │
│   │   ├── heat_fluxes()                [mod_flux.f90]
│   │   └── [exit if converged]
│   │
│   │ ──────────── AFTER ITERATION LOOP ────────────
│   │
│   ├── [if species_flag == 1]
│   │   ├── species_bc()                 [species_solver/mod_species.f90]
│   │   └── solve_species()              [species_solver/mod_species.f90]
│   │         ├── [FVM discretization]   (inline, power-law scheme)
│   │         ├── [boundary transfers]   (inline)
│   │         ├── [assembly + URF]       (inline)
│   │         ├── solution_species_tdma()[species_solver/mod_species.f90]
│   │         ├── enhance_species_speed()[species_solver/mod_species.f90]
│   │         ├── [concentration clip]   (inline, [0,1])
│   │         └── calc_species_residual()[species_solver/mod_species.f90]
│   │
│   ├── update_max_temp()                [mod_defect.f90]
│   │
│   ├── [if mechanical_flag==1 AND mod(step, mech_interval)==0]
│   │   ├── solve_mechanical()           [mechanical/mod_mechanical.f90]
│   │   │     ├── extract_temp_to_fem()    ← T from PHOENIX to FEM grid
│   │   │     ├── update_mech_phase()      [mechanical/mod_mech_material.f90]
│   │   │     ├── compute_dT_gp()          ← dT at Gauss points
│   │   │     ├── Newton loop:
│   │   │     │   ├── compute_residual()   ← 8-color EBE assembly
│   │   │     │   └── solve_mech_cg()      ← Jacobi-preconditioned CG
│   │   │     │         └── ebe_matvec_mech() ← matrix-free or precomputed Ke
│   │   │     └── update_gp_state()        ← J2 return map at each GP
│   │   ├── get_stress_yield()           [mechanical/mod_mechanical.f90]
│   │   ├── write_mech_vtk()            [mechanical/mod_mech_io.f90]
│   │   └── write_mech_history()         [mechanical/mod_mech_io.f90]
│   │
│   ├── CalTime()                        [mod_print.f90]
│   ├── outputres()                      [mod_print.f90]
│   │
│   ├── [velocity zeroing + field update] (inline in main.f90)
│   ├── [if species_flag == 1]
│   │   └── conc_old = concentration     (inline in main.f90)
│   ├── Cust_Out()                       [mod_print.f90]
│   │     ├── write_vtk_vector()         [mod_print.f90]
│   │     └── write_vtk_scalar() ×8-10   [mod_print.f90]
│   ├── write_thermal_history()          [mod_print.f90]
│   └── write_meltpool_history()         [mod_print.f90]
│
│ ══════════════ POST-SIMULATION ══════════════
│
├── compute_defect_determ()              [mod_defect.f90]
├── write_defect_report()                [mod_defect.f90]
│     └── write_defect_vtk() ×2         [mod_defect.f90]
├── EndTime()                            [mod_print.f90]
├── finalize_thermal_history()           [mod_print.f90]
├── finalize_meltpool_history()          [mod_print.f90]
├── [if mechanical_flag == 1]
│     ├── write_mech_timing_report()     [mechanical/mod_mech_io.f90]
│     ├── write_mech_memory_report()     [mechanical/mod_mech_io.f90]
│     ├── finalize_mech_history()        [mechanical/mod_mech_io.f90]
│     ├── finalize_mechanical_io()       [mechanical/mod_mech_io.f90]
│     │     └── generates plot_deformation.py → deformation.gif
│     └── cleanup_mechanical()           [mechanical/mod_mechanical.f90]
├── write_timing_report()                [mod_timing.f90]
└── write_memory_report()                [mod_timing.f90]
```

## SIMPLE Algorithm

Within each iteration, the solver uses the SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm:

1. **Guess** pressure field $p^*$
2. **Solve momentum** (u, v, w) with guessed pressure → $u^*, v^*, w^*$
3. **Solve pressure correction** $p'$ from mass continuity
4. **Correct** velocities: $u = u^* + d_u(p'_W - p'_P)$
5. **Update** pressure: $p = p^* + \alpha_p \cdot p'$
6. **Repeat** until convergence

The enthalpy equation is solved first (before momentum) because material properties depend on temperature.

## Convergence Criteria

| Condition | Criterion |
|-----------|-----------|
| Laser on (heating) | `resorh < 1e-5` AND `0.99 < ratio < 1.01` |
| Laser off (cooling) | `resorh < 1e-6` |
| Max iterations | `niter >= maxit` (forced exit) |

Where `ratio` is the energy balance ratio from `heat_fluxes()`.

## Staggered Grid

```
        j+1  ───────────────────
             |                 |
             |    vVel(i,j+1)  |
             |        ↑        |
             |        |        |
    j    ────┤  uVel ←  P,T,H → uVel(i+1)
             |  (i)   |  (i,j)    |
             |        ↓        |
             |    vVel(i,j)    |
             |                 |
        j-1  ───────────────────
             i-1      i      i+1
```

- Scalar variables (P, T, H, C, fracl) at cell centers
- u-velocity at east/west faces (i±½)
- v-velocity at north/south faces (j±½)
- w-velocity at top/bottom faces (k±½)
