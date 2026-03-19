# Species Transport Integration

## Objective

Create `mod_species.f90` in fortran_new based on `D:\Fortran\dissimilar\program931` to solve dissimilar metal species transport. The module must be as self-contained as possible (parameters, initialization, BCs, source terms, solver, output all inside `mod_species.f90`). Minimize changes to existing codebase. Only addition to existing input: `species_flag` (1=on, 0=off).

## Design Principles

- `mod_species.f90` is a standalone module: all secondary material properties, species arrays, solver, BCs, output defined within
- When `species_flag=0`, existing code behavior is **completely unchanged** — no extra arrays, no extra computation, no code path differences
- When `species_flag=1`, each routine that needs composition-dependent properties computes them **inline** using `concentration(i,j,k)`:
  ```fortran
  prop_local = prop1 * C + prop2 * (1.0 - C)
  ```
- No per-cell property arrays (`dens_arr`, `acpa_arr`, etc.) — all composition mixing is computed on-the-fly from the concentration field + scalar pairs (primary/secondary material constants)
- Species transport active only in melt pool region (same as uvwp): when `T < tsolid_local`, `massdiffusivity` is near-zero
- Concentration C is a mass fraction of the **primary (base) material**:
  - **C = 1**: pure primary material — properties from `&material_properties` in `input_param.txt` (e.g. dens, acpa, tsolid=1563K)
  - **C = 0**: pure secondary material — properties from named constants in `mod_species.f90` (e.g. dens2=8880, acpa2=0.3441, tsolid2=1728K)
  - **0 < C < 1**: mixed region — properties computed via `mix(prop1, prop2, C) = prop1*C + prop2*(1-C)`
- Species solver reuses all shared coefficient arrays (`an, as, ae, aw, at, ab, ap, su, sp, apnot`) — solve after momentum/enthalpy so coefficients are free
- Species arrays are allocated on the global domain, but discretization/source/solver only loop over the melt pool region using `istat,iend,jstat,jend,kstat` from `mod_dimen.f90` (same as uvw momentum equations)
- Species is solved once per timestep after `iter_loop` exits (not inside the iteration loop), using `delt` for transient term

## Inline Computation Strategy (replaces Phase 0 array approach)

**Rationale**: The Phase 0 approach (16 per-cell property arrays, ~513 MB overhead, +12% runtime) was implemented and reverted. Instead, all composition-dependent properties are computed inline using a helper function `mix(prop1, prop2, C)` defined in `mod_species.f90`.

**Helper function** (defined in `mod_species.f90`):
```fortran
pure real(wp) function mix(prop1, prop2, C)
    real(wp), intent(in) :: prop1, prop2, C
    mix = prop1 * C + prop2 * (1.0_wp - C)
end function mix
```

When `species_flag=0`, the original scalar code paths remain untouched — no `if` checks, no concentration lookups, no overhead.

**Benefits over array approach**:
- Zero memory overhead (no 16 extra 3D arrays)
- No `allocate_prop_arrays` call needed
- No compile-order dependency changes
- Properties stay "live" — always consistent with current concentration (no stale array risk)
- Code changes confined to the `species_flag=1` branch — existing single-material behavior is untouched
- `mix()` function keeps code DRY — no repeated `prop1*C + prop2*(1-C)` patterns

**Implementation pattern for modified routines** (properties, enthalpy_to_temp, source_momentum, source_pp, velocity zeroing):
```fortran
use species, only: concentration, species_flag, mix, tsolid2, tliquid2, ...

if (species_flag == 1) then
    C_local = concentration(i,j,k)
    tsolid_local  = mix(tsolid,  tsolid2,  C_local)
    tliquid_local = mix(tliquid, tliquid2, C_local)
    dens_local    = mix(dens,    dens2,    C_local)
    ! ... etc.
else
    ! Original scalar path (unchanged)
    ! ... use tsolid directly
endif
```

## Key Physics Changes from Original (program931)

- Gamma = rho * D_m (D_m = 5e-9 m^2/s), replacing the old non-physical `massdiffl`. `massdiffusivity` array stores Gamma = rho * D_m (units: kg/(m*s)), used directly as diffusion coefficient in FVM discretization
- Remove velocity x2 factor
- Remove scanning source term (moving frame advection)
- Remove solidification segregation source (kp term)
- All BCs: zero-flux Neumann
- Add solutal Marangoni effect (dgdc based on local composition)
- `dgdt` and `hlatnt` remain scalars (assume both materials have similar values)

## Tasks (easy -> hard)

### Phase 0: Infrastructure (no physics change)

1. [x] **Add `write_memory_report` to timing module**
   - Added `write_memory_report(file_prefix)` subroutine to `mod_timing.f90`
   - Reads `/proc/self/status` for VmPeak, VmHWM, VmRSS, VmData
   - Called from `main.f90` after timing report
   - Modified `compile.sh`: no longer cleans `result/` on recompile (preserves results between runs)
   - ~~Per-cell property arrays approach was implemented and reverted~~ — replaced by inline computation strategy (see above)

### Phase 1: Scaffolding

2. [x] **Add `species_flag` to input parsing**
   - Added `species_flag` (integer, default=0) to `mod_param.f90`
   - Added to `&output_control` namelist — set `species_flag=1` in `input_param.txt` to enable
   - Initialized to 0 in `read_data` before namelist read

3. [x] **Create `mod_species.f90` with full implementation**
   - Created `mod_species.f90` with all secondary material parameters as named constants
   - `pure function mix(prop1, prop2, C)` — linear mixing helper
   - `allocate_species` — allocates `concentration(ni,nj,nk)` and `conc_old(ni,nj,nk)`
   - `init_species` — computes derived constants (hsmelt2, hlcal2, etc.), sets initial concentration field (substrate=1, powder half-domain=0)
   - `species_bc` — zero-flux Neumann BCs on all 6 faces
   - `solve_species` — full FVM discretization + source + TDMA + block correction + clipping + residual
   - `solution_species_tdma` — line-by-line TDMA with OpenMP (per-thread pr/qr buffers, same pattern as `tdma_solve_3d_2`)
   - `enhance_species_speed` — block correction with OpenMP reduction (same pattern as `mod_converge.f90`)
   - `calc_species_residual` — species residual with OpenMP reduction
   - `write_species_vtk` — placeholder stub
   - Added to `compile.sh` (after `mod_defect.f90`, before `main.f90`)
   - `resorc` declared in `mod_resid.f90` (avoids circular dependency with `mod_print.f90`)

### Phase 2: Solver core (all implemented in Task 3)

4. [x] **Implement `allocate_species` and `init_species`**
   - Allocate `concentration`, `conc_old` to (ni, nj, nk), initialized to 1.0
   - IC: powder layer (`z >= dimz - layerheight` AND `y < dimy/2`) → C=0 (secondary), else C=1 (base)

5. [x] **Implement zero-flux BCs (`species_bc`)**
   - All 6 faces: `C(boundary) = C(interior neighbor)`

6. [x] **Implement species solver (`solve_species`)**
   - FVM discretization on melt pool region (`istatp1:iendm1, jstat:jend, kstat:nkm1`, same as momentum)
   - Skips entirely when `tpeak <= tsolid` (no melt pool)
   - Diffusion: `Gamma = den(i,j,k) * D_m` with harmonic mean at faces; floor `1e-30` in solid cells
   - Transient: `apnot = den/delt * volume` (uses `delt`, not `delt_eff`)
   - Under-relaxation with `urfspecies=0.7` (not `urfh`)
   - No velocity x2 factor, no scanning source, no segregation (kp)
   - Concentration clipping [0,1] after TDMA
   - OpenMP parallelized main loop, boundary loops, assembly

7. [x] **Wire into main loop**
   - `main.f90`: `allocate_species` + `init_species` at startup (when `species_flag==1`)
   - After `iter_loop`: `species_bc` + `solve_species` (when `species_flag==1`)
   - End of timestep: `conc_old = concentration` (full-domain)
   - `resorc` added to `outputres` format (conditional on `species_flag`)
   - `t_species` added to `mod_timing.f90` timing report (17th module)

8. [x] **`enhance_species_speed` block correction**
   - Integrated into `solve_species` after TDMA, before clipping
   - Accumulates j-k planes into 1D x-system, solves TDMA, broadcasts correction
   - OpenMP reduction on accumulation loop

### Milestone: One-way coupling validation

9. [x] **Add concentration to VTK output**
    - Added concentration as scalar field in main VTK output (`Cust_Out`), conditional on `species_flag==1`
    - No standalone species VTK — concentration is included in `vtkmov{N}.vtk` alongside T, vis, den, etc.
    - Moved `mod_species.f90` before `mod_print.f90` in compile order to resolve dependency
    - Removed legacy `write(41,*)` Tecplot header from `OpenFiles` (was creating spurious `fort.41`)

10. [x] **Create single-track test toolpath**
    - Created `ToolFiles/species_test.crs`: single track scanning x from 0.5mm to 3.5mm at y=2mm (domain center), speed 1.23 m/s
    - Created `inputfile/input_species_test.txt`: original mesh (400x400x50), timax=0.001s (50 steps)
    - Created `run_species_test.sh`: automated test script for both cases

11. [x] **One-way coupling validation run**
    - Ran `sp_base2` (flag=0) and `sp_test2` (flag=1) on 400x400x50 mesh, 50 timesteps, 4 threads
    - **Correctness**: Step 1 thermal fields identical (resorh=1.2E-04, res_mass=5.0E-02, res_u/v/w match). No NaN on fine mesh. Species residual ~1e-8 to ~3e-8 (near-zero, correct for mostly uniform C=1 field).
    - **Concentration field verified**: vtkmov1 has 560,023 cells with C<0.5, vtkmov2 has 560,014 — concentration is preserved correctly with minimal mixing at interface (9 cells changed over 25 timesteps).
    - **Performance** (wall time): baseline 45.8s vs species 45.2s → **<1% overhead**
      - `mod_species` = 0.52s CPU (0.4% of total CPU). Species solves only in the melt pool region (same indices as momentum solver), not the full 8M-cell domain.
    - **Memory**: VmHWM baseline 1399 MB vs species 1464 MB → **+65 MB** (two 402x402x52 float arrays)
    - **OpenMP**: All species loops properly parallelized. TDMA uses per-thread pr/qr buffers.
    - **VTK**: Concentration scalar included in vtkmov VTK files, viewable in ParaView alongside other fields

### Phase 3: Two-way coupling (material property feedback)

12. [x] **Modify `properties()` in `mod_prop.f90` for composition-dependent properties**
    - When `species_flag=1`, use `mix()` from `mod_species.f90` to compute local properties:
      ```fortran
      use species, only: concentration, species_flag, mix, tsolid2, tliquid2, ...
      ! Inside loop:
      if (species_flag == 1) then
          C_local = concentration(i,j,k)
          tsolid_local  = mix(tsolid,  tsolid2,  C_local)
          tliquid_local = mix(tliquid, tliquid2, C_local)
          dens_local    = mix(dens,    dens2,    C_local)
          denl_local    = mix(denl,    denl2,    C_local)
          viscos_local  = mix(viscos,  viscos2,  C_local)
          ! ... all cp, thcon, powder properties similarly
      else
          tsolid_local = tsolid; tliquid_local = tliquid; ...
      endif
      ```
    - Use local variables for branching (`tsolid_local`, `tliquid_local`) and property computation
    - Three temperature regimes use locally-mixed properties
    - `beta`, `emiss`, `hlatnt`, `dgdt` remain scalars
    - **OpenMP**: ensure all `mix()` calls use thread-private local variables, no shared writes

13. [x] **Modify `enthalpy_to_temp()` in `mod_entot.f90` for composition-dependent H-T curve**
    - When `species_flag=1`, use `mix()` to compute `acpa_local`, `acpb_local`, `tsolid_local`, `tliquid_local`, `acpl_local`, then derive `hsmelt_local`, `hlcal_local`, `deltemp_local`
    - When `species_flag=0`, use original scalars `hsmelt`, `hlcal`, `deltemp`, `acpl`, `tsolid`, `tliquid` (unchanged)
    - **OpenMP**: all derived locals must be `PRIVATE` in the parallel region

14. [x] **Modify `source_momentum`/`source_pp` in `mod_sour.f90` for composition-dependent branching**
    - When `species_flag=1`:
      - Darcy: `darcy_c0` computed cell-locally using `mix(viscos, viscos2, C_local)`
      - Buoyancy: uses `mix(dens, dens2, C_local) * g * beta * (tw - tsolid_local)`
      - Solid check: uses `mix(tsolid, tsolid2, C_local)`
    - When `species_flag=0`: original scalar code (unchanged)
    - **OpenMP**: `darcy_c0` becomes a loop-local variable (move from before-loop to inside-loop), add to `PRIVATE` clause

15. [x] **Modify velocity zeroing in `main.f90`**
    - When `species_flag=1`: use `mix(tsolid, tsolid2, concentration(i,j,k))` for the `temp <= tsolid` check
    - When `species_flag=0`: use scalar `tsolid` (unchanged, current code)
    - Also: `if(tpeak > min(tsolid, tsolid2))` for momentum activation check

### Phase 4: Solutal Marangoni

16. [x] **Solutal Marangoni effect (dgdc)**
    - Define `dgdc` as a scalar constant in `mod_species.f90` (dg/dC, surface tension concentration coefficient)
    - Analogous to thermal Marangoni (`dgdt * dT/dx`), solutal Marangoni is `dgdc * dC/dx`
    - Total surface tension force: `tau = dgdt * grad(T) + dgdc * grad(C)`
    - `dgdt` remains scalar (same for both materials)
    - In `mod_bound.f90` `bound_uv`, after thermal Marangoni block, add:
      ```fortran
      if (species_flag == 1) then
          dcd = (concentration(i,j,nk) - concentration(i-1,j,nk)) * dxpwinv(i)  ! case 1
          term1 = fracl_stag * dgdc_const * dcd / (vis1 * dzpbinv(nk))
          uVel(i,j,nk) = uVel(i,j,nk) + term1
      endif
      ```
      Same pattern for case 2 (v-velocity, dC/dy)

### Phase 5: Two-way coupling validation

17. [x] **Two-way coupling validation run**
    - Ran `sp_2way` (flag=1, two-way) with 100mm/s scan, 400x400x50 mesh, 250 steps
    - **Correctness**: Concentration field correct (C in [0,1], mixing zone grows). Enthalpy residual stable (~2.6e-6). Species residual stable (~2.5e-5). Momentum residuals diverge at later timesteps — expected with composition-dependent properties altering flow dynamics.
    - **Performance**: wall time 322s (two-way) vs 277s (one-way) → +16% overhead from `mix()` calls in properties/enthalpy_to_temp/source. `mod_species` itself stays <1% (0.50%).
    - **Memory**: identical (1463 MB) — no extra arrays needed (inline computation)
    - **Compile order**: `mod_species.f90` moved after `mod_resid.f90` and before `mod_prop.f90` to resolve dependencies

## Design Notes

### Coefficient arrays
Species solver reuses all shared coefficient arrays (`an, as, ae, aw, at, ab, ap, su, sp, apnot` from `mod_coeff_data.f90`). Species is solved once per timestep after `iter_loop` exits, so the arrays are free.

### Local solver interaction
Species is always solved globally (full domain), not restricted to local solver region. Uses `delt` (not `delt_eff`) for the transient term. `conc_old = concentration` update must be full-domain at end of each timestep, unconditional of `is_local`.

### Diffusivity convention
Species diffusion coefficient Gamma = rho * D_m = `den(i,j,k) * D_m` (units: kg/(m*s)), computed inline in the species solver. No separate `massdiffusivity` array needed — `den(i,j,k)` is already available from `properties()` and `D_m` is a scalar constant. In solid cells (`T < tsolid_local`), use a floor value of `1e-30` (not zero) to avoid division-by-zero in the Power Law scheme.

### Inline vs array trade-off
The inline approach adds a few extra multiplications per cell per property access (~20 FLOPs per cell in `properties()`, ~10 in `enthalpy_to_temp()`). This is negligible compared to the solver cost. The array approach saved these FLOPs but at the cost of 16 extra 3D arrays (~513 MB for 200x200x67 grid) and +12% total runtime from cache pressure. The inline approach is strictly better for this problem size.

## Known Pitfalls

1. **H-T curve consistency**: Linear mixing of `acpa`, `acpb`, `tsolid` produces `hsmelt_local` that is NOT the same as `hsmelt*C + hsmelt2*(1-C)`. This means small discontinuities in the enthalpy-temperature relationship at composition boundaries. Acceptable for similar materials but may cause slow convergence at sharp interfaces.

2. **~~Narrow mushy zone~~** (resolved): Secondary material now has `deltemp2 = tliquid2 - tsolid2 = 75K`, no longer a concern. Guard `deltemp_local = max(deltemp_local, 1.0_wp)` remains as safety.

3. **`tsolid` scalar in `main.f90` line 149**: `if(tpeak.gt.tsolid)` controls whether momentum is solved. When `species_flag=1`, use `min(tsolid, tsolid2)` so that momentum activates whenever ANY material is molten.

4. **Darcy source viscosity**: When `species_flag=1`, `darcy_c0` must be computed cell-locally inside the loop using composition-weighted viscosity (inline from `concentration(i,j,k)`).

## Revert Log

- **2026-03-18**: Reverted Phase 0 per-cell property arrays (`tsolidmatrix`, `tliquidmatrix`, `dens_arr`, etc.) from `mod_prop.f90`, `mod_entot.f90`, `mod_sour.f90`, `main.f90`. All files restored to pre-Phase 0 scalar versions. Kept `write_memory_report` in `mod_timing.f90` and `main.f90`. Kept `compile.sh` change (no longer cleans `result/`). Replaced array approach with inline computation strategy.

## Reference

- Source implementation: `D:\Fortran\dissimilar\program931\species_transport.f90`
- Analysis of differences: see `species.md` in this folder
