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
- Concentration C: 1 = base material (from `&material_properties`), 0 = secondary material (from `mod_species.f90`)
- Species solver reuses all shared coefficient arrays (`an, as, ae, aw, at, ab, ap, su, sp, apnot`) — solve after momentum/enthalpy so coefficients are free
- Species arrays are allocated on the global domain, but discretization/source/solver only loop over the melt pool region using `istat,iend,jstat,jend,kstat` from `mod_dimen.f90` (same as uvw momentum equations)
- Species is solved once per timestep after `iter_loop` exits (not inside the iteration loop), using `delt` for transient term

## Inline Computation Strategy (replaces Phase 0 array approach)

**Rationale**: The Phase 0 approach (16 per-cell property arrays, ~513 MB overhead, +12% runtime) was implemented and reverted. Instead, all composition-dependent properties are computed inline inside hot loops using local variables:

```fortran
! Example: in properties() when species_flag=1
C_local = concentration(i,j,k)
tsolid_local  = tsolid  * C_local + tsolid2  * (1.0_wp - C_local)
tliquid_local = tliquid * C_local + tliquid2 * (1.0_wp - C_local)
dens_local    = dens    * C_local + dens2    * (1.0_wp - C_local)
! ... etc. for all properties used in this routine
```

When `species_flag=0`, the original scalar code paths remain untouched — no `if` checks, no concentration lookups, no overhead.

**Benefits over array approach**:
- Zero memory overhead (no 16 extra 3D arrays)
- No `allocate_prop_arrays` call needed
- No compile-order dependency changes
- Properties stay "live" — always consistent with current concentration (no stale array risk)
- Code changes confined to the `species_flag=1` branch — existing single-material behavior is untouched

**Implementation pattern for modified routines** (properties, enthalpy_to_temp, source_momentum, source_pp, velocity zeroing):
```fortran
if (species_flag == 1) then
    ! Composition-dependent path: compute local properties from C
    C_local = concentration(i,j,k)
    tsolid_local = tsolid * C_local + tsolid2 * (1.0_wp - C_local)
    ! ... use tsolid_local in branching and computation
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

2. [ ] **Add `species_flag` to input parsing**
   - Add `species_flag` to `&output_control` namelist in `mod_init.f90`
   - Declare `species_flag` in `mod_sim_state.f90` (default = 0)
   - No behavior change yet

3. [ ] **Create `mod_species.f90` skeleton**
   - Module with all secondary material parameters as named constants:
     ```
     ! Dense solid/liquid properties
     dens2=8880, denl2=7800, viscos2=0.003,
     tsolid2=1728, tliquid2=1803, tboiling2=3650,
     acpa2=0.3441, acpb2=400, acpl2=800,
     thconsa2=0.0205, thconsb2=10, thconl2=120
     ! Powder properties
     pden2=7330, pcpa2=0.3508, pcpb2=457.7, pthcona2=0, pthconb2=0.795
     ```
   - Derived quantities computed in `init_species` (not input parameters):
     ```
     hsmelt2 = acpa2*tsolid2**2/2 + acpb2*tsolid2
     hlcal2  = hsmelt2 + cpavg2*(tliquid2 - tsolid2)
     hlfriz2 = hlcal2 + hlatnt
     ```
   - `D_m = 5.0e-9_wp` (molecular mass diffusivity, m^2/s)
   - `urfspecies = 0.7_wp` (dedicated relaxation factor)
   - `dgdc_const` (dg/dC, surface tension concentration coefficient)
   - Allocatable arrays: `concentration(:,:,:)`, `conc_old(:,:,:)`
   - No `massdiffusivity` array — use `den(i,j,k) * D_m` inline (den already computed in `properties()`, D_m is scalar constant)
   - Empty subroutine stubs: `allocate_species`, `init_species`, `solve_species`, `species_bc`, `write_species_vtk`
   - Compiles but does nothing yet

### Phase 2: Solver core

4. [ ] **Implement `allocate_species` and `init_species`**
   - Allocate species-only arrays: `concentration`, `conc_old` to (ni, nj, nk)
   - Initial conditions for single-track dissimilar test:
     - Substrate (z < powder layer top): C = 1 (base material)
     - Powder layer, y >= y_mid: C = 1 (base material)
     - Powder layer, y < y_mid: C = 0 (secondary material)
     - No smooth distance — sharp interface
   - Initialize `conc_old = concentration`

5. [ ] **Implement zero-flux BCs (`species_bc`)**
   - All 6 faces: `C(boundary) = C(interior neighbor)`
   - Called once per timestep, before species discretization

6. [ ] **Implement species solver (`solve_species`)**
   - Port FVM discretization from `program931/species_transport.f90`
   - Loop over melt pool region only: use `istat,iend,jstat,jend,kstat` from `mod_dimen.f90` (same index ranges as uvw)
   - Power Law scheme for convection-diffusion
   - Implicit Euler transient term: `apnot(i,j,k) = den(i,j,k) / delt * volume(i,j,k)` (use `den` from properties, already composition-weighted when species_flag=1)
   - Diffusion coefficient: `Gamma = den(i,j,k) * D_m` computed inline (no array needed). In solid cells (`T < tsolid_local`), use floor `1e-30` instead
   - **Fix: remove velocity x2 factor** — use raw velocities directly
   - **Remove: scanning source term** (no moving frame advection)
   - **Remove: solidification segregation** (no kp term)
   - Reuses all shared coefficient arrays (`an, as, ae, aw, at, ab, ap, su, sp, apnot`)
   - Line-by-line TDMA solver (same pattern as enthalpy)
   - Use `urfspecies` for under-relaxation (NOT `urfh`)
   - Concentration clipping [0, 1] after each TDMA sweep
   - Compute and return species residual (`resorc`)

7. [ ] **Wire into main loop**
   - In `main.f90`:
     - Startup: call `allocate_species` and `init_species` (when `species_flag=1`)
     - Once per timestep, after `iter_loop` exits: call `species_bc`, then `solve_species`
     - End of timestep: update `conc_old = concentration` — must be full-domain, unconditional of `is_local`
   - Add `resorc` to output line
   - Add species solve time to timing report

### Phase 3: Material property coupling (inline approach)

8. [ ] **Modify `properties()` in `mod_prop.f90` for composition-dependent properties**
   - When `species_flag=1`, compute all properties inline from `concentration(i,j,k)`:
     ```fortran
     use species, only: concentration, species_flag, tsolid2, tliquid2, ...
     ! Inside loop:
     if (species_flag == 1) then
         C_local = concentration(i,j,k)
         tsolid_local  = tsolid  * C_local + tsolid2  * (1.0_wp - C_local)
         tliquid_local = tliquid * C_local + tliquid2 * (1.0_wp - C_local)
         dens_local    = dens    * C_local + dens2    * (1.0_wp - C_local)
         denl_local    = denl    * C_local + denl2    * (1.0_wp - C_local)
         viscos_local  = viscos  * C_local + viscos2  * (1.0_wp - C_local)
         ! ... all cp, thcon, powder properties similarly
     else
         tsolid_local = tsolid; tliquid_local = tliquid; ...
     endif
     ```
   - Use local variables for branching (`tsolid_local`, `tliquid_local`) and property computation
   - Three temperature regimes use locally-mixed properties
   - `beta`, `emiss`, `hlatnt`, `dgdt` remain scalars

9. [ ] **Modify `enthalpy_to_temp()` in `mod_entot.f90` for composition-dependent H-T curve**
   - When `species_flag=1`, compute `hsmelt_local`, `hlcal_local`, `deltemp_local` from mixed `acpa`, `acpb`, `tsolid`, `tliquid`, `acpl`
   - When `species_flag=0`, use original scalars `hsmelt`, `hlcal`, `deltemp`, `acpl`, `tsolid`, `tliquid` (unchanged)

10. [ ] **Modify `source_momentum`/`source_pp` in `mod_sour.f90` for composition-dependent branching**
    - When `species_flag=1`:
      - Darcy: `darcy_c0` computed cell-locally using mixed viscosity
      - Buoyancy: uses mixed `dens_local * g * beta * (tw - tsolid_local)`
      - Solid check: uses mixed `tsolid_local`
    - When `species_flag=0`: original scalar code (unchanged)

11. [ ] **Modify velocity zeroing in `main.f90`**
    - When `species_flag=1`: compute `tsolid_local` from `concentration(i,j,k)` for the `temp <= tsolid` check
    - When `species_flag=0`: use scalar `tsolid` (unchanged, current code)
    - Also: `if(tpeak > min(tsolid, tsolid2))` for momentum activation check

### Phase 4: Output and Marangoni

12. [ ] **Implement `write_species_vtk`**
    - Write concentration field as standalone VTK file (same format as defect VTK)
    - ASCII header + binary data, structured grid
    - Called at same frequency as `Cust_Out` (every `outputintervel` steps)
    - Filename: `{case_name}_species{N}.vtk`
    - Also add concentration as a scalar in main VTK output (`Cust_Out`)

13. [ ] **Solutal Marangoni effect (dgdc)**
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

### Phase 5: Testing

14. [ ] **Create single-track test toolpath**
    - Generate a simple single-track toolpath scanning along x-direction at y = half of domain
    - Use `toolpath_generator_rectangle.py` with 1 track (or write manually)
    - Purpose: test species transport with laser scanning across the material interface

15. [ ] **Validation run**
    - Run single-track test with `species_flag=1`
    - Verify: concentration stays in [0,1], mixing occurs only in melt pool, material properties update correctly, VTK output readable in ParaView
    - Compare melt pool shape with `species_flag=0` to verify minimal impact when C=1 everywhere

### Optional / Future

16. [ ] **`enhance_species_speed` block correction**
    - Port 1D block-correction from program931 for faster convergence
    - Not required for correctness, but improves convergence speed

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
