# Species Transport Integration

## Objective

Create `mod_species.f90` in fortran_new based on `D:\Fortran\dissimilar\program931` to solve dissimilar metal species transport. The module must be as self-contained as possible (parameters, initialization, BCs, source terms, solver, output all inside `mod_species.f90`). Minimize changes to existing codebase. Only addition to existing input: `species_flag` (1=on, 0=off).

## Design Principles

- `mod_species.f90` is a standalone module: all secondary material properties, species arrays, solver, BCs, output defined within
- When `species_flag=0`, existing code behavior is unchanged — all property arrays remain constant or T-dependent as before
- When `species_flag=1`, `mod_species.f90` provides functions to update property arrays as composition-dependent
- Species transport active only in melt pool region (same as uvwp): when `T < tsolidmatrix(i,j,k)`, `massdiffusivity` is near-zero
- Concentration C: 1 = base material (from `&material_properties`), 0 = secondary material (from `mod_species.f90`)
- Species solver reuses all shared coefficient arrays (`an, as, ae, aw, at, ab, ap, su, sp, apnot`) — solve after momentum/enthalpy so coefficients are free
- Species arrays are allocated on the global domain, but discretization/source/solver only loop over the melt pool region using `istat,iend,jstat,jend,kstat` from `mod_dimen.f90` (same as uvw momentum equations)
- Species is solved once per timestep after `iter_loop` exits (not inside the iteration loop), using `delt` for transient term

## Key Physics Changes from Original (program931)

- Gamma = rho * D_m (D_m = 5e-9 m^2/s), replacing the old non-physical `massdiffl`. `massdiffusivity` array stores Gamma = rho * D_m (units: kg/(m*s)), used directly as diffusion coefficient in FVM discretization
- Remove velocity x2 factor
- Remove scanning source term (moving frame advection)
- Remove solidification segregation source (kp term)
- All BCs: zero-flux Neumann
- Add solutal Marangoni effect (dgdc based on local composition)
- `dgdt` and `hlatnt` remain scalars (assume both materials have similar values)

## Tasks (easy -> hard)

### Phase 0: Refactor scalars to arrays (no physics change)

1. [ ] **Convert key thermal property scalars to per-cell arrays**
   - Goal: replace scalar properties with arrays in all per-cell computations, initialized to the original scalar values. Results must be bit-identical to the scalar version.
   - Create arrays in `mod_prop.f90`:
     - Phase change: `tsolidmatrix(ni,nj,nk)`, `tliquidmatrix(ni,nj,nk)`
     - Dense solid: `dens_arr(ni,nj,nk)` → `dens`, `acpa_arr`, `acpb_arr`, `thconsa_arr`, `thconsb_arr`
     - Liquid: `denl_arr(ni,nj,nk)` → `denl`, `acpl_arr`, `thconl_arr`
     - Viscosity: `viscmatrix(ni,nj,nk)` → `viscos`
     - Powder: `pden_arr(ni,nj,nk)` → `pden`, `pcpa_arr`, `pcpb_arr`, `pthcona_arr`, `pthconb_arr`
     - All initialized to primary material scalar values
   - Replace scalar usage with array lookups in:
     - `mod_entot.f90` (`enthalpy_to_temp`): compute `hsmelt_local`, `hlcal_local`, `deltemp_local`, `cpavg_local` from arrays as local variables inside the loop. Guard `deltemp_local = max(..., 1.0_wp)`
     - `mod_prop.f90` (`properties`): all three temperature regimes use arrays:
       - `T >= tliquidmatrix`: `denl_arr`, `viscmatrix`, `thconl_arr/acpl_arr`
       - `T <= tsolidmatrix`: `dens_arr`, `acpa_arr/acpb_arr`, `thconsa_arr/thconsb_arr`; powder region uses `pden_arr`, `pcpa_arr/pcpb_arr`, `pthcona_arr/pthconb_arr`
       - Mushy zone: mix of solid/liquid arrays weighted by `fracl`
     - `mod_sour.f90` (`source_momentum`, `source_pp`): use `tsolidmatrix(i,j,k)` for Darcy/buoyancy/solid checks; compute `boufac` cell-locally from `dens_arr * g * beta`; use `viscmatrix(i,j,k)` for Darcy coefficient inside the loop
     - `mod_bound.f90`: no change needed (`dgdt` stays scalar)
     - `main.f90`: velocity zeroing uses `tsolidmatrix(i,j,k)`; melt pool check `if(tpeak > min(tsolid, tsolid2))` (or just `tsolid` for now since arrays = scalar)
   - **Verification**: run existing test case, confirm identical results

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
     hsmelt2 = acp2*tsolid2**2/2 + acpb2*tsolid2
     hlcal2  = hsmelt2 + cpavg2*(tliquid2 - tsolid2)
     hlfriz2 = hlcal2 + hlatnt
     ```
   - `D_m = 5.0e-9_wp` (molecular mass diffusivity, m^2/s)
   - `urfspecies = 0.7_wp` (dedicated relaxation factor)
   - `dgdc_const` (dg/dC, surface tension concentration coefficient)
   - Allocatable arrays: `concentration(:,:,:)`, `conc_old(:,:,:)`, `massdiffusivity(:,:,:)`
   - Note: `tsolidmatrix`, `tliquidmatrix`, `dens_arr`, `denl_arr`, `viscmatrix`, cp/thcon/powder arrays already created in Task 1
   - Empty subroutine stubs: `allocate_species`, `init_species`, `solve_species`, `species_bc`, `update_material_properties`, `update_massdiffusivity`, `write_species_vtk`
   - Compiles but does nothing yet

### Phase 2: Solver core

4. [ ] **Implement `allocate_species` and `init_species`**
   - Allocate species-only arrays: `concentration`, `conc_old`, `massdiffusivity` to (ni, nj, nk)
   - Initial conditions for single-track dissimilar test:
     - Substrate (z < powder layer top): C = 1 (base material)
     - Powder layer, y >= y_mid: C = 1 (base material)
     - Powder layer, y < y_mid: C = 0 (secondary material)
     - No smooth distance — sharp interface
   - Initialize `massdiffusivity`: floor value `1e-30` where `T < tsolid`, else `dens * D_m`
   - Update property arrays from initial C via `update_material_properties` (Task 8)
   - Initialize `conc_old = concentration`

5. [ ] **Implement zero-flux BCs (`species_bc`)**
   - All 6 faces: `C(boundary) = C(interior neighbor)`
   - Called once per timestep, before species discretization

6. [ ] **Implement species solver (`solve_species`)**
   - Port FVM discretization from `program931/species_transport.f90`
   - Loop over melt pool region only: use `istat,iend,jstat,jend,kstat` from `mod_dimen.f90` (same index ranges as uvw)
   - Power Law scheme for convection-diffusion
   - Implicit Euler transient term: `apnot(i,j,k) = densmatrix(i,j,k) / delt * volume(i,j,k)`
   - Diffusion coefficient: use `massdiffusivity(i,j,k)` directly (already stores Gamma = rho * D_m)
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
     - Once per timestep, after `iter_loop` exits: call `species_bc`, `update_massdiffusivity`, then `solve_species`
     - End of timestep: update `conc_old = concentration` — must be full-domain, unconditional of `is_local`
   - Add `resorc` to output line
   - Add species solve time to timing report

### Phase 3: Material property coupling

Since Task 1 already converted scalars to arrays and modified `enthalpy_to_temp`, `properties`, `mod_sour.f90`, `main.f90` to use arrays, these tasks only need to make the arrays composition-dependent.

8. [ ] **Implement `update_material_properties` in `mod_species.f90`**
   - When `species_flag=1`, update all property arrays each timestep based on local C (linear mixing):
     - Phase change: `tsolidmatrix = tsolid*C + tsolid2*(1-C)`, same for `tliquidmatrix`
     - Dense solid: `dens_arr = dens*C + dens2*(1-C)`, `acpa_arr`, `acpb_arr`, `thconsa_arr`, `thconsb_arr`
     - Liquid: `denl_arr = denl*C + denl2*(1-C)`, `acpl_arr`, `thconl_arr`
     - Viscosity: `viscmatrix = viscos*C + viscos2*(1-C)`
     - Powder: `pden_arr = pden*C + pden2*(1-C)`, `pcpa_arr`, `pcpb_arr`, `pthcona_arr`, `pthconb_arr`
     - Note: `beta` remains scalar (similar values for both materials)
     - Note: `hlatnt` and `dgdt` remain scalars (assume similar values for both materials)
   - When `species_flag=0`, arrays stay at scalar initial values (set in Task 1) — no call needed
   - **Call once per timestep, before `iter_loop`** (concentration does not change during iterations)
   - No changes to `properties`, `enthalpy_to_temp`, `mod_sour.f90` needed — they already read from arrays (Task 1)

9. [ ] **Implement `update_massdiffusivity`**
   - `massdiffusivity(i,j,k) = densmatrix(i,j,k) * D_m` when `temp(i,j,k) >= tsolidmatrix(i,j,k)`, else `1e-30` (floor, not zero)
   - Called each timestep before species solve, after `enthalpy_to_temp` updates temperature

### Phase 4: Output and Marangoni

10. [ ] **Implement `write_species_vtk`**
    - Write concentration field as standalone VTK file (same format as defect VTK)
    - ASCII header + binary data, structured grid
    - Called at same frequency as `Cust_Out` (every `outputintervel` steps)
    - Filename: `{case_name}_species{N}.vtk`
    - Also add concentration as a scalar in main VTK output (`Cust_Out`)

11. [ ] **Solutal Marangoni effect (dgdc)**
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

12. [ ] **Create single-track test toolpath**
    - Generate a simple single-track toolpath scanning along x-direction at y = half of domain
    - Use `toolpath_generator_rectangle.py` with 1 track (or write manually)
    - Purpose: test species transport with laser scanning across the material interface

13. [ ] **Validation run**
    - Run single-track test with `species_flag=1`
    - Verify: concentration stays in [0,1], mixing occurs only in melt pool, material properties update correctly, VTK output readable in ParaView
    - Compare melt pool shape with `species_flag=0` to verify minimal impact when C=1 everywhere

### Optional / Future

14. [ ] **`enhance_species_speed` block correction**
    - Port 1D block-correction from program931 for faster convergence
    - Not required for correctness, but improves convergence speed

## Design Notes

### Coefficient arrays
Species solver reuses all shared coefficient arrays (`an, as, ae, aw, at, ab, ap, su, sp, apnot` from `mod_coeff_data.f90`). Species is solved once per timestep after `iter_loop` exits, so the arrays are free.

### Local solver interaction
Species is always solved globally (full domain), not restricted to local solver region. Uses `delt` (not `delt_eff`) for the transient term. `conc_old = concentration` update must be full-domain at end of each timestep, unconditional of `is_local`.

### Diffusivity convention
`massdiffusivity(i,j,k)` stores Gamma = rho * D_m (units: kg/(m*s)). This is the FVM diffusion coefficient used directly in face conductance: `D_face = Gamma * A / dx`. In solid cells, use a floor value of `1e-30` (not zero) to avoid division-by-zero in the Power Law scheme.

### Scalar-to-array strategy (Task 1)
Task 1 converts all key thermal property scalars to per-cell arrays, initialized to the original scalar values. All per-cell computations (`enthalpy_to_temp`, `properties`, `source_momentum`, `source_pp`, velocity zeroing) are modified to read from arrays instead of scalars. With `species_flag=0`, arrays stay at scalar initial values — results are bit-identical to the original code. With `species_flag=1`, `update_material_properties` fills arrays with composition-weighted values — no further code changes needed.

Derived constants (`hsmelt`, `hlcal`, `deltemp`, `cpavg`) are computed as local variables inside the `enthalpy_to_temp` loop from the per-cell arrays. Guard `deltemp_local = max(deltemp_local, 1.0_wp)` as a safety measure.

## Known Pitfalls

1. **H-T curve consistency**: Linear mixing of `acpa`, `acpb`, `tsolid` produces `hsmelt_local` that is NOT the same as `hsmelt*C + hsmelt2*(1-C)`. This means small discontinuities in the enthalpy-temperature relationship at composition boundaries. Acceptable for similar materials but may cause slow convergence at sharp interfaces.

2. **~~Narrow mushy zone~~** (resolved): Secondary material now has `deltemp2 = tliquid2 - tsolid2 = 72K`, no longer a concern. The `max(deltemp_local, 1.0_wp)` guard remains as a safety measure.

3. **`tsolid` scalar remains in `main.f90` line 148**: `if(tpeak.gt.tsolid)` controls whether momentum is solved. When `species_flag=1`, use `min(tsolid, tsolid2)` so that momentum activates whenever ANY material is molten.

4. **Darcy source viscosity**: `mod_sour.f90` line 52 computes `darcy_c0 = 180*viscos/perm_const` using scalar viscosity before the loop. When `species_flag=1`, this should be computed cell-locally inside the loop using composition-weighted viscosity.

## Reference

- Source implementation: `D:\Fortran\dissimilar\program931\species_transport.f90`
- Analysis of differences: see `species.md` in this folder
