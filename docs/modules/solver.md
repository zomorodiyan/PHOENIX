# Solver Modules

FVM discretization, SIMPLE algorithm, TDMA solver, convergence acceleration.

---

## mod_prop.f90 — `property`

Updates material properties (viscosity, diffusivity, density) based on temperature and phase state. Three regimes: solid, mushy, liquid. Powder cells use reduced properties. When `species_flag=1`, uses composition-weighted mixing via `mix()`.

---

## mod_bound.f90 — `boundary`

Boundary conditions for all equations:

- **Enthalpy**: radiation + convection (top), convection (sides/bottom)
- **Momentum**: Marangoni stress at free surface (thermal + solutal), no-slip elsewhere
- **Pressure correction**: zero in solid region

---

## mod_discret.f90 — `discretization`

FVM discretization using power-law scheme:

- `discretize_momentum(idir)` — u/v/w momentum with staggered grid interpolation
- `discretize_enthalpy(...)` — energy equation with `delt` for transient term
- `discretize_pp()` — pressure correction equation

---

## mod_sour.f90 — `source`

Source terms:

- **Momentum**: Darcy resistance (mushy zone), buoyancy (Boussinesq), solid zeroing
- **Enthalpy**: volumetric laser source (Gaussian), latent heat, latent heat advection
- **Pressure**: mass conservation assembly

---

## mod_entot.f90 — `entotemp`

Inverts enthalpy → temperature relationship. Handles solid/mushy/liquid regimes. When `species_flag=1`, uses composition-dependent H-T curve via `mix()`.

---

## mod_solve.f90 — `solver`

- `tdma_solve_3d_2(...)` — OpenMP line-by-line TDMA with double k-sweep and per-thread buffers
- `solution_uvw/enthalpy(...)` — convenience wrappers
- `cleanuvw()` — zeros velocity in solid cells and above boiling point

---

## mod_converge.f90 — `convergence`

Block correction: accumulates j-k plane coefficients into 1D system along x, solves with TDMA. Reduces low-frequency errors.

---

## mod_revise.f90 — `revision`

SIMPLE pressure-velocity correction: updates velocity from pressure correction field.

---

## mod_resid.f90 — `residue`

Residual computation for convergence monitoring: enthalpy, mass, momentum (weighted by reference values).

---

## mod_flux.f90 — `fluxes`

Energy balance: boundary losses on all 6 faces, heat accumulation, energy ratio for convergence checking.
