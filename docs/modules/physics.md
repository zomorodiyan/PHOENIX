# Physics Modules

Laser, toolpath, melt pool detection, species transport, adaptive mesh.

---

## mod_laser.f90 — `laserinput`

Laser beam management. Each timestep: interpolates toolpath to find beam position, computes scan velocity, distributes power as 2D Gaussian on top surface.

---

## mod_toolpath.f90 — `toolpath`

- `read_toolpath()` — reads `.crs` file into `toolmatrix(1000, 5)`: time, x, y, z, laser_flag
- `read_coordinates()` — records beam state to rolling buffer each timestep

---

## mod_dimen.f90 — `dimensions`

Melt pool detection from temperature field:

- Computes `alen` (length), `depth`, `width` via linear interpolation at solidus isotherm
- Sets momentum solver bounds (`istat:iend`, `jstat:jend`, `kstat`) with padding

---

## mod_adaptive_mesh.f90 — `adaptive_mesh_mod`

Movable adaptive structured mesh in X-Y that follows the laser/melt pool. Enabled by `adaptive_flag=1`.

- `amr_init()` — initializes the refined region after grid generation
- `amr_check_remesh(step_idx)` — checks if remeshing is needed (every `remesh_interval` steps)
- `amr_regenerate_grid()` — rebuilds X-Y grid with refined region centered on laser
- `amr_interpolate_all_fields()` — interpolates all fields from old grid to new grid
- `amr_validate_grid()` — sanity checks on new grid

---

## mod_predict.f90 — `prediction`

Field prediction by integer-cell shifting in the scan direction. Shifts enthalpy, velocity, and pressure fields forward by the number of cells the laser traverses in one timestep. Called before the iteration loop on heating steps when `predict_flag=1` and a melt pool exists. Reduces solver iterations by approximately 38% during heating.

- `predict_shift_integer(vx, vy, dt, ...)` — computes integer cell shift from scan velocity, shifts all primary fields in extended melt pool region
- `ishift_field(field, di, dj, ...)` — generic integer-cell shift with direction-aware iteration order (no interpolation)

---

## mod_species.f90 — `species`

Dissimilar metal species transport. See [Species Transport](../species/overview.md).

Key functions: `mix(prop1, prop2, C)`, `solve_species()` (FVM + TDMA + block correction), `species_bc()` (zero-flux Neumann).
