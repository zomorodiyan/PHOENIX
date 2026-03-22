# Module Reference

Detailed documentation for every module in `fortran_new/`.

---

## mod_precision.f90 — `precision`

Defines working precision for the entire codebase.

```fortran
integer, parameter :: wp = selected_real_kind(6, 37)  ! single precision
```

All floating-point variables use `real(wp)`. Change to `selected_real_kind(15, 307)` for double precision.

---

## mod_const.f90 — `constant`

Physical constants and simulation control parameters.

| Constant | Value | Description |
|----------|-------|-------------|
| `g` | 9.8 | Gravitational acceleration (m/s^2) |
| `pi` | 3.14159... | Pi |
| `sigm` | 5.67e-8 | Stefan-Boltzmann constant (W/m^2/K^4) |
| `great` | 1e20 | Large number (for zeroing solid velocity) |
| `small` | 1e-6 | Small number (avoid division by zero) |
| `conv_res_heat` | 1e-5 | Enthalpy convergence threshold (heating) |
| `conv_res_cool` | 1e-6 | Enthalpy convergence threshold (cooling) |
| `vis_solid` | 1e10 | Effective viscosity in solid (Pa*s) |
| `powder_threshold` | 0.5 | solidfield threshold for powder detection |

---

## mod_cfd_utils.f90 — `cfd_utils`

Pure utility functions for CFD calculations.

### `temp_to_enthalpy(T, acpa, acpb, acpl, tsolid, tliquid, hsmelt, hlcal, deltemp)`
Converts temperature to enthalpy using piecewise model:

- Solid ($T \leq T_s$): $H = \frac{a}{2}T^2 + bT$
- Mushy ($T_s < T < T_l$): linear interpolation
- Liquid ($T \geq T_l$): $H = H_l + c_{p,l}(T - T_l)$

### `harmonic_mean(val1, val2, frac)`
Harmonic-mean interpolation for face properties: $\bar{k} = \frac{1}{f/k_1 + (1-f)/k_2}$. Used for viscosity and diffusivity at cell faces.

### `power_law_coeff(diff, flux)`
Power-law discretization scheme: returns $D \cdot \max(0, (1 - 0.1|F/D|)^5) + \max(0, -F)$. Blends central differencing (low Peclet) and upwind (high Peclet).

### `darcy_resistance(viscos, fracl)`
Carman-Kozeny model for mushy zone: $S = \frac{180 \mu}{K_0} \frac{(1-f_l)^2}{f_l + \epsilon}$

---

## mod_param.f90 — `parameters`

Reads all simulation parameters from `input_param.txt`.

### `read_data()`
Parses geometry (zones, CVs, exponents) and 7 namelist groups. Creates result directory `result/<case_name>/` automatically.

See [Input File Reference](../input-reference.md) for parameter details.

---

## mod_geom.f90 — `geometry`

Generates the 3D structured grid.

### `generate_grid()`
Entry point: calls `generate_1d_grid` for x, y, z directions, then computes all geometric quantities.

### `generate_1d_grid(nzones, zones, ncv, powr, vel_grid, scalar_grid, n, nm1)`
Generates 1D non-uniform grid for one coordinate direction using power-law spacing within each zone. Outputs:

- `scalar_grid`: cell-center positions (x, y, z)
- `vel_grid`: face positions (xu, yv, zw)
- `n`, `nm1`: total cells, total cells minus 1

### Geometric quantities computed:
- `volume(i,j,k)`: scalar CV volume
- `volume_u/v/w`: staggered CV volumes
- `areaij/jk/ik`: scalar face areas
- `areauij/uik`, `areavjk/vij`, `areawik/wjk`: staggered face areas
- `dxpwinv/dypsinv/dzpbinv`: inverse distances between scalar nodes (for gradient computation)
- `fracx/fracy/fracz`: interpolation fractions

---

## mod_field_data.f90 — `field_data`

Allocates primary flow field arrays.

### `allocate_field_data(nni, nnj, nnk)`
Allocates all arrays to `(nni, nnj, nnk)`:

| Array | Description |
|-------|-------------|
| `uVel, vVel, wVel` | Current velocity components |
| `unot, vnot, wnot` | Previous time step velocity |
| `pressure, pp` | Pressure and pressure correction |
| `enthalpy, hnot` | Current and previous enthalpy |
| `temp, tnot` | Current and previous temperature |
| `fracl, fraclnot` | Current and previous liquid fraction |
| `solidfield` | Track ID that solidified each cell |
| `localfield` | Local solver region visualization |

---

## mod_coeff_data.f90 — `coeff_data`

Allocates FVM discretization coefficient arrays.

### `allocate_coeff_data(nni, nnj, nnk)`

| Array | Description |
|-------|-------------|
| `vis, diff, den` | Material properties (viscosity, diffusivity, density) |
| `an, as, ae, aw, at, ab` | Neighbor coefficients (north, south, east, west, top, bottom) |
| `ap` | Central coefficient |
| `apnot` | Transient term coefficient |
| `su, sp` | Explicit and implicit source terms |
| `dux, dvy, dwz` | Velocity correction coefficients (for SIMPLE) |

All arrays are **shared** between equations — enthalpy, u, v, w, pp, and species all reuse them sequentially.

---

## mod_sim_state.f90 — `sim_state`

Global simulation state variables.

| Variable | Description |
|----------|-------------|
| `dgdt` | Thermal Marangoni coefficient (= `dgdtp` from input) |
| `deltemp` | Liquidus - solidus temperature difference |
| `cpavg` | Average specific heat in mushy zone |
| `hlcal` | Enthalpy at liquidus temperature |
| `boufac` | Buoyancy factor: `denl * g * beta` |
| `resorm` | Mass conservation residual |
| `refmom` | Reference momentum for residual normalization |
| `beam_pos, beam_posy` | Current laser position |
| `toolmatrix(1000,5)` | Toolpath waypoints |
| `coordhistory(5000,8)` | Beam state history |

---

## mod_init.f90 — `initialization`

Startup initialization wrapper.

### `initialize()`
1. Computes derived constants: `hsmelt`, `hlcal`, `boufac`, `dgdt`, `deltemp`
2. Initializes all fields to preheat conditions (velocity=0, T=tempPreheat)
3. Sets enthalpy boundary values from temperature BCs
4. Initializes `solidfield = 0` (no solidification yet)

---

## mod_laser.f90 — `laserinput`

Laser beam management.

### `laser_beam()`
Called every time step:

1. Interpolates `toolmatrix` to find beam position at current `timet`
2. Advances `PathNum` when passing waypoints
3. Computes scan velocity from position differences
4. Distributes laser power as 2D Gaussian on top surface
5. Identifies initial melt pool search region (imin, imax, jmin, jmax)

---

## mod_toolpath.f90 — `toolpath`

### `read_toolpath()`
Reads `.crs` file into `toolmatrix(1000, 5)`. Format: time, x, y, z, laser_flag.

### `read_coordinates()`
Records current beam state (position, power, velocity) to `coordhistory` rolling buffer each time step.

---

## mod_dimen.f90 — `dimensions`

### `pool_size(ilo, ihi, jlo, jhi, klo, khi)`
Detects melt pool extent from temperature field:

1. Scans domain for cells above `tsolid`
2. Computes pool dimensions: `alen` (x-length), `depth` (z), `width` (y)
3. Sets solver bounds with 2-3 cell padding: `istat`, `iend`, `jstat`, `jend`, `kstat`
4. Pre-computes `istatp1 = istat+1`, `iendm1 = iend-1`

These indices are used by momentum solver, species solver, and velocity cleanup.

---

## mod_local_enthalpy.f90 — `local_enthalpy`

### `get_enthalpy_region(step_idx, is_local, ilo, ihi, jlo, jhi, klo, khi)`
Determines if current step is local or global:

- Global: every `localnum+1` steps, or first step
- Local: otherwise, returns small region around melt pool

### `compute_delt_eff()`
Updates effective time step for all cells: `delt_eff(i,j,k) = delt * (n_skipped(i,j,k) + 1)`.

### `update_skipped(..., is_local)`
Increments skip counter for cells outside current solve region; resets to 0 for solved cells.

---

## mod_prop.f90 — `property`

### `properties(ilo, ihi, jlo, jhi, klo, khi)`
Updates material property arrays based on temperature (and concentration when `species_flag=1`):

**Three regimes:**

| Regime | Condition | Viscosity | Diffusivity | Density |
|--------|-----------|-----------|-------------|---------|
| Liquid | $T \geq T_l$ | `viscos` | `thconl/acpl` | `denl` |
| Solid | $T \leq T_s$ | `vis_solid` (1e10) | `(thconsa*T+thconsb)/(acpa*T+acpb)` | `dens` |
| Mushy | $T_s < T < T_l$ | `viscos` | `fracl*diffl + (1-fracl)*diffs` | `fracl*denl + (1-fracl)*dens` |

Powder cells (top layer, not yet melted) use reduced density and conductivity.

When `species_flag=1`: all properties are computed from `mix(prop_primary, prop_secondary, C)` before applying the three-regime logic.

---

## mod_species.f90 — `species`

Species transport for dissimilar metal mixing. See [Species Transport](../species/overview.md) for details.

### Key functions:

| Function | Description |
|----------|-------------|
| `mix(prop1, prop2, C)` | Pure function: `prop1*C + prop2*(1-C)` |
| `allocate_species()` | Allocates concentration arrays |
| `init_species()` | Sets IC: substrate=1, powder half=0 |
| `species_bc()` | Zero-flux Neumann on all faces |
| `solve_species()` | Full solve: discretize + source + TDMA + block correction + clip + residual |
| `solution_species_tdma()` | OpenMP TDMA with per-thread buffers |
| `enhance_species_speed()` | Block correction in x-direction |
| `calc_species_residual()` | Weighted residual for convergence monitoring |

---

## mod_bound.f90 — `boundary`

### `bound_uv(idir)`
Top surface Marangoni boundary condition (k=nk):

$$u_{surface} = u_{below} + \frac{f_l}{vis} \left(\frac{d\gamma}{dT}\frac{\partial T}{\partial x} + \frac{d\gamma}{dC}\frac{\partial C}{\partial x}\right) \frac{1}{\Delta z^{-1}}$$

The solutal Marangoni term (`dgdc * dC/dx`) is only active when `species_flag=1`.

### `bound_w()`
No-slip condition at melt pool boundaries for w-velocity.

### `bound_pp()`
Zero pressure correction in solid region.

### `bound_enthalpy(ilo, ihi, jlo, jhi, klo, khi, is_local)`
- **Top surface**: radiation + convection loss
- **Other faces**: convection to ambient
- **Local mode**: Dirichlet (frozen values) outside local region

---

## mod_discret.f90 — `discretization`

### `discretize_momentum(idir)`
FVM discretization for u (idir=1), v (idir=2), w (idir=3):

1. Interpolates face velocities from staggered grid
2. Computes convection coefficients: $F = \rho u A$
3. Computes face viscosity via harmonic mean
4. Computes diffusion coefficients: $D = \mu A / \Delta x$
5. Applies power-law scheme for neighbor coefficients
6. Adds transient term: $a_{P}^0 = \rho V / \Delta t$
7. Adds pressure gradient and cross-derivative source terms

### `discretize_enthalpy(ilo, ihi, jlo, jhi, klo, khi)`
Same structure as momentum but uses thermal diffusivity instead of viscosity. Uses `delt_eff` for transient term (local solver compensation).

### `discretize_pp()`
Pressure-correction equation using velocity correction coefficients from momentum solve.

---

## mod_entot.f90 — `entotemp`

### `enthalpy_to_temp(ilo, ihi, jlo, jhi, klo, khi)`
Inverts enthalpy-temperature relationship:

| Region | Condition | Temperature | Liquid Fraction |
|--------|-----------|-------------|-----------------|
| Liquid | $H \geq H_l$ | $T = (H - H_l)/c_{p,l} + T_l$ | $f_l = 1$ |
| Solid | $H \leq H_s$ | $T = (-b + \sqrt{b^2 + 2aH})/a$ | $f_l = 0$ |
| Mushy | $H_s < H < H_l$ | $T = \Delta T \cdot f_l + T_s$ | $f_l = (H - H_s)/(H_l - H_s)$ |

When `species_flag=1`: $H_s$, $H_l$, $\Delta T$, $a$, $b$, $c_{p,l}$ are all computed inline from `mix()` using local concentration.

---

## mod_sour.f90 — `source`

### `source_momentum(idir)`
1. **Darcy resistance** in mushy zone: $S_p = -\frac{180\mu}{K_0} \frac{(1-f_l)^2}{f_l + \epsilon} V$
2. **Buoyancy** (w-momentum only): $S_u = \rho g \beta (T - T_s) V$
3. **Assembly**: $a_P = \sum a_{nb} + a_P^0 - S_P$
4. **Under-relaxation**: $a_P = a_P / \alpha$, $S_U += (1-\alpha) a_P \phi$
5. **Solid zeroing**: If $T \leq T_s$, zero all coefficients, set $a_P = 10^{20}$

When `species_flag=1`: Darcy uses composition-weighted viscosity, buoyancy uses composition-weighted density and tsolid.

### `source_pp()`
Assembly of pressure-correction coefficients. Zeros solid cells.

### `source_enthalpy(ilo, ihi, jlo, jhi, klo, khi)`
1. **Laser volumetric source**: Gaussian distribution in x-y, uniform in z within penetration depth
2. **Latent heat**: $S = -\rho h_{lat} (f_l - f_l^{old}) V / \Delta t$
3. **Latent heat advection**: Upwind flux of liquid fraction through cell faces
4. **Boundary transfers**: Converts boundary conditions to source/sink terms
5. **Assembly and under-relaxation**

---

## mod_resid.f90 — `residue`

### `calc_momentum_residual(vel, resor_out, calc_refmom)`
Weighted residual: $r = \sum |a_{nb}\phi_{nb} + S_U - a_P\phi_P| / \text{refmom}$

Reference momentum: $\text{refmom} = 0.25\pi (\min(L,D,W))^2 \rho u_{max}^2$

### `calc_enthalpy_residual(ilo, ihi, jlo, jhi, klo, khi)`
$r_h = \sum|residual| / \sum|H|$

### `calc_pressure_residual()`
$r_m = \sum|\nabla \cdot \vec{u}| / (\rho \sum|\vec{u}| A)$

---

## mod_solve.f90 — `solver`

### `tdma_solve_3d_2(field, klo, khi, jlo, jhi, ibc, ilo, ihi)`
OpenMP-parallelized line-by-line TDMA:

- Double k-sweep (backward then forward) for convergence
- Per-thread pr/qr buffers to avoid data races
- Forward elimination in i-direction, back substitution

### `solution_uvw(field)` / `solution_enthalpy(ilo,...)`
Convenience wrappers that call TDMA with appropriate bounds.

### `cleanuvw()`
Zeroes velocity at staggered nodes where both adjacent cells are solid ($T \leq T_s$). Also zeroes velocity above boiling point (free surface approximation).

---

## mod_converge.f90 — `convergence`

### `enhance_converge_speed(ilo, ihi, jlo, jhi, klo, khi)`
Block correction: accumulates j-k plane coefficients into a 1D system along x, solves with TDMA, broadcasts correction. Reduces low-frequency errors that the line-by-line TDMA handles slowly.

---

## mod_revise.f90 — `revision`

### `revision_p()`
SIMPLE pressure-velocity correction: $u' = u + d_u(p'_W - p'_P)$ where $d_u = A_{face} / a_P$.

---

## mod_flux.f90 — `fluxes`

### `heat_fluxes()`
Computes energy balance: boundary conduction losses on all 6 faces, heat accumulation (enthalpy change + latent heat), and energy ratio $R = Q_{in} / Q_{out}$. Used for convergence checking ($R \approx 1.0$).

---

## mod_defect.f90 — `defect_field`

Post-simulation defect detection based on maximum temperature history. Identifies lack-of-fusion (incomplete melting) and keyhole porosity (excessive vaporization) within the powder layer.

### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `k_lof` | 0.9 | Lack-of-fusion calibration factor (0–1) |
| `k_kep` | 1.0 | Keyhole porosity scaling factor |

### Arrays

| Array | Shape | Description |
|-------|-------|-------------|
| `max_temp(ni,nj,nk)` | Full domain | Peak temperature reached during simulation (only `k_def_lo:k_def_hi` active) |
| `defect_arr(ni,nj,nk)` | Full domain | Defect classification field |

### `allocate_defect(nni, nnj, nnk)`
Allocates arrays and determines the Z-range for defect analysis. The active range is from the top of the domain (`nkm1`) down to `z(nk) - layerheight`, covering only the powder/build layer.

### `update_max_temp()`
Called every time step. For each cell in the build layer (`k_def_lo:k_def_hi`), updates `max_temp(i,j,k)` if the current temperature exceeds the stored maximum. OpenMP parallelized.

### `compute_defect_determ()`
Called once after the simulation ends. Three-step process:

1. **Classify each cell** based on peak temperature:

    | Condition | Classification | `defect_arr` value |
    |-----------|---------------|-------------------|
    | $T_{max} < T_s$ | Lack of fusion | $-k_{lof}$ (negative) |
    | $T_s \leq T_{max} \leq T_b$ | Sound (defect-free) | $0$ |
    | $T_{max} > T_b$ | Keyhole porosity | $\min\left(k_{kep} \frac{T_{max} - T_b}{T_b},\ 0.99\right)$ |

2. **Build scan region polygon** from toolpath via `compute_scan_range()`
3. **Zero defects outside scanned region** — cells that were never under the laser path are not defects

### `compute_scan_range()`
Builds a convex polygon (counter-clockwise) representing the laser-scanned region from the toolpath data:

- Extracts left and right endpoints of each laser-on track
- Constructs polygon: right boundary (bottom→top) + left boundary (top→bottom)
- Also computes axis-aligned bounding box (`x_scan_min/max`, `y_scan_min/max`)

### `point_in_scan_region(px, py)`
Tests if point `(px, py)` is inside the convex hull polygon using the cross-product winding method. Returns `.true.` if the point is inside (all cross products ≥ 0 for CCW polygon). Used to exclude cells outside the scanned area from defect statistics.

### `write_defect_report()`
Computes volumetric defect fractions and generates output:

1. **Metrics** (over scanned region only):
    - Total reference volume $V_{total}$
    - Lack-of-fusion volume: $V_{lof} = \sum |d_i| \cdot V_i$ where $d_i < 0$
    - Keyhole porosity volume: $V_{kep} = \sum d_i \cdot V_i$ where $0 < d_i \leq 1$
    - Defect fraction: $(V_{lof} + V_{kep}) / V_{total} \times 100\%$

2. **Output files**:
    - `<case>_defect_report.txt` — text summary with fractions, volumes, scan region, parameters
    - Summary printed to `output.txt`
    - Calls `write_defect_vtk` for VTK output

### `write_defect_vtk(fieldname, field)`
Writes a binary VTK structured-grid file for the build layer Z-range only (`k_def_lo:k_def_hi`). Used to output both `maxtemp.vtk` and `defect.vtk`. Same format as `Cust_Out` VTK files (ASCII header + big-endian float32 binary data).

### `maxtemp_stochas()`
Placeholder for future stochastic defect prediction method. Currently a no-op.

---

## mod_timing.f90 — `timing`

### `write_timing_report(itertot, timet_end, wall_elapsed, file_prefix)`
Writes detailed breakdown of CPU time by module (17 categories), sorted by percentage. Includes subroutine-level breakdown for modules >10%.

### `write_memory_report(file_prefix)`
Reads `/proc/self/status` for VmPeak, VmHWM (peak physical RAM), VmRSS, VmData.

---

## mod_print.f90 — `printing`

### `outputres()`
Writes one block per time step to `output.txt`: residuals, velocities, pool dimensions, energy balance, progress.

### `Cust_Out()`
Writes binary VTK file every `outputintervel` steps. Fields: Velocity (vector), T, vis, diff, den, solidID, localfield. When `species_flag=1`: adds concentration and tsolid_field.

### `init_thermal_history()` / `write_thermal_history(timet)` / `finalize_thermal_history()`
Tracks temperature at 10 monitoring points throughout simulation. Generates Python plotting script and PNG at end.

### `init_meltpool_history()` / `write_meltpool_history(timet)` / `finalize_meltpool_history()`
Tracks melt pool length, depth, width, volume, and peak temperature. Only recorded during global solver timesteps. Generates Python plot script and PNG at end.

### `write_vtk_scalar(unit, filename, name, field)` / `write_vtk_vector(...)`
Binary VTK writers: ASCII header (SCALARS/VECTORS + LOOKUP_TABLE) followed by big-endian float32 data.

---

## mod_microstructure.f90 — `microstructure_mod`

Solidification microstructure prediction. Enabled by `micro_flag=1`. Only updates during global solver timesteps.

### `allocate_microstructure(nni, nnj, nnk)`
Allocates arrays: `cool_rate_micro`, `therm_grad`, `solid_rate`, `pdas_arr`, `sdas_arr`, `micro_solidified`.

### `update_microstructure(dt)`
Called each global timestep. Detects solidification events (`fraclnot > 0` → `fracl ≤ 0`). At solidification: computes cooling rate `|dT/dt|`, thermal gradient `G = |∇T|` (central differences), solidification rate `R = |dT/dt| / G`, PDAS `λ₁ = a₁ G^n₁ R^n₂`, SDAS `λ₂ = a₂ |dT/dt|^n₃`. Only records the first solidification event per cell.

### `report_microstructure()`
Post-simulation. Computes min/max/mean statistics over solidified cells in the scan region. Writes `micro_report.txt` and a single `microstructure.vtk` with 5 scalar fields.

---

## mod_crack_risk.f90 — `crack_risk_mod`

Crack risk prediction from thermal strain in the Brittle Temperature Range (BTR). Enabled by `crack_flag=1`. Only updates during global solver timesteps.

### `allocate_crack_risk(nni, nnj, nnk)`
Allocates arrays: `cool_rate_solid`, `strain_rate_solid`, `btr_time`, `crack_risk_arr`, `crack_solidified`.

### `update_crack_risk(dt)`
Called each global timestep. Detects solidification (stores cooling rate and strain rate = β × |dT/dt|). Accumulates BTR strain: when `T_solidus - ΔT_BTR < T < T_solidus`, adds `β × |dT/dt| × dt` to the crack susceptibility index (CSI).

### `compute_crack_report()`
Post-simulation. Computes statistics, identifies high-risk cells (CSI > 1%). Writes `crack_report.txt` and a single `crack_risk.vtk` with 4 scalar fields.
