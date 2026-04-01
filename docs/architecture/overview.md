# Architecture Overview

## Governing Equations

PHOENIX solves the coupled thermo-fluid equations for an incompressible Newtonian fluid with phase change in a 3D Cartesian domain.

### Continuity (Mass Conservation)

$$\nabla \cdot (\rho \vec{u}) = 0$$

Density $\rho$ varies spatially (solid, liquid, mushy, powder have different densities), so this is not equivalent to $\nabla \cdot \vec{u} = 0$. The pressure correction equation enforces zero net mass flux through each control volume.

### Momentum (Navier-Stokes)

$$\rho \frac{\partial \vec{u}}{\partial t} + \rho (\vec{u} \cdot \nabla)\vec{u} = -\nabla p + \nabla \cdot (\mu \nabla \vec{u}) + \vec{S}_u$$

Source terms $\vec{S}_u$:

- **Darcy resistance** (mushy zone): $S_D = -\frac{180\mu}{K_0} \frac{(1-f_l)^2}{f_l + \epsilon} \vec{u}$, where $K_0 = 10^{-10}$ m$^2$, $\epsilon = 10^{-3}$
- **Buoyancy** (w-direction only, Boussinesq): $S_b = \rho g \beta (T - T_s) \hat{k}$
- **Solid penalty**: When $T \leq T_s$, all coefficients are zeroed and $a_P = 10^{20}$ to enforce $\vec{u} = 0$

### Energy (Enthalpy Formulation)

$$\rho \frac{\partial H}{\partial t} + \rho (\vec{u} \cdot \nabla) H = \nabla \cdot \left(\frac{k}{c_p} \nabla H\right) + S_H$$

Source terms $S_H$:

- **Laser volumetric heating**: $q_{vol} = \frac{P \eta f}{\pi r_s^2 d_s} \exp\left(-\frac{f}{r_s^2}\left[(x-x_b)^2 + (y-y_b)^2\right]\right)$ for $z \geq z_{top} - d_s$
- **Latent heat**: $S_L = -\rho h_{lat} \frac{\partial f_l}{\partial t}$ (enthalpy-porosity method)

The enthalpy-temperature relationship is piecewise:

| Region | Condition | $H(T)$ | $T(H)$ |
|--------|-----------|---------|---------|
| Solid | $T \leq T_s$ | $\frac{a}{2}T^2 + bT$ | $\frac{-b + \sqrt{b^2 + 2aH}}{a}$ |
| Mushy | $T_s < T < T_l$ | $H_s + \bar{c}_p(T - T_s)$ | $T_s + \Delta T \cdot f_l$ |
| Liquid | $T \geq T_l$ | $H_l + c_{p,l}(T - T_l)$ | $T_l + (H - H_l)/c_{p,l}$ |

where $a$ = `acpa`, $b$ = `acpb`, $\bar{c}_p$ = `cpavg`, $f_l = (H - H_s)/(H_l - H_s)$.

### Species Transport

$$\rho \frac{\partial C}{\partial t} + \rho (\vec{u} \cdot \nabla) C = \nabla \cdot (\rho D_m \nabla C)$$

where $C$ is the mass fraction of primary material, $D_m = 5 \times 10^{-9}$ m$^2$/s. Solved once per timestep after the iteration loop. See [Species Transport](../species/overview.md) for details.

### Material Properties

Temperature-dependent (and composition-dependent when `species_flag=1`):

| Property | Liquid | Solid | Powder |
|----------|--------|-------|--------|
| Density $\rho$ | `denl` | `dens` | `pden` |
| Viscosity $\mu$ | `viscos` | $10^{10}$ (rigid) | $10^{10}$ |
| Diffusivity $\alpha$ | $k_l/c_{p,l}$ | $(a_k T + b_k)/(a_c T + b_c)$ | $(a_{pk} T + b_{pk})/(a_{pc} T + b_{pc})$ |

In the mushy zone ($T_s < T < T_l$), properties are linearly interpolated by liquid fraction $f_l$.

## Boundary Conditions

### Top Surface ($z = z_{max}$, free surface)

- **Thermal**: Combined radiation + convection loss

$$q_{top} = h_{top}(T - T_{amb}) + \epsilon \sigma (T^4 - T_{amb}^4)$$

- **Velocity**: Marangoni stress (thermal + solutal)

$$u_{surface} = u_{below} + \frac{f_l}{\mu / \Delta z} \left(\frac{d\gamma}{dT}\frac{\partial T}{\partial x} + \frac{d\gamma}{dC}\frac{\partial C}{\partial x}\right)$$

$$v_{surface} = v_{below} + \frac{f_l}{\mu / \Delta z} \left(\frac{d\gamma}{dT}\frac{\partial T}{\partial y} + \frac{d\gamma}{dC}\frac{\partial C}{\partial y}\right)$$

The solutal Marangoni term ($d\gamma/dC$) is only active when `species_flag=1`.

- **w-velocity**: $w_{surface} = 0$ (flat free surface approximation)
- **Pressure**: Zero-gradient

### Bottom Surface ($z = 0$)

- **Thermal**: Convection to ambient: $q_{bottom} = h_{bottom}(T - T_{bottom})$
- **Velocity**: No-slip ($\vec{u} = 0$, solid substrate)

### Side Walls ($x = 0$, $x = L_x$, $y = 0$, $y = L_y$)

- **Thermal**: Convection to ambient: $q_{side} = h_{side}(T - T_{side})$
- **Velocity**: No-slip ($\vec{u} = 0$)

### Species (all boundaries)

- **Zero-flux Neumann**: $\frac{\partial C}{\partial n} = 0$ on all 6 faces

## Program Flow

```
main.f90
│
├── Initialization
│   ├── read_data()           ← Parse input_param.txt
│   ├── read_toolpath()       ← Load .crs toolpath
│   ├── generate_grid()       ← Build non-uniform 3D mesh
│   ├── allocate_fields()     ← Allocate all field arrays
│   ├── initialize()          ← Set initial conditions
│   ├── allocate_species()    ← [if species_flag=1]
│   └── amr_init()            ← [if adaptive_flag=1]
│
├── Time Stepping Loop (timet < timax)
│   │
│   ├── laser_beam()          ← Update beam position
│   ├── read_coordinates()    ← Record beam state
│   ├── [if adaptive_flag=1]  ← AMR check + regenerate grid
│   │
│   ├── [if predict_flag=1 AND laser on AND tpeak > tsolid]
│   │   ├── predict_shift_integer() ← Shift fields by integer cells in scan direction
│   │   └── enthalpy_to_temp()      ← Recompute T, fracl from shifted enthalpy
│   │
│   ├── Iteration Loop (niter < maxit)
│   │   │
│   │   ├── properties()       ← Update vis, diff, den from T (and C)
│   │   ├── bound_enthalpy()   ← Enthalpy BCs (radiation, convection)
│   │   ├── discretize_enthalpy() ← FVM coefficients for energy eq
│   │   ├── source_enthalpy()  ← Laser + latent heat source terms
│   │   ├── calc_enthalpy_residual()
│   │   ├── enhance_converge_speed() ← Block correction
│   │   ├── solution_enthalpy() ← TDMA solve for enthalpy
│   │   ├── enthalpy_to_temp() ← H → T, fracl conversion
│   │   ├── pool_size()        ← Find melt pool bounds
│   │   │
│   │   └── [if melt pool exists]
│   │       ├── cleanuvw()     ← Zero velocity in solid
│   │       ├── u-momentum: bound → discretize → source → residual → TDMA
│   │       ├── v-momentum: bound → discretize → source → residual → TDMA
│   │       ├── w-momentum: bound → discretize → source → residual → TDMA
│   │       └── pressure:   bound → discretize → source → residual → TDMA → revision
│   │
│   ├── [after iter_loop]
│   │   ├── species_bc()       ← [if species_flag=1]
│   │   └── solve_species()    ← [if species_flag=1] FVM + TDMA for concentration
│   │
│   ├── update_max_temp()      ← Defect analysis accumulation
│   ├── outputres()            ← Print residuals to output.txt
│   ├── Cust_Out()             ← Write VTK (every outputintervel steps)
│   └── conc_old = concentration ← [if species_flag=1]
│
└── Post-simulation
    ├── compute_defect_determ()  ← Defect classification
    ├── write_defect_report()    ← Defect VTK + report
    ├── write_timing_report()    ← Performance breakdown
    └── write_memory_report()    ← Memory usage
```

## Numerical Method

| Aspect | Method |
|--------|--------|
| Spatial discretization | Finite Volume Method (FVM) on structured grid |
| Grid type | Staggered (velocities at faces, scalars at centers) |
| Convection scheme | Power Law (blends central and upwind) |
| Pressure-velocity coupling | SIMPLE algorithm |
| Linear solver | Line-by-line TDMA with block correction |
| Time integration | Implicit Euler (first-order) |
| Phase change | Enthalpy method with Darcy resistance in mushy zone |
| Parallelization | OpenMP (shared memory, per-thread TDMA buffers) |

## Module Dependency Graph

Compilation order reflects dependencies (each module depends only on modules compiled before it):

```
mod_precision       ← Foundation: working precision (single/double)
  └── mod_const     ← Physical constants, convergence thresholds
      └── mod_cfd_utils  ← Utility functions (harmonic mean, power law, etc.)
      └── mod_param      ← Input parsing, material properties
          └── mod_geom        ← Grid generation, geometric quantities
          └── mod_field_data  ← Velocity, pressure, enthalpy, temperature arrays
          └── mod_coeff_data  ← FVM coefficient arrays (an,as,ae,aw,at,ab,ap,su,sp)
          └── mod_sim_state   ← Global state (residuals, beam position, toolpath)
              └── mod_init         ← Initialization wrapper
              └── mod_laser        ← Laser beam positioning + heat distribution
              └── mod_dimen        ← Melt pool size detection
              └── mod_adaptive_mesh  ← Adaptive mesh (AMR) for laser tracking
              └── mod_resid        ← Residual calculations
              └── mod_species      ← Species transport (dissimilar metals)
              └── mod_prop         ← Temperature/composition-dependent properties
              └── mod_bound        ← Boundary conditions (Marangoni, radiation)
              └── mod_discret      ← FVM discretization (momentum, enthalpy, pp)
              └── mod_entot        ← Enthalpy ↔ temperature conversion
              └── mod_predict      ← Field prediction (integer-cell shift)
              └── mod_sour         ← Source terms (laser, latent heat, Darcy, buoyancy)
              └── mod_flux         ← Energy balance verification
              └── mod_revise       ← Pressure-velocity correction
              └── mod_solve        ← TDMA solver + velocity cleanup
              └── mod_print        ← Output (text, VTK, thermal history)
              └── mod_converge     ← Block correction acceleration
              └── mod_toolpath     ← Toolpath file reading
              └── mod_timing       ← Performance reporting
              └── mod_defect       ← Defect detection and output
```

## Adaptive Mesh

When `adaptive_flag=1`, the solver uses a movable adaptive structured mesh that follows the laser/melt pool in X-Y. The refined region (cell size `amr_dx_fine`) is centered on the current laser position, with coarser cells outside.

Key behavior:

1. **Initialization**: `amr_init()` sets up the initial refined region after grid generation.
2. **Remesh check**: Every `remesh_interval` timesteps, `amr_check_remesh()` determines if the laser has moved far enough to warrant regridding.
3. **Grid regeneration**: `amr_regenerate_grid()` rebuilds the X-Y grid with the refined region centered on the current laser position.
4. **Field interpolation**: `amr_interpolate_all_fields()` maps all field data from the old grid to the new grid.
5. **Validation**: `amr_validate_grid()` checks grid integrity (monotonicity, positive volumes, correct domain extent).

This allows using fine resolution only where needed (near the melt pool) while keeping the total cell count manageable.
