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

### Mechanical

- **Bottom face** ($k=1$): $\mathbf{u} = \mathbf{0}$ (clamped substrate base)
- **All other faces**: Traction-free (natural BC)

## Program Flow

See [Solver Algorithm — Main Program Call Flow](solver.md) for the detailed call graph with module annotations.

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
              └── mod_mech_material  ← Mechanical material properties, J2 return map
                  └── mod_mechanical   ← EBE FEM solver (Newton + CG)
                  └── mod_mech_io      ← Mechanical VTK, timing, history, GIF
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

## Mechanical Solver

When `mechanical_flag=1`, an EBE (Element-By-Element) FEM solver computes residual stress and deformation from the thermal field. One-way coupling: $T \rightarrow \Delta\varepsilon^{th} \rightarrow \mathbf{u} \rightarrow \boldsymbol{\sigma}$.

### Governing Equation

$$\nabla \cdot \boldsymbol{\sigma} = \mathbf{0}, \quad \Delta\boldsymbol{\sigma} = \mathbf{C} : (\Delta\boldsymbol{\varepsilon} - \alpha_V \Delta T \, \mathbf{I})$$

with J2 plasticity (von Mises yield, temperature-dependent yield strength, radial return mapping).

### FEM Grid

The mechanical grid is coarsened from the thermal grid by `mech_mesh_ratio`. Example: 200x200x50 thermal with ratio=2 gives 100x100x25 FEM nodes. 8-node hexahedral elements with 2x2x2 Gauss quadrature.

### Solution Algorithm

1. Extract temperature from PHOENIX to FEM grid
2. Determine phase (POWDER/LIQUID/SOLID) per node
3. Newton-Raphson iteration ($\text{tol} = 10^{-4}$, max 10):
    - Assemble residual via 8-color EBE
    - Solve $\mathbf{K}\Delta\mathbf{u} = -\mathbf{R}$ via Jacobi-preconditioned CG
4. Update Gauss point stress/strain state with J2 return map
5. Smooth stresses to nodes, compute von Mises

### Parallel Mode

With `bash run.sh <case> N1 N2`, thermal and mechanical run as separate OS processes:

- Thermal writes binary temperature files every `mech_interval` steps
- Mechanical process polls, reads, and solves independently
- ~1.6x speedup on 24-core machine (10+10 threads optimal)

See [Mechanical Solver](../mechanics/overview.md) for full details.
