# Mechanical Solver

## Overview

The mechanical solver computes residual stress and deformation in laser powder bed fusion using an Element-By-Element (EBE) finite element method with J2 elastoplasticity. It provides one-way coupling from the thermal solver: temperature drives thermal strain, which produces displacement and stress through quasi-static mechanical equilibrium.

The solver is implemented in `fortran_new/mechanical_solver/` as three modules:

| Module | File | Purpose |
|--------|------|---------|
| `mech_material` | `mod_mech_material.f90` | Material properties, phase logic, J2 return map |
| `mechanical_solver` | `mod_mechanical.f90` | FEM solver (Newton-Raphson + CG), EBE assembly |
| `mech_io` | `mod_mech_io.f90` | VTK output, timing, history, deformation GIF |

## Enabling the Mechanical Solver

In `input_param.txt`:

```
&mechanical_params mechanical_flag=1, mech_interval=25, mech_output_interval=3, mech_mesh_ratio=2 /
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mechanical_flag` | 0 | Enable mechanical solver (0=off, 1=on) |
| `mech_interval` | 10 | Solve every N thermal steps |
| `mech_output_interval` | 5 | Write VTK every N mechanical solves |
| `mech_mesh_ratio` | 2 | FEM grid coarsening ratio (1=same as thermal, 2=half) |

When `mechanical_flag=0`, the mechanical solver is completely inert.

## Physics

### Governing Equation

Quasi-static mechanical equilibrium (inertia neglected):

$$\nabla \cdot \boldsymbol{\sigma} = \mathbf{0}$$

with the constitutive relation for incremental stress:

$$\Delta\boldsymbol{\sigma} = \mathbf{C} : (\Delta\boldsymbol{\varepsilon} - \Delta\boldsymbol{\varepsilon}^{th})$$

where $\mathbf{C}$ is the isotropic elasticity tensor and the thermal strain increment is:

$$\Delta\varepsilon^{th}_{ij} = \alpha_V \Delta T \, \delta_{ij}$$

### Isotropic Elasticity

The 6x6 elasticity matrix in Voigt notation:

$$\mathbf{C} = \begin{bmatrix} \lambda+2\mu & \lambda & \lambda & 0 & 0 & 0 \\ \lambda & \lambda+2\mu & \lambda & 0 & 0 & 0 \\ \lambda & \lambda & \lambda+2\mu & 0 & 0 & 0 \\ 0 & 0 & 0 & \mu & 0 & 0 \\ 0 & 0 & 0 & 0 & \mu & 0 \\ 0 & 0 & 0 & 0 & 0 & \mu \end{bmatrix}$$

with Lame parameters:

$$\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}, \quad \mu = \frac{E}{2(1+\nu)}$$

### J2 Plasticity (von Mises)

The yield condition uses the von Mises criterion with temperature-dependent yield strength:

$$f = \|\mathbf{s}\| - \sigma_y(T) \leq 0$$

where $\mathbf{s}$ is the deviatoric stress and $\sigma_y(T)$ decreases linearly from $\sigma_{y,0}$ at $T_{ref}$ to 0 at $T_{solidus}$:

$$\sigma_y(T) = \sigma_{y,0} \cdot \frac{T_s - T}{T_s - T_{ref}}$$

When $f > 0$, a radial return mapping projects the trial stress back onto the yield surface:

$$\boldsymbol{\sigma} = s_{mean}\mathbf{I} + \mathbf{s}_{trial} \cdot \frac{\sigma_y}{\|\mathbf{s}_{trial}\|}$$

### Phase-Dependent Properties

Each FEM node is classified into one of three mechanical phases:

| Phase | Condition | Young's Modulus | Thermal Expansion |
|-------|-----------|-----------------|-------------------|
| **SOLID** (2) | Substrate or solidified ($T < T_s$ and $sf > 0.5$) | $E_{solid}$ = 70 GPa | $\alpha_V$ = 1e-5 /K |
| **LIQUID** (1) | Currently molten ($T \geq T_s$) | $E_{soft}$ = 0.7 GPa | 0 |
| **POWDER** (0) | In powder layer, never melted | $E_{soft}$ = 0.7 GPa | 0 |

Powder and liquid use a soft modulus ($E_{solid}/100$) rather than zero to maintain numerical stability.

## Material Parameters

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| $E_{solid}$ | 70 | GPa | Young's modulus (solid/substrate) |
| $E_{soft}$ | 0.7 | GPa | Young's modulus (powder/liquid) |
| $\nu$ | 0.3 | - | Poisson's ratio |
| $\sigma_{y,0}$ | 250 | MPa | Yield stress at reference temperature |
| $T_{ref}$ | 293 | K | Reference temperature for yield |
| $\alpha_V$ | 1e-5 | 1/K | Volumetric thermal expansion |

These are defined as `parameter` constants in `mod_mech_material.f90`.

## Solver Details

### FEM Discretization

- **Element type**: 8-node hexahedral (trilinear brick)
- **Quadrature**: 2x2x2 Gauss points per element
- **DOFs**: 3 per node (ux, uy, uz), 24 per element

### FEM Grid

The mechanical grid is coarsened from the thermal grid by `mech_mesh_ratio`:

- FEM node $i$ maps to thermal cell center $x(2 + (i-1) \cdot r)$ where $r$ = `mech_mesh_ratio`
- Last node forced to domain boundary for full coverage
- Example: 200x200x50 thermal with ratio=2 gives 100x100x25 FEM nodes

### EBE Solver

The Element-By-Element approach avoids assembling a global stiffness matrix:

1. **Newton-Raphson** outer loop (tolerance $10^{-4}$, max 10 iterations)
2. **Jacobi-preconditioned CG** inner solver (tolerance $10^{-4}$, max 20000 iterations)
3. **8-color element coloring** for thread-safe OpenMP parallelism in scatter operations
4. **Matrix-free matvec** for non-uniform grids (AMR): computes $\mathbf{B}^T\mathbf{C}\mathbf{B}\mathbf{u}$ per Gauss point without forming $\mathbf{K}_e$

### Boundary Conditions

- **Bottom face** ($k=1$): $\mathbf{u} = \mathbf{0}$ (clamped substrate base)
- **All other faces**: traction-free (natural BC)

### Solution Algorithm

Each mechanical solve:

1. Extract temperature and solidfield from PHOENIX to FEM grid
2. Determine mechanical phase (POWDER / LIQUID / SOLID) per node
3. Compute incremental temperature change $\Delta T$ at Gauss points
4. Newton iteration:
    - Compute residual $\mathbf{R} = \mathbf{B}^T \boldsymbol{\sigma} \cdot J$
    - Solve $\mathbf{K} \, \Delta\mathbf{u} = -\mathbf{R}$ via CG
    - Update displacement $\mathbf{u} \leftarrow \mathbf{u} + \Delta\mathbf{u}$
5. Update Gauss point stress/strain state (with J2 return map)
6. Smooth stresses to nodes, compute von Mises

## AMR Compatibility

When adaptive mesh refinement is enabled, the mechanical solver handles grid changes:

### Grid Update (`update_mech_grid`)

Called after each AMR remesh. Updates:

1. **FEM coordinates**: refresh `fem_x`, `fem_y` from new PHOENIX grid
2. **Precomputed Ke**: recompute per-layer stiffness matrices
3. **Field interpolation** (nearest-neighbor to avoid diffusion):
    - Displacement: `ux_mech`, `uy_mech`, `uz_mech`
    - Gauss point stress: `sig_gp`
    - Yield function: `f_plus`
    - Reference temperature: `T_old_mech`
4. **Strain recomputation**: `eps_gp = B_{new} \cdot u_{interp}` (ensures consistency)

### Matrix-Free Matvec

For non-uniform grids (where dx/dy vary per element), the CG matvec uses a matrix-free approach:

- **Uniform grid**: uses precomputed per-layer $\mathbf{K}_e$ (fast)
- **Non-uniform grid**: computes $\mathbf{f}_e = \sum_{g=1}^{8} \mathbf{B}_g^T \, \mathbf{C} \, \mathbf{B}_g \, \mathbf{u}_e \cdot J_g$ per element (~5x faster than assembling $\mathbf{K}_e$ on-the-fly)

Grid uniformity is detected automatically (1% tolerance) and cached in a `grid_uniform` flag.

### Nearest-Neighbor Interpolation

All mechanical fields use nearest-neighbor (not bilinear) interpolation after AMR remesh. Bilinear interpolation acts as a low-pass filter, causing cumulative stress/displacement diffusion over repeated AMR cycles. Nearest-neighbor preserves field magnitudes exactly with only sub-element spatial offset.

## Thermal-Mechanical Coupling

### Serial Mode (default)

```bash
bash run.sh <case> <threads> &
# Example: bash run.sh baseline 10 &
```

Mechanical solver runs in-loop every `mech_interval` thermal steps. Appears in the thermal timing report.

### Parallel Mode

```bash
bash run.sh <case> <thermal_threads> <mech_threads> &
# Example: bash run.sh baseline 10 10 &
```

Thermal and mechanical run as **separate OS processes** with file-based communication:

- Thermal writes binary input files (`_mech_input_NNNNN.bin`) at each `mech_interval`
- Mechanical process polls for files, reads, solves, writes VTK/history
- Mechanical computes its termination step from `timax`, `delt`, `mech_interval`
- No sentinel files needed

**Binary file contents**: step index, time, grid dimensions, temperature field, solidfield, x/y coordinates (for AMR).

**Thread allocation**: each process gets independent OpenMP threads. Optimal split depends on grid size; on a 24-core machine, 10+10 achieved **1.63x speedup** over serial.

## VTK Output

Mechanical results are written to separate VTK files: `<case>_mech_NNNNN.vtk`

| Field | Type | Description |
|-------|------|-------------|
| `Temperature` | scalar | Temperature at FEM nodes (K) |
| `ux`, `uy`, `uz` | scalar | Displacement components (m) |
| `phase` | integer | 0=powder, 1=liquid, 2=solid |
| `sxx`, `syy`, `szz` | scalar | Normal stress components (Pa) |
| `von_mises` | scalar | Von Mises equivalent stress (Pa) |
| `fplus` | scalar | Yield function (>0 = plastic yielding) |

### Deformation GIF

After simulation, `finalize_mechanical_io` generates a Python script that creates:

- `<case>_deformation.gif`: von Mises stress animation with 10x deformation magnification
- `<case>_deformation_final.png`: high-resolution final frame

Uses matplotlib for rendering (no GPU needed) and imageio for GIF encoding.

## Performance

Typical metrics for a 100x50x25 FEM grid with AMR:

| Metric | Value |
|--------|-------|
| Mechanical fraction (serial) | ~52% of total wall time |
| Avg wall time per solve | ~20 s (10 threads) |
| FEM memory (sig_gp + eps_gp + u + stress) | ~24 MB |
| VTK file size | ~5 MB each |
| Parallel speedup (10+10 threads) | 1.63x |

## History Output

Mechanical history is written every solve step to `<case>_mech_history.txt`:

- 10 monitoring points (same as thermal history)
- Fields: time, T, ux, uy, uz, sxx, syy per point
- Auto-generates `<case>_mech_history.png` with 3 subplots (temperature, displacement, stress)
