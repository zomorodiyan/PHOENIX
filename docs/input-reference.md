# Input File Reference

All simulation parameters are defined in `fortran_new/inputfile/input_param.txt`. The file uses Fortran namelist format.

## File Structure

The input file has two sections:

1. **Geometry block** (free-format, line-by-line)
2. **Namelist blocks** (`&name ... /` format)

## Geometry Parameters

The geometry is defined by zones in each direction. Each zone has a length, number of control volumes, and power-law exponent for grid stretching.

```
nzx                          ! number of x-zones (max 7)
xzone(1), xzone(2), ...     ! length of each zone (m)
ncvx(1), ncvx(2), ...       ! number of CVs in each zone
powrx(1), powrx(2), ...     ! power-law exponents (1=uniform, >1=cluster at end, <0=cluster at start)
```

Same pattern repeated for y (`nzy`) and z (`nzz`).

**Example** (current setup):

| Direction | Zones | Lengths | CVs | Exponents | Total |
|-----------|-------|---------|-----|-----------|-------|
| x | 1 | 4.0 mm | 400 | 1 (uniform) | 400 cells |
| y | 1 | 4.0 mm | 400 | 1 (uniform) | 400 cells |
| z | 2 | 0.5 mm + 0.2 mm | 10 + 40 | -1.5, 1 | 50 cells |

The z-direction uses two zones: a substrate zone (0.5 mm, coarse, clustered toward top) and a powder/build zone (0.2 mm, fine, uniform). Total domain: 4 mm x 4 mm x 0.7 mm.

!!! tip "Grid Stretching"
    Power-law exponent controls cell clustering:

    - `powr = 1`: uniform spacing
    - `powr > 1`: cells cluster toward the end of the zone
    - `powr < 0`: cells cluster toward the start of the zone (e.g., -1.5 clusters z-cells toward the substrate-powder interface)

## Namelist Blocks

### `&process_parameters` — Surface Laser Source

| Parameter | Unit | Description |
|-----------|------|-------------|
| `alaspow` | W | Surface laser power (set to 0 when using volumetric source) |
| `alaseta` | - | Surface laser absorption efficiency |
| `alasrb` | m | Laser beam radius (1/e^2) |
| `alasfact` | - | Gaussian factor (typically 2 for 1/e^2 definition) |

### `&volumetric_parameters` — Volumetric Heat Source

| Parameter | Unit | Description |
|-----------|------|-------------|
| `alaspowvol` | W | Volumetric laser power |
| `alasetavol` | - | Volumetric absorption efficiency |
| `sourcerad` | m | Source radius (Gaussian 1/e^2) |
| `sourcedepth` | m | Source penetration depth |

!!! note
    Either surface (`alaspow`) or volumetric (`alaspowvol`) source is used. Set the unused one to 0. The volumetric source distributes heat as:

    $$q(x,y,z) = \frac{P \cdot \eta \cdot f}{\pi r_s^2 \cdot d_s} \exp\left(-\frac{f}{r_s^2}\left[(x-x_b)^2 + (y-y_b)^2\right]\right)$$

### `&material_properties` — Primary Material (C=1)

| Parameter | Unit | Description | Example |
|-----------|------|-------------|---------|
| `dens` | kg/m^3 | Solid density | 8440 |
| `denl` | kg/m^3 | Liquid density | 7640 |
| `viscos` | Pa*s | Dynamic viscosity (liquid) | 0.007 |
| `tsolid` | K | Solidus temperature | 1563 |
| `tliquid` | K | Liquidus temperature | 1623 |
| `tboiling` | K | Boiling/vaporization temperature | 2650 |
| `hlatnt` | J/kg | Latent heat of fusion | 290000 |
| `acpa` | J/(kg*K^2) | Solid specific heat coefficient (quadratic term) | 0.2441 |
| `acpb` | J/(kg*K) | Solid specific heat coefficient (linear term) | 338.59 |
| `acpl` | J/(kg*K) | Liquid specific heat (constant) | 709.25 |
| `thconsa` | W/(m*K^2) | Solid thermal conductivity coefficient (linear term) | 0.0155 |
| `thconsb` | W/(m*K) | Solid thermal conductivity coefficient (constant term) | 5.0435 |
| `thconl` | W/(m*K) | Liquid thermal conductivity (constant) | 30.078 |
| `beta` | 1/K | Thermal expansion coefficient (for buoyancy) | 5e-5 |
| `emiss` | - | Surface emissivity (for radiative loss) | 0.3 |
| `dgdtp` | N/(m*K) | dg/dT — thermal Marangoni coefficient | -3.8e-4 |

!!! info "Specific Heat Model"
    Solid specific heat is temperature-dependent: $c_p(T) = a \cdot T + b$ where `acpa` = $a$, `acpb` = $b$.
    The enthalpy-temperature curve is:

    - **Solid** ($T \leq T_s$): $H = \frac{a}{2}T^2 + bT$
    - **Mushy** ($T_s < T < T_l$): linear interpolation between $H_s$ and $H_l$
    - **Liquid** ($T \geq T_l$): $H = H_l + c_{p,l}(T - T_l)$

### `&powder_properties` — Powder Layer

| Parameter | Unit | Description | Example |
|-----------|------|-------------|---------|
| `layerheight` | m | Powder layer thickness | 0.04e-3 |
| `pden` | kg/m^3 | Powder density (accounts for porosity) | 4330 |
| `pcpa` | J/(kg*K^2) | Powder specific heat (quadratic) | 0.2508 |
| `pcpb` | J/(kg*K) | Powder specific heat (linear) | 357.7 |
| `pthcona` | W/(m*K^2) | Powder conductivity (linear) | 0 |
| `pthconb` | W/(m*K) | Powder conductivity (constant) | 0.995 |

Powder properties apply to cells in the top `layerheight` of the domain that have not yet been melted (`solidfield <= 0.5`).

### `&numerical_relax` — Solver Parameters

| Parameter | Unit | Description | Example |
|-----------|------|-------------|---------|
| `maxit` | - | Max iterations per time step | 60 |
| `delt` | s | Time step size | 2e-5 |
| `timax` | s | Total simulation time | 0.09 |
| `urfu` | - | Under-relaxation factor for u-velocity | 0.7 |
| `urfv` | - | Under-relaxation factor for v-velocity | 0.7 |
| `urfw` | - | Under-relaxation factor for w-velocity | 0.7 |
| `urfp` | - | Under-relaxation factor for pressure | 0.7 |
| `urfh` | - | Under-relaxation factor for enthalpy | 0.7 |

### `&boundary_conditions` — Thermal Boundaries

| Parameter | Unit | Description | Example |
|-----------|------|-------------|---------|
| `htci` | W/(m^2*K) | Convection coefficient on x-faces (west/east) | 5 |
| `htcj` | W/(m^2*K) | Convection coefficient on y-faces (north/south) | 5 |
| `htck1` | W/(m^2*K) | Convection coefficient on bottom face | 5 |
| `htckn` | W/(m^2*K) | Convection coefficient on top face | 5 |
| `tempWest` | K | Far-field temperature at west boundary | 293 |
| `tempEast` | K | Far-field temperature at east boundary | 293 |
| `tempNorth` | K | Far-field temperature at north boundary | 293 |
| `tempBottom` | K | Far-field temperature at bottom boundary | 293 |
| `tempPreheat` | K | Initial/preheat temperature | 293 |
| `tempAmb` | K | Ambient temperature (for radiation) | 293 |

The top surface uses combined convection + radiation:
$$q_{loss} = h(T - T_{amb}) + \epsilon \sigma (T^4 - T_{amb}^4)$$

### `&local_solver` — Adaptive Solver

| Parameter | Unit | Description | Example |
|-----------|------|-------------|---------|
| `localnum` | - | Number of local steps between global steps | 4 |
| `local_half_x` | m | Half-width of local region in x | 1.0e-3 |
| `local_half_y` | m | Half-width of local region in y | 2.0e-4 |
| `local_depth_z` | m | Depth of local region in z | 2.0e-4 |

The local solver only updates enthalpy in a small region around the melt pool for `localnum` consecutive steps, then solves the full domain. This provides significant speedup (typically 3-5x) with minimal accuracy loss.

### `&output_control` — Output Settings

| Parameter | Unit | Description | Example |
|-----------|------|-------------|---------|
| `outputintervel` | steps | VTK output frequency (every N time steps) | 50 |
| `case_name` | - | Case identifier (sets output directory name) | 'testcase' |
| `toolpath_file` | - | Path to toolpath file (.crs) | './ToolFiles/B26.crs' |
| `species_flag` | - | Enable species transport (0=off, 1=on) | 0 |
| `micro_flag` | - | Enable solidification microstructure prediction (0=off, 1=on) | 0 |
| `crack_flag` | - | Enable crack risk prediction (0=off, 1=on) | 0 |

### `&microstructure_params` — Solidification Microstructure (optional)

Only read when `micro_flag=1`. All parameters have defaults for IN718.

| Parameter | Unit | Description | Default |
|-----------|------|-------------|---------|
| `a1_pdas` | m | PDAS prefactor | 50e-6 |
| `a2_sdas` | m | SDAS prefactor | 10e-6 |
| `n1_pdas` | - | PDAS exponent for G | -0.5 |
| `n2_pdas` | - | PDAS exponent for R | -0.25 |
| `n3_sdas` | - | SDAS exponent for cooling rate | -0.333 |

PDAS: $\lambda_1 = a_1 \cdot G^{n_1} \cdot R^{n_2}$, SDAS: $\lambda_2 = a_2 \cdot \dot{T}^{n_3}$

### `&crack_params` — Crack Risk Prediction (optional)

Only read when `crack_flag=1`.

| Parameter | Unit | Description | Default |
|-----------|------|-------------|---------|
| `delta_t_btr` | K | Brittle Temperature Range width below T_solidus | 100.0 |

CSI = $\int_{BTR} \alpha \cdot |\dot{T}| \, dt$ where BTR = $[T_s - \Delta T_{BTR}, \, T_s]$
