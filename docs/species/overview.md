# Species Transport

## Overview

The species transport module (`mod_species.f90`) enables simulation of dissimilar metal mixing in additive manufacturing. It solves a convection-diffusion equation for concentration and provides composition-dependent material properties through two-way coupling.

## Concentration Convention

| Value | Meaning | Properties Source |
|-------|---------|-------------------|
| **C = 1** | Pure primary (base) material | `&material_properties` in `input_param.txt` |
| **C = 0** | Pure secondary material | Named constants in `mod_species.f90` |
| **0 < C < 1** | Mixed region | `mix(prop1, prop2, C) = prop1*C + prop2*(1-C)` |

## Enabling Species Transport

In `input_param.txt`:

```
&output_control ... species_flag=1 /
```

When `species_flag=0` (default), the species module is completely inert — no arrays allocated, no computation, no overhead.

## Physics

### Transport Equation

$$\frac{\partial (\rho C)}{\partial t} + \nabla \cdot (\rho \vec{u} C) = \nabla \cdot (\Gamma \nabla C)$$

where $\Gamma = \rho D_m$ is the mass diffusion coefficient, $D_m = 5 \times 10^{-9}$ m$^2$/s.

In solid regions ($T < T_s$), $\Gamma$ is set to a floor value ($10^{-30}$) to effectively freeze concentration.

### Two-Way Coupling

When `species_flag=1`, material properties depend on local concentration:

- **Thermal**: `tsolid`, `tliquid`, `acpa`, `acpb`, `acpl`, `thconsa`, `thconsb`, `thconl`
- **Mechanical**: `dens`, `denl`, `viscos`
- **Powder**: `pden`, `pcpa`, `pcpb`, `pthcona`, `pthconb`

All computed inline via `mix()` — no extra arrays needed.

Properties that remain scalar (same for both materials):

- `beta` (thermal expansion), `emiss` (emissivity)
- `hlatnt` (latent heat), `dgdt` (thermal Marangoni)

### Solutal Marangoni

Surface tension depends on concentration: $\gamma = \gamma_0 + \frac{d\gamma}{dT}(T - T_{ref}) + \frac{d\gamma}{dC}(C - C_{ref})$

The solutal Marangoni stress drives flow from low-C to high-C regions (when `dgdc < 0`):

$$\tau_C = \frac{d\gamma}{dC} \frac{\partial C}{\partial x}$$

Currently `dgdc_const = -0.3` N/m (typical for liquid metal binary alloys).

### Initial Conditions

- **Substrate** (below powder layer): C = 1 (base material)
- **Powder layer, y < y_mid**: C = 1 (base material)
- **Powder layer, y >= y_mid**: C = 0 (secondary material)

This represents a half-and-half powder bed configuration for dissimilar metal testing, with the secondary material in the upper half of the y-domain.

## Secondary Material Properties

Defined as named constants in `mod_species.f90`:

| Property | Primary | Secondary | Unit |
|----------|---------|-----------|------|
| `dens/dens2` | 8440 | 8880 | kg/m^3 |
| `denl/denl2` | 7640 | 7800 | kg/m^3 |
| `viscos/viscos2` | 0.007 | 0.003 | Pa*s |
| `tsolid/tsolid2` | 1563 | 1728 | K |
| `tliquid/tliquid2` | 1623 | 1803 | K |
| `acpa/acpa2` | 0.2441 | 0.3441 | J/(kg*K^2) |
| `acpb/acpb2` | 338.59 | 400 | J/(kg*K) |
| `acpl/acpl2` | 709.25 | 800 | J/(kg*K) |
| `thconl/thconl2` | 30.078 | 120 | W/(m*K) |
| `D_m` | — | 5e-9 | m^2/s |
| `dgdc_const` | — | -0.3 | N/m |

## Solver Details

- Discretization: Power-law scheme (same as enthalpy)
- Solver: Line-by-line TDMA with block correction
- Domain: Melt pool region only (`istatp1:iendm1, jstat:jend, kstat:nkm1`)
- Timing: Once per timestep, after iteration loop exits
- Under-relaxation: `urfspecies = 0.7`
- Clipping: concentration forced to [0, 1] after each solve
- Transient: uses `delt`

## Performance

Typical overhead with species enabled:

| Metric | Impact |
|--------|--------|
| Species solver CPU | < 1% (solves only in melt pool) |
| Two-way coupling overhead | ~16% (mix() calls in properties, enthalpy_to_temp, source) |
| Memory | +65 MB (two concentration arrays on 400x400x50 grid) |

## VTK Output

When `species_flag=1`, two additional scalar fields appear in `vtkmov*.vtk`:

- **concentration**: mass fraction of primary material [0, 1]
- **tsolid_field**: composition-dependent solidus temperature (K)
