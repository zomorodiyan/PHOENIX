# Species Transport Analysis

Analysis of `program931/species_transport.f90` — dissimilar metal welding concentration solver.

## Governing Equation

Transient convection-diffusion equation for concentration C (mass fraction of material 1):

```
d(rho*C)/dt + div(rho*u*C) = div(D*grad(C)) + S
```

where D = `massdiffl` (liquid mass diffusivity parameter, units: kg/(m*s)), and C ranges from 0 (pure material 2) to 1 (pure material 1).

## Discretization

**Finite Volume Method** with **Power Law scheme** for combined convection-diffusion:

```
a_nb = D_f * max(0, (1 - 0.1*|F/D|)^5) + max(0, -F)   (upwind neighbor)
a_nb = D_f * max(0, (1 - 0.1*|F/D|)^5) + max(0,  F)   (downwind neighbor)
```

where F = rho*u*A (convective flux) and D = massdiffusivity*A/dx (diffusive conductance).

**Transient term**: Implicit Euler — `acpnot = (rho/dt)*Volume`, added to `su` as `acpnot*C_old`.

**Solver**: Line-by-line TDMA sweeping in x-direction, with double sweeps in j and k (2x2 = 4 total sweeps per call). An additional `enhance_species_speed` subroutine applies a 1D block-correction (whole-field correction) in the x-direction to accelerate convergence.

**Under-relaxation**: Uses `urfh` (enthalpy relaxation factor), shared with the energy equation — should ideally have its own `urfspecies`.

## Mass Diffusivity

```fortran
massdiffusivity(i,j,k) = massdiffl                         ! liquid
if (temp_hyw(i,j,k) <= tsolidmatrix(i,j,k)) massdiffusivity = 1e-11   ! solid (effectively zero)
```

- Uses composition-dependent solidus `tsolidmatrix(i,j,k)` for the solid/liquid switch
- Uses `temp_hyw` (likely a smoothed/old temperature) rather than current `temperature` — this may be intentional for stability

## Source Terms

### Active: Scanning source (x-direction advection due to moving frame)

```fortran
if (scanvel > small) then
    term1 = areajk(j,k) * rhoscan       ! rhoscan = rho * scanvel
    su(i,j,k) += concentration(i-1,j,k) * term1
    sp(i,j,k) -= term1
end if
```

This adds a pseudo-convection in the x-direction representing the moving laser/workpiece frame. The upwind scheme (using `i-1`) implies the scan moves in the +x direction.

Note: `rhoscan` is set to 0 when `delt < 100` (transient mode), so this source is only active in steady-state mode. In practice, since typical `delt` values are O(1e-5), **the scanning source is disabled in transient simulations**.

### Commented out: Solidification segregation

```fortran
! su += -kp*C*(rho/dt)*V*(fl_old - fl) - (rho/dt)*V*((1-fl)*C - (1-fl_old)*C_old)
```

with `kp = 0.16` (equilibrium partition coefficient). This term would model solute rejection during solidification (Scheil-type). Currently **disabled**.

## Initial Conditions

Two-material setup separated along the y-axis:

```fortran
! y < -smoothdistance:  C = 0        (material 2)
! -smoothdistance < y < smoothdistance:  C = (y + smoothdistance) / (2*smoothdistance)  (linear ramp)
! y >= smoothdistance:  C = 1        (material 1)
```

Material properties are then interpolated by concentration (linear mixing rule):

| Property | Formula |
|---|---|
| tsolid | tsolid*C + tsolid2*(1-C) |
| tliquid | tliquid*C + tliquid2*(1-C) |
| tvap | tvap*C + tvap2*(1-C) |
| hlatnt | hlatnt*C + hlatnt2*(1-C) |
| viscosity | viscos*C + viscos2*(1-C) |
| dgdt | dgdtp*C + dgdtp2*(1-C) |
| cp (solid) | acp*C + acp2*(1-C) |
| thermal conductivity (liquid) | thconl*C + thconl2*(1-C) |

## Boundary Conditions

All boundaries use **zero-flux Neumann** (dC/dn = 0):

```fortran
! Top/bottom (k):    C(k=nk) = C(k=nk-1),   C(k=1) = C(k=2)
! North/south (j):   C(j=nj) = C(j=nj-1),   C(j=1) = C(j=2)
! East/west (i):     C(i=ni) = C(i=ni-1),    C(i=1) = C(i=2)
```

This is physically appropriate — no species flux through domain walls.

## Comparison with fortran_new (PHOENIX version)

| Feature | dissimilar (program931) | fortran_new (mod_species.f90) |
|---|---|---|
| Mass diffusivity | `massdiffl` parameter | Hardcoded `5.25e-3` |
| Solidus check | `tsolidmatrix(i,j,k)` (composition-dependent) | `tsolid` (constant) |
| Temperature field | `temp_hyw` | `temp` |
| Velocity factor | No multiplication (1x) | **Velocities multiplied by 2** |
| Scanning source | Active (`rhoscan`) | **Commented out** |
| Segregation source | Commented out (`kp=0.16`) | Active with `kp=1` (no segregation) |
| Top BC (k=nk) | Zero-flux | **Solubility BC**: `C = 2*(13-C_interior)/(dz*5.25e-3) + C_interior` when T > tliquid |
| Bottom BC (k=1) | Zero-flux `C(1)=C(2)` | **Fixed `C=0`** |
| OpenMP | Not parallelized (species subroutine) | Parallelized |
| Under-relaxation | Shared `urfh` | Shared `urfh` |

## Issues and Potential Bugs

### 1. Velocity x2 factor (fortran_new only)
The legacy/fortran_new version multiplies all velocities by 2 before computing convective fluxes. This has **no physical justification** and will produce incorrect convection. The dissimilar version correctly uses raw velocities. **Status: Fixed in dissimilar.**

### 2. Scanning source disabled in transient mode
`rhoscan` is set to 0 when `delt < 100` (i.e., all transient runs). The scanning source term is therefore **never active** in typical simulations. If the moving-frame formulation is intended for transient problems, this logic needs revisiting.

### 3. Segregation source commented out
The solidification segregation term with `kp=0.16` is fully commented out. Without it, the solver treats the mushy zone as a simple mixture — no solute rejection during solidification. This is a simplification that may be acceptable for initial runs but will miss important physics (microsegregation, constitutional undercooling).

### 4. Shared under-relaxation factor
Species uses `urfh` (the enthalpy relaxation factor). Species transport typically converges differently from energy — a dedicated `urfspecies` parameter would allow independent tuning. The code even has a comment acknowledging this: `! for species transport, this coefficient (urfspecies) can be changed accordingly.`

### 5. Constant density in convective fluxes
The convective fluxes use a single `dens` value rather than a composition-weighted density. For dissimilar metals with significantly different densities, this could introduce mass conservation errors. The `densmatrix(i,j,k)` array exists but is not used in species transport.

### 6. Mass diffusivity in solid
Set to `1e-11` — effectively zero but not exactly zero. This is fine numerically but creates a sharp diffusivity discontinuity at the solidus. The mushy zone (between solidus and liquidus) uses the full liquid diffusivity, which may overestimate mixing in partially solidified regions.

### 7. No concentration bounds enforcement
There is no clipping of concentration to [0, 1]. With numerical errors, especially from the `enhance_species_speed` block correction, values could drift outside the physical range. This could cause issues with concentration-weighted material properties (negative or >1 values are unphysical).

### 8. Property update timing
Material properties (tsolidmatrix, viscosity, etc.) are initialized based on the initial concentration field in `mod_init.f90`. It is not clear from the species transport module alone whether these properties are re-updated each time step as concentration evolves. If not, the simulation uses stale material properties. **Needs verification in the main loop.**

### 9. dgdc computation
```fortran
dgdc(i,j,k) = (dgconst - dgconst2)  ! *concentration...
```
The concentration-dependent part is commented out, making dgdc a constant. This means the solutal Marangoni effect uses a fixed concentration gradient coefficient regardless of local composition — a simplification.
