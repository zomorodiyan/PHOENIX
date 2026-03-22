# Task: Crack Risk Prediction from Temperature Field

## Objective
Compute crack risk indicators from the temperature field during simulation — no mechanical solver needed. Predicts **where** cracks are most likely to form based on thermal strain rate and time spent in the Brittle Temperature Range (BTR).

### Motivation
- Cracking (solidification cracking, strain-age cracking) is a major defect mode in AM
- Full thermo-mechanical FEA is expensive and requires a separate solver
- The temperature field alone provides surprisingly good crack risk indicators:
  - High cooling rate near solidus → high thermal strain rate → solidification cracking
  - Long time in BTR → ductility-dip cracking
- These indicators have been validated in literature (Kou criterion, BTR analysis)

---

## Physics Background

### Thermal Strain Rate
```
dε/dt = α × |dT/dt|
```
where α = thermal expansion coefficient (`beta` in input), dT/dt = cooling rate.
High thermal strain rate during solidification → material contracts faster than liquid can feed → cracking.

### Brittle Temperature Range (BTR)
The temperature range [T_BTR_low, T_solidus] where the material has near-zero ductility:
- Solid fraction is high (>0.9) but grain boundaries still have thin liquid films
- Cannot accommodate strain → intergranular cracking
- `T_BTR_low` is alloy-specific, typically `T_solidus - ΔT_BTR` (ΔT_BTR ≈ 50-150 K for Ni-alloys)

### Crack Susceptibility Index (CSI)
```
CSI = ε_BTR = ∫(over BTR) α × |dT/dt| dt
```
Total accumulated thermal strain while the cell is in the BTR. Higher CSI → higher crack risk.
This is dimensionless and directly comparable across process conditions.

---

## Design Decisions

### Arrays (allocated in mod_crack_risk.f90)
- `cool_rate_solid(ni,nj,nk)` — cooling rate at moment of solidification (K/s)
- `strain_rate_solid(ni,nj,nk)` — thermal strain rate at solidification (1/s)
- `btr_time(ni,nj,nk)` — cumulative time spent in BTR (s)
- `crack_risk(ni,nj,nk)` — crack susceptibility index CSI (dimensionless)
- `solidified(ni,nj,nk)` — logical mask: has this cell solidified?
- `in_btr(ni,nj,nk)` — logical mask: is cell currently in BTR?
- Memory: ~4 real + 2 logical arrays on layer region ≈ minimal overhead

### Input parameters
Add `crack_flag` to `&output_control` namelist:
```
&output_control  ..., crack_flag=1 /
```

Add new namelist `&crack_params`:
```
&crack_params  delta_t_btr=100.0 /
```
- `crack_flag`: 0 = disabled (default), 1 = enabled
- `delta_t_btr`: BTR width below T_solidus (K), default 100 K for IN718

### Z range
Same as defect module: only the build layer (`k_def_lo:k_def_hi`).

### Solidification detection
A cell solidifies when `fracl` transitions from >0 to 0 (liquid fraction drops to zero). Only record the **first** solidification event (ignore remelting cycles for now).

---

## Tasks

### Task 1 — Create mod_crack_risk.f90

New module following mod_defect.f90 pattern:

```fortran
module crack_risk_mod
    use precision
    use geometry
    use initialization
    use parameters
    use sim_state
    use constant

    implicit none

    real(wp) :: delta_t_btr = 100.0_wp   ! BTR width (K)
    real(wp), allocatable :: cool_rate_solid(:,:,:)
    real(wp), allocatable :: strain_rate_solid(:,:,:)
    real(wp), allocatable :: btr_time(:,:,:)
    real(wp), allocatable :: crack_risk_arr(:,:,:)
    logical, allocatable  :: solidified(:,:,:)
    logical, allocatable  :: in_btr(:,:,:)

contains

    subroutine allocate_crack_risk(nni, nnj, nnk)
    subroutine update_crack_risk(dt)
    subroutine compute_crack_report()
    subroutine write_crack_vtk(fieldname, field)

end module
```

**Subroutines:**

1. `allocate_crack_risk(nni, nnj, nnk)`:
   - Allocate all arrays, initialize to zero/false
   - Use same k_def_lo:k_def_hi as defect module

2. `update_crack_risk(dt)`:
   - Called every timestep after enthalpy solve
   - For each cell in layer:
     - Compute `dTdt = (temp - tnot) / dt`
     - **Solidification detection**: if `fraclnot > 0` and `fracl == 0` and `.not. solidified`:
       - Store `cool_rate_solid = |dTdt|`
       - Store `strain_rate_solid = beta * |dTdt|`
       - Set `solidified = .true.`
     - **BTR tracking**: if `temp < tsolid` and `temp > tsolid - delta_t_btr`:
       - `btr_time += dt`
       - `crack_risk += beta * |dTdt| * dt` (accumulated strain in BTR)
       - Set `in_btr = .true.`
     - Else: `in_btr = .false.`

3. `compute_crack_report()`:
   - Post-simulation: compute statistics over scanned region
   - Write `crack_report.txt`: max/mean CSI, max cooling rate, high-risk fraction
   - Write VTK files: `crack_csi`, `cool_rate_solid`, `strain_rate_solid`, `btr_time`

4. `write_crack_vtk(fieldname, field)`:
   - Same pattern as `write_defect_vtk` in mod_defect.f90

### Task 2 — Add crack_flag to Input System

**mod_param.f90:**
- Add `integer crack_flag` and `real(wp) delta_t_btr`
- Add to namelist: `&output_control ..., crack_flag /`
- Add new namelist: `&crack_params delta_t_btr /`
- Default: `crack_flag = 0`, `delta_t_btr = 100.0`

**input_param.txt:**
- Add `crack_flag=0` to `&output_control`
- Add new line: `&crack_params delta_t_btr=100.0 /`

### Task 3 — Integrate into main.f90

```fortran
use crack_risk_mod

! After allocate_defect:
if (crack_flag == 1) call allocate_crack_risk(ni, nj, nk)

! After update_max_temp in time_loop:
if (crack_flag == 1) call update_crack_risk(delt)

! After write_defect_report post-simulation:
if (crack_flag == 1) call compute_crack_report()
```

### Task 4 — Add to compile.sh

Insert `mod_crack_risk.f90` after `mod_defect.f90` in the compile order.

### Task 5 — Testing

- Enable `crack_flag=1` in input_param.txt
- Run a test case (e.g., 2x2mm toolpath)
- Verify:
  - `cool_rate_solid` values are O(10⁴–10⁶ K/s) for LPBF
  - `crack_risk` (CSI) highest at melt pool boundaries
  - `btr_time` is short (O(10⁻⁴ s)) for fast-moving laser
  - VTK output loads correctly in ParaView
  - Zero overhead when `crack_flag=0`

### Task 6 — Documentation

- Update `docs/input-reference.md` with `crack_flag` and `&crack_params`
- Update `docs/output.md` with new VTK fields and crack report
- Update `docs/architecture/modules.md` with `mod_crack_risk.f90`
- Update `docs/index.md` key features list

---

## Notes
- **Independent of microstructure module**: `crack_flag` and `micro_flag` are fully independent. Both compute cooling rate internally from `(temp - tnot) / delt` — no shared arrays or dependencies. Can enable either one, both, or neither.
- No feedback into solver — purely diagnostic output
- `beta` (thermal expansion coefficient) already available in `mod_param.f90`
- No Young's modulus needed — CSI is a strain-based indicator, not stress
- Future extension: add Kou criterion (|dT/d(fs^1/2)|) using fracl data
- Future extension: combine with species transport for composition-dependent BTR

---

## Completion Log (2026-03-21)

All tasks implemented and validated:

- **mod_crack_risk.f90** created with `allocate_crack_risk`, `update_crack_risk`, `compute_crack_report`, `write_crack_vtk`
- **crack_flag** added to `&output_control` (default 0)
- **&crack_params** namelist added (delta_t_btr)
- Integrated into main.f90, compile.sh
- VTK output: crack_csi, cool_rate_solid, strain_rate_solid, btr_time

### Validation Results (2x2mm, IN718, 300W volumetric)
- Max CSI: 0.0695 (7% strain) — physically reasonable
- Mean CSI: 0.0098 (~1%)
- Max cooling rate at solidification: 5.32 MK/s ✓
- Max strain rate: 265.8 /s ✓ (beta × cooling rate)
- High-risk cells (CSI > 1%): 77,120 / 297,003 = 26%
- delta_T_BTR = 100 K, beta = 5e-5
