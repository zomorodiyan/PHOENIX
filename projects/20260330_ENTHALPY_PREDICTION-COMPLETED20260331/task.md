# Enthalpy Prediction for Heating Steps

## Objective

During laser-on (heating) timesteps, shift all fields (enthalpy, velocity, pressure) along the scan direction by integer number of cells within the melt pool region. This provides a near-perfect initial guess for the iterative solver, reducing iteration count by ~38%.

## Status: ACCEPTED (integer-cell shift)

## Execution Tracking

- **`log.md`**: Execution log with system timestamps.
- **`results.md`**: Before/after comparison.

---

## Approaches Tested

### 1. Bilinear interpolation shift — REJECTED
- Shifts field by `v*dt` using bilinear interpolation
- Introduces numerical diffusion at solid-liquid interface → more mushy zone cells → heavier source terms
- Result: 7% SLOWER (enthalpy-only), 5.5% SLOWER (local all-field)

### 2. Integer-cell shift, local region — ACCEPTED
- Shifts all fields by `nint(v*dt/dx)` cells using index copy (no interpolation)
- Preserves sharp interfaces exactly
- Only within melt pool region (from `pool_size` istat:iend + extension)
- Result: **31% faster wall time, 38% fewer iterations**

---

## Implementation

### Task 1: `predict_flag` in `mod_param.f90`
- Added `integer :: predict_flag = 0` (0=off, 1=integer-cell shift)
- Added to `namelist / output_control /`

### Task 2: `mod_predict.f90` (separate module for future extensibility)
- `predict_shift_integer(vx, vy, dt, is, ie, js, je, ks, ke)`: shifts h, u, v, w, p
- `ishift_field(field, di, dj, is, ie, js, je, ks, ke)`: integer-cell shift with direction-aware iteration order
- Future: could be replaced with ML-based prediction

### Task 3: Call from `main.f90`
- After `pool_size`, before `iter_loop`
- Only when `predict_flag==1` and laser is on and tpeak > tsolid
- Shift region extended by `nint(v*dt/dx)` cells beyond melt pool in scan direction
- `enthalpy_to_temp` called after shift on extended region

### Task 4: Validation on B26 toolpath
- 400×400×52, AMR, 20 threads, timax=0.092
- See `results.md` for full comparison

---

## Files Modified

| File | Action |
|------|--------|
| `mod_param.f90` | Added `predict_flag` to `&output_control` namelist |
| `mod_predict.f90` | **NEW** — prediction subroutines (separate module for future ML integration) |
| `mod_entot.f90` | Remove prediction code (moved to mod_predict.f90) |
| `main.f90` | Call prediction before iter_loop |
| `input_param.txt` | Added `predict_flag=0` (default off) |
| `compile.sh` | Added `mod_predict.f90` |

## Key Design Decisions

1. **Integer shift, not interpolation**: avoids numerical diffusion at solidification front
2. **All fields shifted together**: maintains consistency between enthalpy and velocity
3. **Local region only**: avoids disturbing ambient far field
4. **Extended region**: covers cells ahead of laser in scan direction
5. **Separate module**: `mod_predict.f90` can be replaced with ML model in future
6. **Flag-controlled**: `predict_flag=0` disables prediction completely (default)
