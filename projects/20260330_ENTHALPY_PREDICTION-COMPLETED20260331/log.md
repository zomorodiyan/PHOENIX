# Enthalpy Prediction — Execution Log

### [2026-03-30 21:46:00] Implementation
- Added `predict_enthalpy_shift(vx, vy, dt)` to `mod_entot.f90`
- Bilinear interpolation from upstream position `(x - vx*dt, y - vy*dt)`
- Called from `main.f90` before iter_loop, only during laser-on steps
- Compilation: clean build

### [2026-03-30 21:48:31] First attempt: enthalpy-only shift, 20 threads
- Very slow — 20 threads on small grid (170×60×52) has too much OpenMP overhead
- Only 15 thermal history records in 2+ minutes
- Killed and switched to 4 threads

### [2026-03-30 21:52:00] Second attempt: all-field shift (enthalpy + u,v,w,p)
- Added `predict_field_shift` for generic 3D fields
- Result: WORSE — shifted velocity creates artificial divergence
- 52 iter at t=1.78ms vs baseline 36 iter
- Conclusion: do NOT shift velocity/pressure

### [2026-03-30 21:55:00] Third attempt: enthalpy-only shift, 4 threads
- Reverted to enthalpy-only shift
- Completed successfully

### [2026-03-30 22:00:00] Results (4 threads — UNFAIR comparison)
- Iterations: 36,001 → 35,583 (-1.2%)
- Wall time: 387.5s → 311.3s (-19.6%) — appeared faster but baseline was 20 threads

### [2026-03-31 00:00:00] Re-run with 20 threads (fair comparison)
- Baseline: 387.5s wall, 7594s CPU (20 threads)
- Prediction: 414.9s wall, 8089s CPU (20 threads)
- Result: **7% SLOWER** — shift creates enthalpy/velocity inconsistency, heavier source terms

### [2026-03-31 00:10:00] Code reverted
- Removed predict_enthalpy_shift and predict_field_shift from mod_entot.f90
- Removed prediction call from main.f90
- Compilation: clean build
- Project marked as REJECTED
