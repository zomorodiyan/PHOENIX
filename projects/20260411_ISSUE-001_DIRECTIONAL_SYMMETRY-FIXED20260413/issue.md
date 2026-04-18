## ISSUE-001 — ishift_field incorrect when dj > 0

Status: Verified
Severity: Medium
Discovered: 2026-04-11
Fixed: 2026-04-13
Verified: 2026-04-13
Owner: Zomorodiyan
Files: `fortran_new/thermal_fluid_solver/mod_predict.f90`

### Physical/Simulation Impact

Prediction quality became direction-dependent when `dj > 0`, degrading solver convergence for +y scanning. This primarily affects numerical convergence efficiency rather than final physics for current production runs where `scanvely = 0` keeps `dj = 0`.

### Reproduction

- Case/setup: two single-track toolpaths, one +y and one -y.
- Input parameter deltas from baseline: `delt = 2x10^-5 s`, `timax = 8x10^-5 s`, 400x400 domain, AMR enabled.
- Command(s) used: standard compile/run workflow for each case (`scan_plus_y`, `scan_minus_y`, `fixed_plus_y`, `fixed_minus_y`).
- Observed behavior: buggy code required substantially more inner iterations for +y than -y.
- Expected behavior: comparable convergence behavior for symmetric opposite scan directions.

### Root Cause

The j-loop in `ishift_field` always iterated forward. When `dj > 0`, source cell `j - dj` had already been overwritten before being read, causing a cascading copy artifact. The i-loop already handled direction by switching traversal based on the sign of `di`; j-direction logic was missing the same guard.

### Fix

Mirror i-direction traversal logic for j-direction traversal:

```fortran
! Before (incorrect)
do j = js, je

! After (fixed)
if (dj >= 0) then
    do j = je, js, -1
else
    do j = js, je
endif
```

### Verification

- Validation method: compare inner iteration counts per timestep for +y and -y scans before/after fix.
- Before/after metric(s):

| Run | Scan direction | Code | Total iterations | Avg per timestep |
|-----|----------------|------|-----------------|-----------------|
| `scan_plus_y`  | +y | buggy  | 115 | 28.8 |
| `scan_minus_y` | −y | buggy  |  92 | 23.0 |
| `fixed_plus_y` | +y | fixed  |  91 | 22.8 |
| `fixed_minus_y`| −y | fixed  |  92 | 23.0 |

- Acceptance criteria: remove directional convergence bias.
- Result: buggy code showed ~25% more iterations for +y (28.8 vs 23.0 avg/timestep); fixed code restored near-symmetry (22.8 vs 23.0 avg/timestep).

### Regression Protection

```bash
cd projects/20260411_ISSUE-001_DIRECTIONAL_SYMMETRY-FIXED20260413
bash run_regression.sh
```

Runs +x/-x and +y/-y scan pairs and fails if total iteration counts differ by more than 2.

### Links

- Issue folder: `projects/20260411_ISSUE-001_DIRECTIONAL_SYMMETRY-FIXED20260413/`
- Task log: [task.md](task.md)
- Results: [results.md](results.md)
- Work log: [log.md](log.md)
