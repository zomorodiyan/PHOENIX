## ISSUE-002 — Top surface convective loss missing; radiation double-counted

Status: Open
Severity: High
Discovered: 2026-04-13
Fixed: N/A
Verified: N/A
Owner: Zomorodiyan
Files: `fortran_new/thermal_fluid_solver/mod_bound.f90`

### Physical/Simulation Impact

Top-surface heat-loss physics is incorrect: convection is dropped and radiation is double-counted, which can materially bias thermal fields, melt-pool evolution, and cooling rates.

### Reproduction

- Case/setup: any case using top boundary heat-loss evaluation in `bound_enthalpy`.
- Input parameter deltas from baseline: none required.
- Command(s) used: code inspection of top-surface boundary update.
- Observed behavior: top update uses `hlossradia` twice.
- Expected behavior: use `hlossradia + hlossconvec`.

### Root Cause

In top-surface enthalpy and accumulated top-loss updates, the second loss term was mistakenly copied as `hlossradia` instead of `hlossconvec`, even though `hlossconvec` is calculated.

### Fix

Replace the second `hlossradia` with `hlossconvec` in both top-surface lines:

```fortran
! Lines 92-93 — current (incorrect)
enthalpy(i,j,nk) = enthalpy(i,j,nkm1) + (heatin(i,j) - hlossradia - hlossradia)/ctmp1
ahtoploss = ahtoploss + (hlossradia + hlossradia)*areaij(i,j)

! Fixed
enthalpy(i,j,nk) = enthalpy(i,j,nkm1) + (heatin(i,j) - hlossradia - hlossconvec)/ctmp1
ahtoploss = ahtoploss + (hlossradia + hlossconvec)*areaij(i,j)
```

### Verification

- Validation method: compare boundary-loss decomposition and integrated top losses before/after patch on a short controlled run.
- Before/after metric(s): `ahtoploss` components (radiative vs convective), top temperature history.
- Acceptance criteria: non-zero convective contribution at top when `temp > tempAmb`, and no duplicated radiative term.
- Result: Pending.

### Regression Protection

Pending — add a boundary-loss sanity check that validates each face uses the intended `(radiation + convection)` structure.

### Links

- Issue folder: `projects/20260413_ISSUE-002_TOP_SURFACE_HEAT_LOSS/`
