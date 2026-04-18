## ISSUE-003 — Bottom surface `hlossconvec` silently includes radiation

Status: Open
Severity: Medium
Discovered: 2026-04-13
Fixed: N/A
Verified: N/A
Owner: Zomorodiyan
Files: `fortran_new/thermal_fluid_solver/mod_bound.f90`

### Physical/Simulation Impact

Primarily a correctness/maintainability risk due to inconsistent variable semantics across boundaries. Can hide future physics bugs and makes audits of energy-loss terms error-prone.

### Reproduction

- Case/setup: any case; issue visible by code inspection.
- Input parameter deltas from baseline: none required.
- Command(s) used: inspect bottom boundary-loss expression in `bound_enthalpy`.
- Observed behavior: `hlossconvec` at bottom includes both convection and radiation.
- Expected behavior: convection and radiation as separate named terms, consistent with other faces.

### Root Cause

Bottom-face implementation combined two physical mechanisms into one variable while retaining a convection-only name, diverging from the naming/term separation used elsewhere.

```fortran
! Line 100 — current (misleading)
hlossconvec = htck1*(temp(i,j,1)-tempAmb) + emiss*sigm*(temp(i,j,1)**4-tempAmb**4)
```

### Fix

Pending. Refactor bottom-face loss expression so variable naming matches physical content — either split into `hlossconvec` + `hlossradia`, or rename combined variable explicitly.

### Verification

- Validation method: ensure symbolic consistency of boundary terms and compare integrated bottom loss against separated-term implementation.
- Before/after metric(s): bottom-face radiative/convective reporting consistency; energy-balance diagnostics.
- Acceptance criteria: variable names map one-to-one with physical mechanism.
- Result: Pending.

### Regression Protection

Pending — add a boundary-term consistency audit: convection variables must not include radiative terms unless explicitly labeled as combined losses.

### Links

- Issue folder: `projects/20260413_ISSUE-003_BOTTOM_SURFACE_HLOSSCONVEC/`
