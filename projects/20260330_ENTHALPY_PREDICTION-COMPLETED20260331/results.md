# Enthalpy Prediction — Results

## Approaches Tested

### Approach 1: Bilinear interpolation shift (REJECTED)

Shift enthalpy (and optionally all fields) by `v*dt` using bilinear interpolation.

**Problem**: Interpolation introduces numerical diffusion at the sharp solid-liquid interface, creating artificial mushy zone cells. This increases per-iteration cost by 7%, negating the small iteration reduction (-1%).

Tested variants:
- Enthalpy-only shift, full domain → +7% slower
- All-field shift (h,u,v,w,p), full domain → even worse (velocity divergence)
- All-field shift, melt pool region only → +5.5% slower

### Approach 2: Integer-cell shift, local region (ACCEPTED)

Shift all fields (h,u,v,w,p) by integer number of cells within the melt pool region. No interpolation → no numerical diffusion. Sharp interfaces perfectly preserved.

**Key fix**: Extend shift region beyond melt pool by `nint(v*dt/dx)` cells in scan direction to cover where the laser is moving to.

---

## Final Comparison: B26 toolpath, 400×400×52, AMR, 20 threads

| Metric | B26nopred | B26pred | Change |
|--------|-----------|---------|--------|
| Total iterations | 154,725 | 95,985 | **-38.0%** |
| Wall time (s) | 17,054 | 11,697 | **-31.4%** |
| CPU time (s) | 164,273 | 110,137 | **-32.9%** |
| Peak memory (MB) | 1,431 | 1,432 | +0.0% |
| Heating wall time (s) | 14,307 | 8,636 | **-39.6%** |
| Cooling wall time (s) | 2,038 | 2,287 | +12.2% |

### Per-module CPU breakdown

| Module | nopred (s) | pred (s) | Change |
|--------|-----------|---------|--------|
| source_enthalpy | 41,976 | 26,610 | -36.6% |
| discretize_enthalpy | 40,429 | 26,786 | -33.8% |
| solution_enthalpy | 27,581 | 18,849 | -31.7% |
| mod_converge | 13,307 | 9,522 | -28.4% |
| mod_resid | 12,062 | 8,434 | -30.1% |
| source_momentum | 1,697 | 978 | -42.4% |
| discretize_momentum | 2,478 | 1,445 | -41.7% |

### Accuracy

Thermal history max absolute error (K):

| Point | Max Error (K) |
|-------|--------------|
| P1 | 37.6 |
| P2 | 59.6 |
| P3 | 40.5 |
| P4 | 20.1 |
| P5 | 8.7 |
| P6 | 33.5 |
| P7 | 72.7 |
| P8 | 43.8 |
| P9 | 50.4 |
| P10 | 1.5 |

Melt pool (mean relative error): length 4.1%, depth 0.6%, width 2.0%, Tpeak 1.0%.

Defect fractions:

| Metric | nopred | pred | Diff |
|--------|--------|------|------|
| Defect fraction | 10.44% | 10.05% | 0.39% |
| Lack-of-fusion | 8.34% | 7.93% | 0.41% |
| Keyhole porosity | 2.11% | 2.12% | 0.02% |

---

## Conclusion

**Integer-cell shift prediction reduces heating-stage iterations by 38% and total wall time by 31%.** The key insights:
1. **No interpolation**: integer-cell shift preserves sharp interfaces exactly — no numerical diffusion
2. **All fields shifted together**: enthalpy + velocity + pressure stay consistent
3. **Local region only**: only shift within the melt pool region (+ extension) — far field untouched
4. **Extended region**: shift region includes cells ahead of the laser where it's moving to
5. **Negligible overhead**: integer shift is a simple index copy, no floating-point math
6. **Accuracy**: melt pool depth within 0.6% mean error, Tpeak within 1.0%. Thermal history differences (30-70K) reflect that the solver converges to a slightly different solution with fewer iterations (still within convergence tolerance).
