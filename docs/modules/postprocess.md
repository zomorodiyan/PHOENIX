# Post-Processing Modules

Defect prediction, output, timing.

---

## mod_defect.f90 — `defect_field`

Post-simulation defect detection from peak temperature history. Supports uniform defect mesh when `adaptive_flag==1`.

- `update_max_temp()` — called every timestep, tracks maximum temperature per cell
- `compute_defect_determ()` — classifies cells: lack-of-fusion ($T_{max} < T_s$), sound, keyhole ($T_{max} > T_b$)
- `write_defect_report()` — volumetric defect fractions, VTK output (`maxtemp.vtk`, `defect.vtk`)
- `compute_scan_range()` / `point_in_scan_region()` — builds convex hull of scanned region from toolpath

---

## mod_print.f90 — `printing`

All output routines:

- `outputres()` — per-timestep text log (residuals, pool size, energy balance)
- `Cust_Out()` — VTK snapshots at `outputintervel` frequency (T, velocity, vis, diff, den, solidID, fracl, [concentration, tsolid_field])
- `init/write/finalize_thermal_history()` — temperature at 10 monitoring points + Python plot
- `update_thermal_history_indices()` — recomputes monitoring point grid indices after AMR remesh
- `init/write/finalize_meltpool_history()` — melt pool geometry time-series + Python plot

---

## mod_timing.f90 — `timing`

- `write_timing_report()` — CPU time breakdown by module (17 categories)
- `write_memory_report()` — peak memory from `/proc/self/status`
