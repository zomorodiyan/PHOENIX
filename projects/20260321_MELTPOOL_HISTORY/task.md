# Task: Melt Pool Geometry History Logger

## Objective
Create a time-series log of melt pool dimensions (length, depth, width, volume, aspect ratios) in a clean CSV-like format, with an auto-generated Python plotting script. Follows the same pattern as `thermal_history.txt`.

### Motivation
- Melt pool dimensions are already computed per timestep in `mod_dimen.f90` (`pool_size`) but buried in verbose `output.txt`
- No easy way to extract and plot melt pool geometry vs. time
- Aspect ratios (depth/width) are key indicators of keyhole vs. conduction mode melting
- Melt pool volume is not tracked at all

---

## Design

### Output file: `<case>_meltpool_history.txt`
Header + one row per timestep:
```
# time(s)  length(m)  depth(m)  width(m)  volume(m3)  D/W  Tpeak(K)  laser_power
```

### Melt pool volume
Sum of `fracl(i,j,k) * volume(i,j,k)` over all cells where `fracl > 0`. Computed in the logger, not in `pool_size`.

### Aspect ratio
`D/W = depth / width` (0 when width = 0). Keyhole mode: D/W > 0.5. Conduction mode: D/W < 0.5.

### Python plotting script
Auto-generated `<case>_plot_meltpool.py` that plots:
- Length, depth, width vs. time (3 subplots)
- D/W aspect ratio vs. time
- Melt pool volume vs. time

### Implementation
Add 3 subroutines to `mod_print.f90` (where thermal_history already lives):
- `init_meltpool_history()` — open file, write header, generate Python script
- `write_meltpool_history(time)` — write one row per timestep
- `finalize_meltpool_history()` — close file

No new module needed. No new flag — always enabled (like thermal_history).

---

## Tasks

### Task 1 — Add subroutines to mod_print.f90
- `init_meltpool_history()`: open unit, write header
- `write_meltpool_history(time)`: compute volume from fracl, write row with time, alen, depth, width, volume, D/W, tpeak, laser power
- `finalize_meltpool_history()`: close file, generate plot script

### Task 2 — Integrate into main.f90
- Call `init_meltpool_history()` after `init_thermal_history`
- Call `write_meltpool_history(timet)` after `write_thermal_history(timet)`
- Call `finalize_meltpool_history()` after `finalize_thermal_history`

### Task 3 — Test and validate
- Run simulation, verify meltpool_history.txt is generated
- Verify values match output.txt (alen, depth, width should match)
- Run the plot script

---

## Notes
- No new flag, no new module — minimal changes
- Follows exact pattern of thermal_history (init/write/finalize)
- Volume computation uses existing `fracl` and `volume` arrays
