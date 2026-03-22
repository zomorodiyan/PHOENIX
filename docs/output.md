# Output & Visualization

All output files are written to `fortran_new/result/<case_name>/`.

## Output Files

| File | Description | When |
|------|-------------|------|
| `<case>_output.txt` | Iteration-by-iteration text log | Every time step |
| `<case>_vtkmov{N}.vtk` | Binary VTK structured grid | Every `outputintervel` steps |
| `<case>_timing_report.txt` | CPU time breakdown by module | End of simulation |
| `<case>_memory_report.txt` | Peak memory usage (VmHWM, VmRSS) | End of simulation |
| `<case>_defect_report.txt` | Defect fraction summary | End of simulation |
| `<case>_maxtemp.vtk` | Max temperature field (VTK) | End of simulation |
| `<case>_defect.vtk` | Defect classification field (VTK) | End of simulation |
| `<case>_thermal_history.txt` | Temperature at 10 monitoring points | Every time step |
| `<case>_thermal_history.png` | Temperature evolution plot | End of simulation |
| `<case>_meltpool_history.txt` | Melt pool length, depth, width, volume, Tpeak | Every global step |
| `<case>_meltpool_history.png` | Melt pool geometry evolution plot | End of simulation |
| `<case>_micro_report.txt` | Microstructure statistics (G, R, PDAS, SDAS) | End of simulation (micro_flag=1) |
| `<case>_microstructure.vtk` | Microstructure fields (VTK, multi-scalar) | End of simulation (micro_flag=1) |
| `<case>_crack_report.txt` | Crack risk statistics (CSI, strain rate) | End of simulation (crack_flag=1) |
| `<case>_crack_risk.vtk` | Crack risk fields (VTK, multi-scalar) | End of simulation (crack_flag=1) |

## output.txt Format

Each time step produces a 4-line block:

```
  time  iter  time/iter  tot_iter  res_enth  res_mass  res_u  res_v  res_w  [res_spec]
 2.10E-05   60    0.006        60   1.2E-04   5.0E-02  1.3E-02  ...

  Tmax        umax       vmax         wmax       length       depth     width
  1594.47   0.00E+00    0.00E+00    0.00E+00   ...

  north  south  top  toploss  bottom  west  east  hout  accu  hin  heatvol  ratio
    0.0    0.0    0.0    0.0    0.0    0.0    0.0   0.0   0.0   0.0   0.0   0.00

progress%  beam_posx  beam_posy  beam_posz  power  scanspeed  speedx  speedy
     0.70  5.258E-04  2.000E-03  7.000E-04  300.0    1.230    1.230    0.000
```

| Field | Description |
|-------|-------------|
| `res_enth` | Enthalpy equation residual |
| `res_mass` | Mass conservation (pressure correction) residual |
| `res_u/v/w` | Momentum residuals |
| `res_spec` | Species residual (only when `species_flag=1`) |
| `Tmax` | Peak temperature in domain (K) |
| `umax/vmax/wmax` | Peak velocities (m/s) |
| `length/depth/width` | Melt pool dimensions (m) |
| `ratio` | Energy balance ratio (should be ~1.0) |
| `progress%` | Simulation progress (%) |

## VTK Files

Binary structured-grid VTK files, viewable in ParaView.

### Scalar Fields

| Field Name | Description | Unit | Always |
|------------|-------------|------|--------|
| `T` | Temperature | K | Yes |
| `vis` | Viscosity | Pa*s | Yes |
| `diff` | Thermal diffusivity | m^2/s | Yes |
| `den` | Density | kg/m^3 | Yes |
| `solidID` | Solidification track ID | - | Yes |
| `local` | Local solver region marker | - | Yes |
| `fracl` | Liquid fraction (0=solid, 1=liquid) | - | Yes |
| `concentration` | Species mass fraction (C=1: primary) | - | species_flag=1 |
| `tsolid_field` | Composition-dependent solidus | K | species_flag=1 |

### Vector Fields

| Field Name | Description | Unit |
|------------|-------------|------|
| `Velocity` | Flow velocity (interpolated to cell centers) | m/s |

### Opening in ParaView

1. File → Open → select `<case>_vtkmov*.vtk`
2. Click "Apply"
3. Select scalar/vector field from dropdown
4. Use animation controls to step through time

!!! tip "Useful ParaView Filters"
    - **Threshold**: Show only melt pool (`T > 1563`)
    - **Slice**: Cut through domain to see cross-sections
    - **Glyph**: Visualize velocity vectors
    - **Calculator**: Compute derived quantities (e.g., `T - 1563` for superheat)

## Timing Report

```
============================================
  PHOENIX Module Timing Report
============================================
  Total iterations (itertot): 17055
  Total CPU time:     2513.425 s
  Total wall time:      765.089 s
--------------------------------------------
  Module              |   Time(s)   |   Ratio(%)
--------------------------------------------
             mod_prop |     440.355  |   17.52%
             mod_sour |     406.017  |   16.15%
             ...
           mod_species|       4.355  |    0.17%
```

## Memory Report

```
============================================
  PHOENIX Memory Report
============================================
  VmPeak:  ...   (peak virtual memory)
  VmHWM:   ...   (peak physical RAM — use this to size your machine)
  VmRSS:   ...   (current physical RAM at report time)
  VmData:  ...   (heap + stack, excludes shared libraries)
```

## Defect Report

```
  === Defect Analysis (maxtemp_determ) ===
  Defect fraction:           X.XXXXXX %
  Lack-of-fusion fraction:   X.XXXXXX %
  Keyhole porosity fraction: X.XXXXXX %
```

Spatial distribution available in `<case>_defect.vtk` and `<case>_maxtemp.vtk`.

## Microstructure Report (micro_flag=1)

Computed at the moment each cell solidifies (liquid fraction drops to 0). Only updated during global solver timesteps.

```
  === Microstructure Analysis ===
  Solidified cells: 16013
  G (K/m):   min= 1.15E+06  max= 2.93E+07
  R (m/s):   min= 3.51E-03  max= 1.16E+00
  PDAS (m):  min= 2.01E-08  max= 8.06E-08
  SDAS (m):  min= 6.65E-08  max= 3.39E-07
```

VTK file `<case>_microstructure.vtk` contains 5 scalar fields:

| Field | Description | Unit |
|-------|-------------|------|
| `cooling_rate` | Cooling rate at solidification | K/s |
| `thermal_gradient` | Temperature gradient magnitude G | K/m |
| `solidification_rate` | Solidification front velocity R | m/s |
| `PDAS` | Primary dendrite arm spacing | m |
| `SDAS` | Secondary dendrite arm spacing | m |

## Crack Risk Report (crack_flag=1)

Tracks thermal strain accumulated while material passes through the Brittle Temperature Range (BTR). Only updated during global solver timesteps.

```
  === Crack Risk Analysis ===
  Max CSI:            6.20E-02
  Mean CSI:           2.26E-03
  Max cooling rate:   3.44E+06 K/s
  High-risk frac:     1.842253 %
```

VTK file `<case>_crack_risk.vtk` contains 4 scalar fields:

| Field | Description | Unit |
|-------|-------------|------|
| `crack_csi` | Crack susceptibility index (accumulated thermal strain in BTR) | - |
| `cooling_rate_solid` | Cooling rate at solidification | K/s |
| `strain_rate_solid` | Thermal strain rate at solidification (β × dT/dt) | 1/s |
| `btr_time` | Cumulative time spent in BTR | s |

## Melt Pool History

Time-series log of melt pool geometry, recorded every global solver timestep.

| Column | Description | Unit |
|--------|-------------|------|
| time | Simulation time | s |
| length | Melt pool length (along scan) | m |
| depth | Melt pool depth (from surface) | m |
| width | Melt pool width (transverse) | m |
| volume | Melt pool volume (from liquid fraction) | m³ |
| Tpeak | Peak temperature in domain | K |
| laser_on | Laser state (1=on, 0=off) | - |

Auto-generated plot: `<case>_meltpool_history.png` (4-panel: length, depth/width, volume, Tpeak).
