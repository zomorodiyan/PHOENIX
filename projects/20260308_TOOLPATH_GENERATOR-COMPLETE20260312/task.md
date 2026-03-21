# Task: Toolpath Generator for Rectangle Fill

## Objective
Create `toolpath_generator_rectangle.py` in `fortran_new/ToolFiles/` that generates a laser scan toolpath filling a rectangular region in X-Y. Output format matches `.crs` files (see `B26.crs` for reference).

---

## Output Format (`.crs`)
Space-delimited columns, one waypoint per line:
| Column | Description |
|--------|-------------|
| 1 | Time (s) |
| 2 | X coordinate (m) |
| 3 | Y coordinate (m) |
| 4 | Z coordinate (m) |
| 5 | Laser on/off (1 = on, 0 = off) |

- First line: `0.0  <start_x>  <start_y>  <start_z>  0` (initial position, laser off)
- Laser-on segments: each scan track is a pair of waypoints (start, end) with laser = 1
- Laser-off segments: turn-around between tracks, with laser = 0
- Last line: laser off (dwell after final track)

---

## Input Parameters

| Parameter | Description |
|-----------|-------------|
| `start_x, start_y, start_z` | Starting corner coordinates (m). Z is constant throughout. |
| `size_x` | Rectangle extent in X direction (m) |
| `size_y` | Rectangle extent in Y direction (m) |
| `scan_axis` | Scan along `"x"` or `"y"` |
| `bidirectional` | `True` = alternating direction; `False` = unidirectional (return jump between tracks) |
| `hatch_spacing` | Nominal hatch spacing (m). Adjusted so the last track aligns exactly with the rectangle boundary (see note below). |
| `scan_speed` | Laser scan speed (m/s) |
| `turnaround_time` | Time with laser off between tracks (s) |
| `output_filename` | Name of the output `.crs` file |
| `rotation_angle` | Counter-clockwise rotation angle (degrees) about the rectangle center (see Task 3) |

### Hatch Spacing Adjustment
The number of tracks is `round(size_perp / hatch_spacing) + 1`, where `size_perp` is the rectangle size perpendicular to scan direction. The actual hatch spacing is then `size_perp / (n_tracks - 1)`, ensuring the first and last tracks sit exactly on the rectangle edges.

---

## Task 1 — Core Toolpath Generation (no rotation)
Generate the toolpath for `rotation_angle = 0`:
- Compute adjusted hatch spacing and number of tracks
- For each track, compute start/end coordinates and time
- Handle both directional and bi-directional modes
- Handle both scan_axis = `"x"` and `"y"`
- Write output to `.crs` file matching the reference format

---

## Task 2 — Visualization (PNG output)
After generating the toolpath, produce a PNG plot:
- Draw each scan track as an arrow showing scan direction
- Mark the starting point with a green dot only (no text label)
- Display all input parameters as a text block at the top of the plot (start coordinates, size_x, size_y, scan_axis, bidirectional, hatch_spacing (nominal and actual), scan_speed, turnaround_time, rotation_angle, n_tracks)
- Save as `<output_filename>.png` alongside the `.crs` file

---

## Task 3 — Rotation Support
Add rotation_angle support. The scan direction remains along the **global** X or Y axis before and after rotation — only the rectangle rotates, not the scan lines.

### Geometry
- Rectangle center: `C = (start_x + size_x/2, start_y + size_y/2)` before rotation
- The 4 corners of the rectangle rotate CCW by `rotation_angle` about `C`
- The starting point also rotates with the rectangle
- After rotation, the rectangle becomes a rotated quadrilateral in global X-Y

### Starting Point After Rotation
- The starting point is **not** the rotated original corner. Instead:
  - For `scan_axis = "x"`: starting point = the **bottommost vertex** (min Y) of the rotated rectangle
  - For `scan_axis = "y"`: starting point = the **leftmost vertex** (min X) of the rotated rectangle
- This ensures scanning always begins from the extremal point in the perpendicular direction

### Scan Strategy
- Scan tracks are always axis-aligned (global X or Y), regardless of rotation
- Each track is a horizontal line (if `scan_axis = "x"`) or vertical line (if `scan_axis = "y"`) that intersects the rotated rectangle
- **Track lengths vary**: each track starts and ends at the intersection of the scan line with the rotated rectangle edges
- Perpendicular coverage: scan lines span the full perpendicular extent of the rotated rectangle (its bounding range in Y for x-scan, or X for y-scan)
- Hatch spacing is measured in the perpendicular direction (global Y for x-scan, global X for y-scan)
- Hatch spacing adjustment: `n_tracks` is based on the **perpendicular extent** of the rotated rectangle, adjusted so the first and last tracks touch the outermost vertices

### Implementation
1. Compute the 4 rotated corners of the rectangle
2. Determine the perpendicular extent (min/max Y for x-scan, or min/max X for y-scan)
3. Compute `n_tracks` and adjusted hatch spacing over this extent
4. For each scan line, compute intersections with the 4 edges of the rotated rectangle to get track start/end points
5. Sort/order intersection points to determine scan direction (bidirectional or unidirectional)
6. The first track should be the one closest to the rotated starting point
7. Update the PNG visualization to show the rotated rectangle outline and axis-aligned scan tracks

---

## Task 4 — Testing & Validation
Generate representative test cases and their PNGs to verify correctness:

| Test | size_x | size_y | scan_axis | bidirectional | rotation |
|------|--------|--------|-----------|---------------|----------|
| A | 3 mm | 3 mm | X | yes | 0 deg |
| B | 3 mm | 3 mm | Y | yes | 0 deg |
| C | 3 mm | 3 mm | X | no | 0 deg |
| D | 3 mm | 3 mm | X | yes | 45 deg |
| E | 3 mm | 3 mm | X | yes | 90 deg |
| F | 3 mm | 3 mm | Y | no | 30 deg |
| G | 5 mm | 2 mm | X | yes | 0 deg |
| H | 2 mm | 5 mm | X | yes | 0 deg |
| I | 5 mm | 2 mm | X | yes | 45 deg |
| J | 5 mm | 2 mm | Y | yes | 30 deg |
| K | 2 mm | 5 mm | Y | no | 60 deg |
| L | 3 mm | 3 mm | X | yes | 0 deg | hatch_spacing=0.2mm |
| M | 3 mm | 3 mm | X | yes | 45 deg | hatch_spacing=0.05mm |
| N | 5 mm | 2 mm | X | yes | 30 deg | turnaround_time=0.002s |
| O | 3 mm | 3 mm | Y | yes | 0 deg | hatch_spacing=0.15mm, turnaround_time=0.001s |

Verify:
- Track count and adjusted hatch spacing are correct
- First/last tracks align with rectangle boundaries
- Timing (track duration, turnaround) is consistent
- Rotation transforms all coordinates correctly
- PNGs visually confirm the scan patterns

---

## Task 5 — Integration Test with PHOENIX Simulation
Use the toolpath generator to create a toolpath and run the PHOENIX simulation with it.
No Fortran source code modifications allowed — only replace the toolpath file and adjust `input_param.txt`.

### Setup
- The simulation reads `./ToolFiles/B26.crs` (hardcoded in `mod_toolpath.f90`)
- Back up the original `B26.crs` before replacing
- Generate a new toolpath with parameters matching the simulation domain (4mm x 4mm x 0.7mm)
- Adjust `timax` in `input_param.txt` to match the new toolpath total time

### Test Run
- Generate toolpath: `start=(0.5, 0.5, 0.6975) mm`, `size_x=3mm`, `size_y=3mm`, `scan_axis=x`, bidirectional, `hatch_spacing=0.1mm`, `scan_speed=1.23 m/s`, `turnaround_time=0.5ms`, `rotation=0 deg`
- Replace `B26.crs` with generated toolpath
- Update `timax` to match toolpath end time
- Compile and run simulation
- Verify simulation completes without errors

---

## Notes
