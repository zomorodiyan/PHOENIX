# Test Report: Toolpath Generator for Rectangle Fill

**Date:** 2026-03-08
**Script:** `fortran_new/ToolFiles/toolpath_generator_rectangle.py`

---

## Test Summary

| Test | size_x | size_y | scan | dir | rot | hatch | turnaround | n_tracks | waypoints | Status |
|------|--------|--------|------|-----|-----|-------|------------|----------|-----------|--------|
| A | 3 mm | 3 mm | X | bidir | 0° | 0.10 mm | 0.50 ms | 31 | 64 | PASS |
| B | 3 mm | 3 mm | Y | bidir | 0° | 0.10 mm | 0.50 ms | 31 | 64 | PASS |
| C | 3 mm | 3 mm | X | unidir | 0° | 0.10 mm | 0.50 ms | 31 | 64 | PASS |
| D | 3 mm | 3 mm | X | bidir | 45° | 0.10 mm (actual 0.1010) | 0.50 ms | 43 | 84 | PASS |
| E | 3 mm | 3 mm | X | bidir | 90° | 0.10 mm | 0.50 ms | 31 | 64 | PASS |
| F | 3 mm | 3 mm | Y | unidir | 30° | 0.10 mm | 0.50 ms | 42 | 82 | PASS |
| G | 5 mm | 2 mm | X | bidir | 0° | 0.10 mm | 0.50 ms | 21 | 44 | PASS |
| H | 2 mm | 5 mm | X | bidir | 0° | 0.10 mm | 0.50 ms | 51 | 104 | PASS |
| I | 5 mm | 2 mm | X | bidir | 45° | 0.10 mm (actual 0.1010) | 0.50 ms | 50 | 98 | PASS |
| J | 5 mm | 2 mm | Y | bidir | 30° | 0.10 mm (actual 0.1006) | 0.50 ms | 54 | 106 | PASS |
| K | 2 mm | 5 mm | Y | unidir | 60° | 0.10 mm (actual 0.1006) | 0.50 ms | 54 | 106 | PASS |
| L | 3 mm | 3 mm | X | bidir | 0° | 0.20 mm | 0.50 ms | 16 | 34 | PASS |
| M | 3 mm | 3 mm | X | bidir | 45° | 0.05 mm (actual 0.0499) | 0.50 ms | 86 | 170 | PASS |
| N | 5 mm | 2 mm | X | bidir | 30° | 0.10 mm (actual 0.1008) | 2.00 ms | 43 | 84 | PASS |
| O | 3 mm | 3 mm | Y | bidir | 0° | 0.15 mm | 1.00 ms | 21 | 44 | PASS |

**All common parameters:** start=(0.50, 0.50, 0.70) mm, scan_speed=1.20 m/s (unless noted)

---

## Verification Checks

### 1. Hatch Spacing Adjustment
- **Square, no rotation (A):** 3mm / (31-1) = 0.1000 mm — exact match
- **Square, 45° rotation (D):** perpendicular extent = 3*sqrt(2) = 4.243 mm, 4.243/(43-1) = 0.1010 mm — correct
- **90° rotation (E):** same as 0° (square symmetry), 31 tracks — correct
- **Wider hatch (L):** 3mm / (16-1) = 0.2000 mm — exact match
- **Finer hatch (M):** 4.243/(86-1) = 0.0499 mm — correct

### 2. Track Count Consistency
- Non-rotated rectangle: n_tracks = round(size_perp/hatch) + 1
  - 3mm/0.1mm + 1 = 31 (A, B, C, E) — correct
  - 2mm/0.1mm + 1 = 21 (G) — correct
  - 5mm/0.1mm + 1 = 51 (H) — correct
- Rotated: perpendicular extent increases, more tracks needed — verified in D, F, I, J, K

### 3. Starting Point
- x-scan: starting point = bottommost vertex (min Y of rotated rectangle)
- y-scan: starting point = leftmost vertex (min X of rotated rectangle)
- Verified in all rotated cases (D, E, F, I, J, K, M, N): green dot at extremal vertex

### 4. Scan Direction
- Bidirectional (A, D): arrows alternate left/right — verified in PNGs
- Unidirectional (C): all arrows point same direction — verified in PNG
- Y-scan (B, O): vertical arrows — verified in PNGs

### 5. Rotation Geometry
- Scan lines remain axis-aligned after rotation (global X or Y)
- Track lengths vary: shorter near tips, longest at widest point — verified in D, I, N
- Rectangle outline correctly shows rotated quadrilateral

### 6. Different Hatch Spacing
- L (0.2mm): 16 tracks, wider spacing visible in PNG
- M (0.05mm): 86 tracks, dense fill visible in PNG
- O (0.15mm): 21 tracks — correct

### 7. Different Turnaround Time
- N (2.0ms): same geometry as similar cases, longer turnaround reflected in .crs timing
- O (1.0ms): turnaround time correctly applied

---

## Task 5: PHOENIX Integration Test

- Generated toolpath: start=(0.5, 0.5, 0.6975) mm, size_x=3mm, size_y=3mm, scan_axis=x, bidir, hatch=0.1mm, speed=1.23 m/s, turnaround=0.5ms
- Replaced `B26.crs` with generated toolpath (original backed up as `B26_original.crs`)
- Toolpath total time: 0.0911s (within timax=0.092s)

### Result: **PASS**
- Simulation completed: 4601 steps, 2 hr 13 min 49 s wall time
- Laser scanning range correctly detected: X=[0.5, 3.5] mm, Y=[0.5, 3.5] mm — matches generated toolpath
- Defect analysis completed:
  - Defect fraction: 7.91%
  - Lack-of-fusion: 5.69%
  - Keyhole porosity: 2.22%
- 230 VTK output files generated
- No errors — floating-point warnings (IEEE_UNDERFLOW, IEEE_DENORMAL) are normal for CFD

---

## Output Files
- `.crs` files: `fortran_new/ToolFiles/test_output/test_{A-O}.crs`
- `.png` files: `fortran_new/ToolFiles/test_output/test_{A-O}.png`
- Integration test toolpath: `fortran_new/ToolFiles/B26.crs` (generated)
- Original backup: `fortran_new/ToolFiles/B26_original.crs`
