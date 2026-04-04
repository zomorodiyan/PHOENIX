# Thermal-Mechanical — Results

## Single-Track (single_track.crs)

- **Grid**: 200×200×52 thermal, 51×51×13 FEM (mech_mesh_ratio=4)
- **timax**: 0.002 s (100 thermal steps, 10 mechanical solves)
- **Threads**: 4

| Mech Step | Newton res | max von Mises (MPa) | Yield elements |
|-----------|-----------|---------------------|---------------|
| 1 | 3.06E-02 | 30.3 | 5 |
| 5 | 3.39E-02 | 49.4 | 4 |
| 10 | 2.55E-02 | 126.9 | 11 |

| Metric | Value |
|--------|-------|
| Mechanical CPU time | 2764 s (90% of total) |
| FEM memory | 24 MB |
| VTK files | 10 |

## Multi-Track (center_rot0.crs, 11 tracks)

- **Grid**: 200×200×52 thermal, 51×51×13 FEM (mech_mesh_ratio=4)
- **timax**: 0.015 s (750 thermal steps, 75 mechanical solves)
- **Threads**: 4

| Mech Step | Newton res | max von Mises (MPa) | Yield elements |
|-----------|-----------|---------------------|---------------|
| 1 | 3.27E-02 | 27.1 | 4 |
| 25 | 3.46E-02 | 131.1 | 19 |
| 50 | 3.12E-02 | 133.3 | 48 |
| 75 | 1.70E-02 | 132.4 | 54 |

| Metric | Value |
|--------|-------|
| Total wall time | 1248 s (20.8 min) |
| Mechanical CPU | 1197 s (24.6% of total) |
| Thermal CPU | 3669 s (75.4% of total) |
| FEM memory | 24 MB |
| VTK files | 15 mech + 16 thermal |
| Mech history PNG | Generated |

### Observations

1. **Stress evolution**: von Mises stress grows from 27 MPa to ~132 MPa over 11 tracks, plateauing below the 250 MPa yield stress
2. **Yield elements**: growing from 4 to 54 as more material solidifies and accumulates thermal strain
3. **Newton convergence**: all 75 solves converged with relative residual ~0.02-0.03
4. **Threading**: 4 threads optimal for the 51×51×13 FEM grid. 20 threads caused massive OpenMP overhead (>16 min for a single solve vs ~10s with 4 threads)
5. **mech_mesh_ratio=4**: reduces FEM from 200³ to 51³ nodes, making the solver tractable while maintaining sufficient resolution for stress fields

## Bug Fixes During Implementation

1. **Tangent/residual inconsistency**: ebe_matvec used reference Ke, residual used per-element B → fixed with precomputed per-Z-layer Ke
2. **Array size mismatch**: VTK/history passed thermal-sized temp array to FEM-sized routines → fixed with `T_fem_last`
3. **FEM grid sizing**: memory report showed thermal grid instead of FEM grid → fixed with public Nnx/Nny/Nnz
4. **Per-element Ke in CG**: computing 24×24 Ke every CG iteration was too slow → precomputed per-Z-layer (only dz varies)
