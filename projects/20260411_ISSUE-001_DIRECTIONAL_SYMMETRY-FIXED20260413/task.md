# Directional Symmetry Regression Check

## Objective

Create a small repeatable regression check for ISSUE-001 that runs four short simulations in the positive and negative x and y directions with the same scan speed, then verifies that each positive/negative pair has matching iteration behavior within a small tolerance.

## Status

COMPLETED

## Scope

- Generate four single-track CRS inputs: +x, -x, +y, -y
- Run short simulations with a fixed small timestep window
- Parse the final total iteration count from each run
- Fail if either directional pair drifts beyond tolerance
- Document the run and results for future reruns

## Planned Checks

1. Build or reuse `cluster_main`.
2. Run four 4-step cases with the same `scan_speed` and `delt`.
3. Compare total iterations for +x vs -x and +y vs -y.
4. Record the pass/fail result in `results.md`.

## Files

- `run_directional_symmetry.sh`: regression runner
- `scan_plus_x.crs`, `scan_minus_x.crs`, `scan_plus_y.crs`, `scan_minus_y.crs`: short single-track inputs
- `log.md`: timestamped execution log
- `results.md`: final comparison summary
