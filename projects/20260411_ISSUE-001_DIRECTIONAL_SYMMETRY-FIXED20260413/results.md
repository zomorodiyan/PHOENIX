# Directional Symmetry Regression Check - Results

## Run Command

`bash ./run_directional_symmetry.sh 4`

## Final Comparison

| Case | Direction | Total iterations |
|------|-----------|------------------|
| `reg_plus_x` | +x | 90 |
| `reg_minus_x` | -x | 89 |
| `reg_plus_y` | +y | 91 |
| `reg_minus_y` | -y | 92 |

## Directional Checks

- +x / -x difference: 1 iteration
- +y / -y difference: 1 iteration
- Tolerance: 2 iterations
- Result: PASS

## Conclusion

The 4-case regression harness completed successfully on 4 timesteps per case. Both directional pairs stayed within tolerance, so the ISSUE-001 fix is protected by a repeatable check.

