# Issue Template

Use this template for all new issues to keep reporting, debugging, and validation consistent.

---

## ISSUE-XXX — Short, specific title

Status: Open | In Progress | Fixed | Verified | Closed
Severity: High | Medium | Low
Discovered: YYYY-MM-DD
Fixed: YYYY-MM-DD (if applicable)
Verified: YYYY-MM-DD (if applicable)
Owner: Zomorodiyan
Files: `path/to/file1`, `path/to/file2`

### Physical/Simulation Impact

State whether this can change physical conclusions, numerical stability, convergence behavior, or only code clarity/maintainability.

### Reproduction

- Case/setup:
- Input parameter deltas from baseline:
- Command(s) used:
- Observed behavior:
- Expected behavior:

### Root Cause

One short paragraph describing the exact mechanism and where it occurs.

### Fix

Describe what changed and why this is the minimal correct fix.

### Verification

- Validation method:
- Before/after metric(s):
- Acceptance criteria:
- Result:

### Regression Protection

Document a repeatable check and how to run it:

```bash
cd issues/ISSUE-XXX_SHORT_NAME
bash run_regression.sh
```

Script skeleton (`run_regression.sh`):

```bash
#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
FORTRAN_DIR="$ROOT_DIR/fortran_new"

# Compile if needed
if [[ ! -x "$FORTRAN_DIR/cluster_main" ]]; then
  (cd "$FORTRAN_DIR" && bash compile.sh)
fi

# Run case(s)
# (cd "$FORTRAN_DIR" && OMP_NUM_THREADS=4 ./cluster_main > output.txt 2>&1)

# Check result
# Fail loudly if regression detected:
# echo "FAIL: ..." >&2; exit 1
echo "PASS"
```

### Links

- Issue folder: `issues/ISSUE-XXX_SHORT_NAME/`
- Related project: `projects/YYYYMMDD_PROJECT_NAME/`
