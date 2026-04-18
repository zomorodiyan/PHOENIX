#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
FORTRAN_DIR="$ROOT_DIR/fortran_new"
BASE_INPUT_FILE="$FORTRAN_DIR/inputfile/input_param.txt"
PROJECT_INPUT_TEMPLATE="$SCRIPT_DIR/input_param.issue001.template.txt"
PROJECT_INPUT_WORK="$SCRIPT_DIR/input_param.issue001.work.txt"
RESULT_DIR="$SCRIPT_DIR/run_outputs"
SOLVER_RESULT_DIR="$FORTRAN_DIR/result"
OMP_THREADS="${1:-4}"
TIMAX="8.0e-5"
DELTA_T="2.0e-5"
TOLERANCE="2"

mkdir -p "$RESULT_DIR"

if [[ ! -f "$PROJECT_INPUT_TEMPLATE" ]]; then
  echo "Missing project input template: $PROJECT_INPUT_TEMPLATE" >&2
  exit 1
fi

if [[ ! -x "$FORTRAN_DIR/cluster_main" ]]; then
  (cd "$FORTRAN_DIR" && bash compile.sh)
fi

backup_file="$(mktemp)"
cp "$BASE_INPUT_FILE" "$backup_file"
restore_input() {
  cp "$backup_file" "$BASE_INPUT_FILE"
  rm -f "$PROJECT_INPUT_WORK"
  rm -f "$backup_file"
}
trap restore_input EXIT

set_input() {
  local case_name="$1"
  local toolpath_file="$2"
  python3 - "$PROJECT_INPUT_WORK" "$case_name" "$toolpath_file" "$TIMAX" "$DELTA_T" <<'PY'
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
case_name = sys.argv[2]
toolpath_file = sys.argv[3]
timax = sys.argv[4]
delt = sys.argv[5]
text = path.read_text()
text = re.sub(r"case_name='[^']*'", f"case_name='{case_name}'", text)
text = re.sub(r"toolpath_file='[^']*'", f"toolpath_file='{toolpath_file}'", text)
text = re.sub(r"(timax=)[^,]+", rf"\g<1>{timax}", text, count=1)
text = re.sub(r"(delt=)[^,]+", rf"\g<1>{delt}", text, count=1)
path.write_text(text)
PY
}

extract_tot_iter() {
  awk '/^[[:space:]]*[0-9.]+E[+-][0-9]+/ { total = $4 } END { if (total == "") exit 1; print total }' "$1"
}

run_case() {
  local case_name="$1"
  local toolpath_file="$2"
  local output_file="$RESULT_DIR/${case_name}.out"

  cp "$PROJECT_INPUT_TEMPLATE" "$PROJECT_INPUT_WORK"
  set_input "$case_name" "$toolpath_file"
  cp "$PROJECT_INPUT_WORK" "$BASE_INPUT_FILE"
  (cd "$FORTRAN_DIR" && OMP_NUM_THREADS="$OMP_THREADS" ./cluster_main > "$output_file" 2>&1)
}

summarize_case() {
  local case_name="$1"
  local output_file="$SOLVER_RESULT_DIR/${case_name}/${case_name}_output.txt"
  local total_iter
  total_iter="$(extract_tot_iter "$output_file")"
  printf '%-14s %s\n' "$case_name" "$total_iter"
}

check_pair() {
  local left_case="$1"
  local right_case="$2"
  local left_total right_total diff
  left_total="$(extract_tot_iter "$SOLVER_RESULT_DIR/${left_case}/${left_case}_output.txt")"
  right_total="$(extract_tot_iter "$SOLVER_RESULT_DIR/${right_case}/${right_case}_output.txt")"
  if (( left_total > right_total )); then
    diff=$((left_total - right_total))
  else
    diff=$((right_total - left_total))
  fi
  if (( diff > TOLERANCE )); then
    echo "FAIL: $left_case ($left_total) vs $right_case ($right_total) differ by $diff iterations" >&2
    return 1
  fi
  echo "PASS: $left_case ($left_total) vs $right_case ($right_total) differ by $diff iterations"
}

cases=(
  "reg_plus_x $SCRIPT_DIR/scan_plus_x.crs"
  "reg_minus_x $SCRIPT_DIR/scan_minus_x.crs"
  "reg_plus_y $SCRIPT_DIR/scan_plus_y.crs"
  "reg_minus_y $SCRIPT_DIR/scan_minus_y.crs"
)

for entry in "${cases[@]}"; do
  case_name="${entry%% *}"
  toolpath_file="${entry#* }"
  run_case "$case_name" "$toolpath_file"
done

printf '\nFinal total iterations\n'
for entry in "${cases[@]}"; do
  case_name="${entry%% *}"
  summarize_case "$case_name"
done

printf '\nDirectional checks\n'
check_pair reg_plus_x reg_minus_x
check_pair reg_plus_y reg_minus_y

echo "Directional symmetry regression check completed successfully."
