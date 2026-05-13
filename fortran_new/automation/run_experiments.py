#!/usr/bin/env python3
"""PHOENIX ML-ready experiment runner.

Creates isolated per-case workspaces, generates case-specific toolpaths,
runs cases in parallel batches, and collects KPI + VTK artifacts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple


DEFAULT_VTK_PATTERNS = ["*_defect.vtk", "*_maxtemp.vtk", "*_vtkmov*.vtk"]


@dataclass
class CaseRecord:
    index: int
    case_name: str
    workspace: Path
    params: Dict[str, Any]
    status: str = "pending"
    pid: str = ""
    retries_used: int = 0
    start_time_utc: str = ""
    end_time_utc: str = ""
    elapsed_s: str = ""
    exit_code: str = ""
    failure_reason: str = ""
    stdout_path: str = ""
    output_dir: str = ""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_config(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in {".yaml", ".yml"}:
        raise ValueError("YAML is not supported in this v1. Use JSON config.")
    return json.loads(text)


def resolve_path(path_str: str, repo_root: Path, fortran_dir: Path, config_dir: Path) -> Path:
    p = Path(path_str)
    if p.is_absolute():
        return p
    cands = [
        (Path.cwd() / p),
        (fortran_dir / p),
        (repo_root / p),
        (config_dir / p),
    ]
    for c in cands:
        if c.exists():
            return c.resolve()
    return (fortran_dir / p).resolve()


def sanitize(v: Any) -> str:
    txt = str(v)
    txt = txt.replace("-", "m").replace(".", "p")
    txt = re.sub(r"[^A-Za-z0-9_]+", "_", txt)
    return txt.strip("_")[:40]


def build_case_name(i: int, params: Dict[str, Any]) -> str:
    parts = [f"k_{k}_{sanitize(params[k])}" for k in sorted(params.keys())]
    return f"case_{i:04d}_" + "__".join(parts)


def cartesian_cases(sweep: Dict[str, List[Any]]) -> Iterable[Dict[str, Any]]:
    keys = list(sweep.keys())
    vals = [sweep[k] for k in keys]
    for combo in itertools.product(*vals):
        yield dict(zip(keys, combo))


def replace_namelist_value(text: str, key: str, value: Any) -> str:
    if isinstance(value, str):
        rep = f"{key}='{value}'"
    elif isinstance(value, bool):
        rep = f"{key}={1 if value else 0}"
    else:
        rep = f"{key}={value}"
    pattern = rf"(?i)\b{re.escape(key)}\s*=\s*('[^']*'|[^,\s/]+)"
    new_text, n = re.subn(pattern, rep, text)
    if n == 0:
        raise ValueError(f"Could not find key '{key}' in input template")
    return new_text


def run_cmd(cmd: List[str], cwd: Path, env: Optional[Dict[str, str]] = None) -> None:
    proc = subprocess.run(cmd, cwd=str(cwd), env=env, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")


def build_once(fortran_dir: Path) -> None:
    run_cmd(["bash", "compile.sh"], cwd=fortran_dir)


def generate_toolpath(
    fortran_dir: Path,
    case_workspace: Path,
    cfg_toolpath: Dict[str, Any],
    case_params: Dict[str, Any],
) -> Path:
    template = cfg_toolpath.get("template", {})
    sweep_keys = cfg_toolpath.get("sweep_keys", {"scan_speed": "scan_speed", "hatch_spacing": "hatch_spacing"})
    speed_key = sweep_keys.get("scan_speed", "scan_speed")
    hatch_key = sweep_keys.get("hatch_spacing", "hatch_spacing")
    if speed_key not in case_params or hatch_key not in case_params:
        raise ValueError("Sweep must include scan_speed and hatch_spacing (or mapped keys)")

    tool_dir = case_workspace / "ToolFiles"
    tool_dir.mkdir(parents=True, exist_ok=True)
    out_path = tool_dir / "case_toolpath.crs"

    cmd = [
        "python3",
        str(fortran_dir / "ToolFiles" / "toolpath_generator_rectangle.py"),
        "--output",
        str(out_path),
        "--scan_speed",
        str(case_params[speed_key]),
        "--hatch_spacing",
        str(case_params[hatch_key]),
    ]

    for key in (
        "start_x",
        "start_y",
        "start_z",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "scan_axis",
        "turnaround_time",
        "rotation_angle",
        "domain_x",
        "domain_y",
    ):
        if key in template:
            cmd.extend([f"--{key}", str(template[key])])
    if template.get("unidirectional", False):
        cmd.append("--unidirectional")
    else:
        cmd.append("--bidirectional")

    run_cmd(cmd, cwd=fortran_dir)
    return out_path


def materialize_case(
    fortran_dir: Path,
    run_id: str,
    input_template: Path,
    case_record: CaseRecord,
    toolpath_cfg: Dict[str, Any],
    static_params: Dict[str, Any],
) -> None:
    ws = case_record.workspace
    (ws / "inputfile").mkdir(parents=True, exist_ok=True)
    (ws / "result").mkdir(parents=True, exist_ok=True)

    text = input_template.read_text(encoding="utf-8")
    toolpath_path = generate_toolpath(fortran_dir, ws, toolpath_cfg, case_record.params)

    merged = dict(static_params)
    merged.update(case_record.params)
    merged["case_name"] = case_record.case_name
    merged["toolpath_file"] = str(toolpath_path)

    for key, value in merged.items():
        text = replace_namelist_value(text, key, value)

    (ws / "inputfile" / "input_param.txt").write_text(text, encoding="utf-8")

    bin_src = fortran_dir / "cluster_main"
    if not bin_src.exists():
        raise FileNotFoundError("cluster_main not found. Compile first.")
    shutil.copy2(bin_src, ws / "cluster_main")

    case_record.stdout_path = str(ws / f"{case_record.case_name}.stdout.log")
    case_record.output_dir = str(ws / "result" / case_record.case_name)


def write_manifest(path: Path, rows: List[CaseRecord]) -> None:
    fields = [
        "index",
        "case_name",
        "status",
        "pid",
        "retries_used",
        "start_time_utc",
        "end_time_utc",
        "elapsed_s",
        "exit_code",
        "failure_reason",
        "workspace",
        "stdout_path",
        "output_dir",
        "params_json",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "index": r.index,
                    "case_name": r.case_name,
                    "status": r.status,
                    "pid": r.pid,
                    "retries_used": r.retries_used,
                    "start_time_utc": r.start_time_utc,
                    "end_time_utc": r.end_time_utc,
                    "elapsed_s": r.elapsed_s,
                    "exit_code": r.exit_code,
                    "failure_reason": r.failure_reason,
                    "workspace": str(r.workspace),
                    "stdout_path": r.stdout_path,
                    "output_dir": r.output_dir,
                    "params_json": json.dumps(r.params, sort_keys=True),
                }
            )


def parse_timing_report(path: Path) -> Tuple[Optional[float], Optional[float]]:
    if not path.exists():
        return None, None
    text = path.read_text(encoding="utf-8", errors="ignore")
    cpu = re.search(r"Total CPU time:\s+([\d.]+)\s+s", text)
    wall = re.search(r"Total wall time:\s+([\d.]+)\s+s", text)
    return (float(cpu.group(1)) if cpu else None, float(wall.group(1)) if wall else None)


def parse_defect_report(path: Path) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if not path.exists():
        return None, None, None
    text = path.read_text(encoding="utf-8", errors="ignore")
    def _x(label: str) -> Optional[float]:
        m = re.search(rf"{label}:\s+([0-9.]+)\s*%", text, flags=re.IGNORECASE)
        return float(m.group(1)) if m else None
    return _x("Defect fraction"), _x("Lack-of-fusion fraction"), _x("Keyhole porosity fraction")


def parse_meltpool(path: Path) -> Dict[str, Optional[float]]:
    out = {
        "mp_len_max_m": None,
        "mp_depth_max_m": None,
        "mp_width_max_m": None,
        "mp_vol_max_m3": None,
        "tpeak_max_k": None,
    }
    if not path.exists():
        return out
    vals: List[Tuple[float, float, float, float, float]] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 6:
            continue
        try:
            vals.append((float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4]), float(parts[5])))
        except ValueError:
            continue
    if not vals:
        return out
    out["mp_len_max_m"] = max(v[0] for v in vals)
    out["mp_depth_max_m"] = max(v[1] for v in vals)
    out["mp_width_max_m"] = max(v[2] for v in vals)
    out["mp_vol_max_m3"] = max(v[3] for v in vals)
    out["tpeak_max_k"] = max(v[4] for v in vals)
    return out


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def collect_kpis_and_ml(
    cases: List[CaseRecord],
    results_dir: Path,
    vtk_patterns: List[str],
    copy_mode: str,
) -> None:
    kpi_rows: List[Dict[str, Any]] = []
    ml_rows: List[Dict[str, Any]] = []
    failures: List[Dict[str, Any]] = []

    ml_root = results_dir / "ml_data"
    ml_root.mkdir(parents=True, exist_ok=True)

    for c in cases:
        if c.status != "success":
            failures.append({"case_name": c.case_name, "status": c.status, "reason": c.failure_reason, "stdout_path": c.stdout_path})
            continue
        out_dir = Path(c.output_dir)
        prefix = out_dir / f"{c.case_name}"
        cpu_s, wall_s = parse_timing_report(prefix.with_name(prefix.name + "_timing_report.txt"))
        defect_pct, lof_pct, keyhole_pct = parse_defect_report(prefix.with_name(prefix.name + "_defect_report.txt"))
        melt = parse_meltpool(prefix.with_name(prefix.name + "_meltpool_history.txt"))
        row: Dict[str, Any] = {
            "case_name": c.case_name,
            "status": c.status,
            "cpu_time_s": cpu_s,
            "wall_time_s": wall_s,
            "defect_pct": defect_pct,
            "lof_pct": lof_pct,
            "keyhole_pct": keyhole_pct,
            "output_dir": str(out_dir),
        }
        row.update(melt)
        for k, v in sorted(c.params.items()):
            row[f"param_{k}"] = v
        kpi_rows.append(row)

        case_ml = ml_root / c.case_name
        case_ml.mkdir(parents=True, exist_ok=True)
        for pat in vtk_patterns:
            for src in out_dir.glob(pat):
                dst = case_ml / src.name
                if dst.exists():
                    dst.unlink()
                if copy_mode == "symlink":
                    os.symlink(src.resolve(), dst)
                else:
                    shutil.copy2(src, dst)
                ml_rows.append(
                    {
                        "case_name": c.case_name,
                        "source_path": str(src),
                        "stored_path": str(dst),
                        "bytes": dst.stat().st_size,
                        "sha256": sha256_file(dst),
                    }
                )

    if kpi_rows:
        with (results_dir / "kpi_summary.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=sorted(kpi_rows[0].keys()))
            w.writeheader()
            w.writerows(kpi_rows)
    else:
        (results_dir / "kpi_summary.csv").write_text("", encoding="utf-8")

    with (results_dir / "failures.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["case_name", "status", "reason", "stdout_path"])
        w.writeheader()
        w.writerows(failures)

    with (results_dir / "ml_index.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["case_name", "source_path", "stored_path", "bytes", "sha256"])
        w.writeheader()
        w.writerows(ml_rows)


def launch_case(case: CaseRecord, execution: Dict[str, Any]) -> subprocess.Popen:
    thermal_threads = int(execution.get("thermal_threads", 4))
    mech_threads = int(execution.get("mech_threads", 0))
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(thermal_threads)
    env["PHOENIX_THERMAL_THREADS"] = str(thermal_threads)
    env["PHOENIX_MECH_THREADS"] = str(mech_threads)
    env["PHOENIX_RUN_MODE"] = "single"
    out = open(case.stdout_path, "w", encoding="utf-8")
    cmd = ["./cluster_main"]
    return subprocess.Popen(cmd, cwd=str(case.workspace), env=env, stdout=out, stderr=subprocess.STDOUT, start_new_session=True)


def execute_cases(cases: List[CaseRecord], manifest_path: Path, execution: Dict[str, Any]) -> None:
    max_concurrent = int(execution["max_concurrent_runs"])
    retries = int(execution.get("retries", 0))
    timeout_s = int(execution.get("timeout_s", 0))
    running: Dict[int, Tuple[subprocess.Popen, float]] = {}
    pending = [c for c in cases if c.status in {"pending", "failed_retry"}]

    try:
        while pending or running:
            while pending and len(running) < max_concurrent:
                c = pending.pop(0)
                c.status = "running"
                c.start_time_utc = now_utc()
                p = launch_case(c, execution)
                c.pid = str(p.pid)
                running[c.index] = (p, time.time())
                write_manifest(manifest_path, cases)

            finished: List[int] = []
            for idx, (proc, t0) in running.items():
                rc = proc.poll()
                elapsed = time.time() - t0
                c = cases[idx]
                if timeout_s > 0 and rc is None and elapsed > timeout_s:
                    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                    c.status = "failed"
                    c.failure_reason = f"timeout>{timeout_s}s"
                    c.exit_code = "124"
                    c.end_time_utc = now_utc()
                    c.elapsed_s = f"{elapsed:.2f}"
                    finished.append(idx)
                    continue
                if rc is None:
                    continue
                c.end_time_utc = now_utc()
                c.elapsed_s = f"{elapsed:.2f}"
                c.exit_code = str(rc)
                if rc == 0:
                    c.status = "success"
                else:
                    if c.retries_used < retries:
                        c.retries_used += 1
                        c.status = "failed_retry"
                        pending.append(c)
                    else:
                        c.status = "failed"
                        c.failure_reason = f"exit_code={rc}"
                finished.append(idx)

            for idx in finished:
                running.pop(idx, None)
            if finished:
                write_manifest(manifest_path, cases)
            time.sleep(1.0)
    except KeyboardInterrupt:
        for proc, _ in running.values():
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            except Exception:
                pass
        raise


def git_hash(repo_root: Path) -> str:
    try:
        out = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(repo_root), text=True)
        return out.strip()
    except Exception:
        return ""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="Path to JSON config")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--max-cases", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-compile", action="store_true")
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    fortran_dir = script_dir.parent
    repo_root = fortran_dir.parent
    config_path = Path(args.config).resolve()
    cfg = load_config(config_path)

    sweep = cfg["sweep"]
    execution = cfg["execution"]
    toolpath_cfg = cfg.get("toolpath", {})
    static_params = cfg.get("static_params", {})
    ml_cfg = cfg.get("ml_outputs", {})
    vtk_patterns = ml_cfg.get("vtk_patterns", DEFAULT_VTK_PATTERNS)
    copy_mode = ml_cfg.get("copy_mode", "copy")
    if copy_mode not in {"copy", "symlink"}:
        raise ValueError("ml_outputs.copy_mode must be 'copy' or 'symlink'")

    run_root = fortran_dir / "automation" / "runs" / args.run_id
    result_root = fortran_dir / "automation" / "results" / args.run_id
    run_root.mkdir(parents=True, exist_ok=True)
    result_root.mkdir(parents=True, exist_ok=True)
    manifest_path = result_root / "cases_manifest.csv"

    input_template = resolve_path(
        cfg.get("input_template", "inputfile/input_param.txt"),
        repo_root=repo_root,
        fortran_dir=fortran_dir,
        config_dir=config_path.parent,
    )

    combos = list(cartesian_cases(sweep))
    if args.max_cases > 0:
        combos = combos[: args.max_cases]

    cases: List[CaseRecord] = []
    for i, params in enumerate(combos):
        case_name = build_case_name(i, params)
        ws = run_root / case_name
        cases.append(CaseRecord(index=i, case_name=case_name, workspace=ws, params=params))

    if args.resume and manifest_path.exists():
        old = {}
        with manifest_path.open("r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                old[row["case_name"]] = row
        for c in cases:
            if c.case_name in old:
                o = old[c.case_name]
                c.status = o.get("status", c.status)
                c.pid = o.get("pid", "")
                c.retries_used = int(o.get("retries_used", "0") or 0)
                c.start_time_utc = o.get("start_time_utc", "")
                c.end_time_utc = o.get("end_time_utc", "")
                c.elapsed_s = o.get("elapsed_s", "")
                c.exit_code = o.get("exit_code", "")
                c.failure_reason = o.get("failure_reason", "")
                c.stdout_path = o.get("stdout_path", "")
                c.output_dir = o.get("output_dir", "")
                if c.status == "running":
                    c.status = "failed_retry"
                    c.failure_reason = "resume_interrupted_run"

    if not args.skip_compile and not args.dry_run:
        build_once(fortran_dir)

    for c in cases:
        if args.resume and c.status in {"success"}:
            continue
        materialize_case(fortran_dir, args.run_id, input_template, c, toolpath_cfg, static_params)
        if c.status not in {"success"}:
            c.status = "pending"

    write_manifest(manifest_path, cases)

    metadata = {
        "run_id": args.run_id,
        "created_utc": now_utc(),
        "config_path": str(config_path),
        "git_commit": git_hash(repo_root),
        "execution": execution,
        "ml_outputs": {"vtk_patterns": vtk_patterns, "copy_mode": copy_mode},
        "num_cases": len(cases),
    }
    (result_root / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    if args.dry_run:
        print(f"Dry run complete. Materialized {len(cases)} cases at {run_root}")
        return

    execute_cases(cases, manifest_path, execution)
    collect_kpis_and_ml(cases, result_root, vtk_patterns, copy_mode)
    write_manifest(manifest_path, cases)
    print(f"Completed run_id={args.run_id}. Results: {result_root}")


if __name__ == "__main__":
    main()
