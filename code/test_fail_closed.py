#!/usr/bin/env python3
"""Focused negative regressions for previously demonstrated false PASS paths."""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from mpmath import mp

ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable


def run(*arguments):
    return subprocess.run(
        [PYTHON, *map(str, arguments)], cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )


def require(condition: bool, message: str, result=None):
    if not condition:
        detail = "" if result is None else "\n" + result.stdout
        raise SystemExit(f"FAIL-CLOSED REGRESSION FAILED: {message}{detail}")


def det_at_origin(data: dict):
    mp.dps = 80
    A = [[mp.mpf(str(value)) for value in row] for row in data["A"]]
    B = [[mp.mpf(str(value)) for value in row] for row in data["B"]]
    c = [mp.mpf(str(value)) for value in data["c"]]
    weights = [1 / (2 + mp.exp(value) + mp.exp(-value)) for value in c]
    matrix = [[sum(A[row][i] * weights[i] * B[i][column] for i in range(4)) for column in range(2)] for row in range(2)]
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def main() -> int:
    scratch = ROOT / '_local'
    scratch.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="njc-fail-closed-", dir=scratch) as temporary:
        temp = Path(temporary)

        witness = json.loads((ROOT / "data" / "witness.json").read_text(encoding="utf-8"))
        for row in range(2):
            witness["A"][row][1] = -witness["A"][row][1]
        bad_witness = temp / "bad-witness.json"
        bad_witness.write_text(json.dumps(witness, indent=2) + "\n", encoding="utf-8", newline="\n")
        analytic_output = temp / "analytic.json"
        result = run(ROOT / "code" / "certify_analytic_determinant.py", "--witness", bad_witness, "--json-out", analytic_output)
        require(det_at_origin(witness) < 0, "analytic plant no longer has a negative determinant")
        require(result.returncode != 0, "analytic checker accepted omitted nonnegative-minor premises", result)
        require(json.loads(analytic_output.read_text())["verdict"] == "FAIL", "analytic failure report is not FAIL")
        print("[PASS] omitted analytic premises fail closed")

        empty_state = temp / "empty-state.json"
        empty_state.write_text(json.dumps({"processed": 0, "certified": 0, "paths": [], "worst": "+inf", "stack": []}), encoding="utf-8", newline="\n")
        result = run(
            ROOT / "code" / "certify_parameter_radius.py", "--delta", "0", "--skip-tail",
            "--state-file", empty_state, "--json-out", temp / "empty-state-report.json",
        )
        require(result.returncode != 0, "empty compact resume state manufactured PASS", result)
        print("[PASS] forged compact state fails closed")

        tail_state = temp / "tail-state.json"
        tail_state.write_text(json.dumps({"next_k": 200000, "worst": "1", "worstk": 0}), encoding="utf-8", newline="\n")
        result = run(
            ROOT / "code" / "certify_parameter_radius.py", "--delta", "0", "--skip-interior",
            "--tail-state-file", tail_state, "--json-out", temp / "tail-state-report.json",
        )
        require(result.returncode != 0, "forged tail resume state manufactured PASS", result)
        print("[PASS] forged tail state fails closed")

        pause_state = temp / "pause-state.json"
        pause_report = temp / "pause-report.json"
        result = run(
            ROOT / "code" / "certify_parameter_radius.py", "--delta", "0", "--skip-tail",
            "--state-file", pause_state, "--chunk-seconds", "0.000001", "--json-out", pause_report,
        )
        require(result.returncode == 2, "PAUSED did not return process status 2", result)
        require(json.loads(pause_report.read_text())["verdict"] == "PAUSED", "pause report is not PAUSED")
        print("[PASS] PAUSED is not process success")

        saved_state = temp / "saved-valid-state.json"
        shutil.copy2(pause_state, saved_state)
        result = run(
            ROOT / "code" / "certify_parameter_radius.py", "--delta", "0", "--skip-tail",
            "--state-file", pause_state, "--chunk-seconds", "0.000001", "--json-out", pause_report,
        )
        require(result.returncode in (0, 2), "valid compact checkpoint did not resume", result)
        require(json.loads(pause_report.read_text())["verdict"] in ("PASS", "PAUSED"), "valid resume returned failure")
        print("[PASS] valid context-bound compact state resumes")

        result = run(
            ROOT / "code" / "certify_parameter_radius.py", "--delta", "0.000001", "--skip-tail",
            "--state-file", saved_state, "--json-out", temp / "cross-context-report.json",
        )
        require(result.returncode != 0, "checkpoint copied across parameter contexts was accepted", result)
        print("[PASS] cross-context compact state fails closed")

        report_copy = temp / "reports"
        shutil.copytree(ROOT / "results" / "witness_common_5e-6", report_copy)
        compact_path = report_copy / "compact.json"
        compact = json.loads(compact_path.read_text())
        compact["interior"]["passed"] = False
        compact["interior"]["unresolved"] = 1
        compact_path.write_text(json.dumps(compact, indent=2) + "\n", encoding="utf-8", newline="\n")
        arb_path = report_copy / "compact_arb.json"
        arb_report = json.loads(arb_path.read_text())
        arb_report["compact_report_sha256"] = hashlib.sha256(compact_path.read_bytes()).hexdigest()
        arb_path.write_text(json.dumps(arb_report, indent=2) + "\n", encoding="utf-8", newline="\n")
        result = run(
            ROOT / "code" / "rebuild_open_family_aggregate.py", report_copy,
            "--config", ROOT / "config" / "witness_common.json",
        )
        require(result.returncode != 0, "contradictory compact report passed aggregate rebuild", result)
        print("[PASS] contradictory proof report fails semantic validation")

        omitted_copy = temp / "reports-omitted"
        shutil.copytree(ROOT / "results" / "witness_common_5e-6", omitted_copy)
        (omitted_copy / "tail_independent.json").unlink()
        result = run(
            ROOT / "code" / "rebuild_open_family_aggregate.py", omitted_copy,
            "--config", ROOT / "config" / "witness_common.json",
        )
        require(result.returncode != 0, "omitted expected component passed aggregate rebuild", result)
        print("[PASS] omitted aggregate component fails closed")

        lines = (ROOT / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines()
        incomplete_manifest = temp / "SHA256SUMS-incomplete.txt"
        incomplete_manifest.write_text("\n".join(lines[1:]) + "\n", encoding="utf-8", newline="\n")
        result = run(ROOT / "code" / "verify_sha256_manifest.py", "--manifest", incomplete_manifest)
        require(result.returncode != 0, "incomplete manifest passed verification", result)
        print("[PASS] incomplete manifest fails closed")

    print("VERDICT: FOCUSED FAIL-CLOSED REGRESSIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
