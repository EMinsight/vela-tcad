#!/usr/bin/env python3
"""Freeze the five non-transient Sentaurus BVmethods branches and Vela ledger."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import analyze_sentaurus_bvmethods as analyzer  # noqa: E402


METHODS = {
    "ABA_poisson": 3,
    "ABA_coupled": 4,
    "resistor": 5,
    "voltage2current": 6,
    "continuation": 7,
}
EXPECTED_BV_V = {
    "ABA_poisson": 5.305525632989282,
    "ABA_coupled": 6.377494277837012,
    "resistor": 6.379791636301563,
    "voltage2current": 6.38318420057198,
    "continuation": 6.383727168968036,
}


def write_text_lf(path: Path, text: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(text)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_continuation_bv(path: Path, target: float = 1.0e-4) -> dict[str, Any]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = [row for row in csv.DictReader(stream)
                if row.get("converged", "1") == "1"]
    bias_name = "inner_voltage_V" if rows and "inner_voltage_V" in rows[0] else "bias_V"
    current_name = (
        "current_total_A_per_um" if rows and "current_total_A_per_um" in rows[0]
        else "current_total")
    biases = [float(row[bias_name]) for row in rows]
    currents = [abs(float(row[current_name])) for row in rows]
    crossing = analyzer.first_upward_crossing(biases, currents, target)
    reference = EXPECTED_BV_V["continuation"]
    relative = None if crossing is None else abs(crossing - reference) / reference
    return {
        "curve": path.name,
        "points": len(rows),
        "target_current_A_per_um": target,
        "vela_bv_V": crossing,
        "sentaurus_bv_V": reference,
        "relative_error": relative,
        "relative_error_limit": 0.02,
        "status": "pass" if relative is not None and relative <= 0.02 else "fail",
    }


def freeze(
    sentaurus_dir: Path,
    fixture: Path,
    iic_summary: Path,
    full_physics_summary: Path,
    continuation_csv: Path | None = None,
) -> dict[str, Any]:
    curves = fixture / "reference_curves"
    curves.mkdir(parents=True, exist_ok=True)
    methods = []
    hashes: dict[str, str] = {}
    for method, node in METHODS.items():
        source = sentaurus_dir / f"n{node}_des.plt"
        datasets, rows = analyzer.read_plt(source)
        curve = curves / f"{method}.csv"
        analyzer.write_curve(curve, datasets, rows)
        result = analyzer.analyze_method(
            method, datasets, rows,
            analyzer.DEFAULT_CURRENT_THRESHOLD_A_PER_UM,
            analyzer.DEFAULT_ION_INTEGRAL_THRESHOLD)
        expected = EXPECTED_BV_V[method]
        extracted = result.get("bv_V")
        result.update({
            "node": node,
            "source_artifact": source.name,
            "curve": f"reference_curves/{curve.name}",
            "expected_bv_V": expected,
            "reference_matches": extracted is not None and abs(extracted - expected) <= 1.0e-9,
        })
        methods.append(result)
        hashes[f"reference_curves/{curve.name}"] = sha256(curve)
        hashes[f"runtime/{source.name}"] = sha256(source)

    for path in sorted((fixture / "source").iterdir()):
        if path.is_file():
            hashes[f"source/{path.name}"] = sha256(path)

    iic = json.loads(iic_summary.read_text(encoding="utf-8"))
    full = json.loads(full_physics_summary.read_text(encoding="utf-8"))
    continuation = (
        extract_continuation_bv(continuation_csv)
        if continuation_csv is not None and continuation_csv.is_file()
        else {
            "status": "pending",
            "reason": "no converged Vela BVmethods NMOS arclength curve supplied",
            "sentaurus_bv_V": EXPECTED_BV_V["continuation"],
        }
    )
    continuation_diagnostic = fixture / "continuation_diagnostic_20260818.json"
    if continuation["status"] != "pass" and continuation_diagnostic.is_file():
        diagnostic = json.loads(
            continuation_diagnostic.read_text(encoding="utf-8"))
        continuation["diagnostic"] = continuation_diagnostic.name
        continuation["classification"] = diagnostic["conclusion"]
        hashes[continuation_diagnostic.name] = sha256(continuation_diagnostic)
    result = {
        "schema": "vela.bvmethods.nontransient_freeze.v1",
        "scope": list(METHODS),
        "excluded": ["transient"],
        "sentaurus_version": "O-2018.06-SP2",
        "sentaurus_reference": {
            "status": "pass" if all(item["reference_matches"] for item in methods) else "fail",
            "methods": methods,
        },
        "vela_acceptance": {
            "aba_poisson": {
                "status": "reference_operator_only",
                "reason": "Poisson ABA is retained as an ionization-path operator reference",
            },
            "aba_coupled_iic": iic,
            "external_resistor": full["external_resistor_cross_check"],
            "voltage_to_current": full["sentaurus_full_model_acceptance"],
            "continuation": continuation,
        },
        "status": (
            "pass" if all(item["reference_matches"] for item in methods)
            and iic["path_iic"]["status"] == "pass"
            and full["status"] == "PASS"
            and continuation["status"] == "pass"
            else "pass_with_continuation_pending"
            if all(item["reference_matches"] for item in methods)
            and iic["path_iic"]["status"] == "pass"
            and full["status"] == "PASS"
            and continuation["status"] == "pending"
            else "fail"
        ),
        "sha256": hashes,
    }
    output = fixture / "bvmethods_nontransient_validation_20260817.json"
    write_text_lf(output, json.dumps(result, indent=2) + "\n")
    write_markdown(
        fixture / "bvmethods_nontransient_validation_20260817.md", result)
    return result


def write_markdown(path: Path, result: dict[str, Any]) -> None:
    lines = [
        "# BVmethods non-transient freeze", "",
        f"Status: **{result['status'].upper()}**", "",
        "Transient node 8 is explicitly outside this scope.", "",
        "| Method | Node | Rows | Sentaurus BV | Reference frozen | Vela mapping |",
        "| --- | ---: | ---: | ---: | --- | --- |",
    ]
    mapping = {
        "ABA_poisson": "operator reference",
        "ABA_coupled": "path/current IIC",
        "resistor": "external resistor",
        "voltage2current": "voltage to current",
        "continuation": "pseudo-arclength",
    }
    for item in result["sentaurus_reference"]["methods"]:
        lines.append(
            f"| {item['method']} | {item['node']} | {item['rows']} | "
            f"{item['bv_V']:.9g} V | {'yes' if item['reference_matches'] else 'no'} | "
            f"{mapping[item['method']]} |")
    continuation = result["vela_acceptance"]["continuation"]
    lines.extend([
        "", "## Vela acceptance", "",
        f"- Path IIC: `{result['vela_acceptance']['aba_coupled_iic']['path_iic']['status']}`.",
        f"- Current IIC: `{result['vela_acceptance']['aba_coupled_iic']['current_iic']['status']}`.",
        f"- External resistor: `{result['vela_acceptance']['external_resistor']['status']}`.",
        f"- Voltage to current: `{result['vela_acceptance']['voltage_to_current']['status']}`.",
        f"- Continuation: `{continuation['status']}`.", "",
    ])
    if continuation["status"] != "pass":
        lines.append(
            "The Sentaurus continuation reference is frozen, but Vela continuation "
            "is not claimed complete until a converged NMOS arclength curve is supplied.")
        lines.append("")
        if continuation.get("diagnostic"):
            lines.append(
                f"The bounded numerical trials are recorded in "
                f"`{continuation['diagnostic']}`; no physical parameter was changed.")
            lines.append("")
    write_text_lf(path, "\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument(
        "--fixture", type=Path,
        default=ROOT / "reference_tcad" / "bvmethods_sentaurus2018")
    parser.add_argument("--iic-summary", type=Path, required=True)
    parser.add_argument("--full-physics-summary", type=Path, required=True)
    parser.add_argument("--continuation-csv", type=Path)
    args = parser.parse_args()
    result = freeze(
        args.sentaurus_dir.resolve(), args.fixture.resolve(),
        args.iic_summary.resolve(), args.full_physics_summary.resolve(),
        args.continuation_csv.resolve() if args.continuation_csv else None)
    print(json.dumps({"status": result["status"], "scope": result["scope"]}))
    return 0 if result["status"] != "fail" else 1


if __name__ == "__main__":
    raise SystemExit(main())
