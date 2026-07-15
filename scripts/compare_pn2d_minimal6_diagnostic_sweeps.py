#!/usr/bin/env python3
"""Compare exact PN2D Minimal6 diagnostic-sweep checkpoints without interpolation."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts.pn2d_minimal6_diagnostics.schemas import DISCLAIMER, validate_bv_comparison_v1

SCHEMA = "vela.pn2d_minimal6_bv_comparison.v1"
EPSILON = 1.0e-12
OBSERVABLES = (
    "anode_current_A_per_um", "cathode_current_A_per_um", "max_field_V_per_m",
    "native_source_integral_s_inv_per_cm", "reconstructed_source_integral_s_inv_per_cm",
)


def _sha_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical_sha(payload: Any) -> str:
    return _sha_bytes(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8"))


def _finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def ratio_record(numerator: float | None, denominator: float | None) -> dict[str, Any]:
    """Return a typed ratio; zero and unavailable values are never coerced."""
    if numerator is None or denominator is None or not _finite(numerator) or not _finite(denominator):
        return {"classification": "unavailable", "value": None}
    if denominator == 0.0:
        return {"classification": "zero_denominator", "value": None}
    if numerator == 0.0:
        return {"classification": "zero_numerator", "value": 0.0}
    return {"classification": "available", "value": float(numerator / denominator)}


def _accepted(manifest: dict[str, Any], solver: str) -> list[dict[str, Any]]:
    rows = manifest.get("accepted_checkpoints", [])
    if not isinstance(rows, list):
        raise ValueError(f"{solver} manifest has invalid accepted_checkpoints")
    checked: list[dict[str, Any]] = []
    seen: set[tuple[str, float]] = set()
    for row in rows:
        if row.get("solver") != solver or row.get("status") != "accepted":
            raise ValueError(f"{solver} manifest has invalid accepted checkpoint")
        actual, target = row.get("actual_bias_V"), row.get("target_bias_V")
        if not _finite(actual) or not _finite(target) or abs(float(actual) - float(target)) > EPSILON:
            raise ValueError(f"{solver} checkpoint is not an exact target bias")
        topology = row.get("topology")
        if not isinstance(topology, str) or not topology:
            raise ValueError(f"{solver} checkpoint lacks topology")
        key = (topology, float(target))
        if key in seen:
            raise ValueError(f"{solver} has duplicate exact checkpoint {key}")
        seen.add(key)
        observables = row.get("observables")
        if not isinstance(observables, dict) or any(not _finite(observables.get(name)) for name in OBSERVABLES):
            raise ValueError(f"{solver} checkpoint lacks finite observables")
        checked.append(row)
    return checked


def _index(rows: list[dict[str, Any]]) -> dict[tuple[str, float], dict[str, Any]]:
    return {(str(row["topology"]), float(row["target_bias_V"])): row for row in rows}


def _current(row: dict[str, Any]) -> float:
    return float(row["observables"]["anode_current_A_per_um"])


def _ratio_rows(numerator: dict[str, Any], denominator: dict[str, Any], observable: str, *, absolute: bool = False) -> dict[str, Any]:
    left = float(numerator["observables"][observable])
    right = float(denominator["observables"][observable])
    return ratio_record(abs(left) if absolute else left, abs(right) if absolute else right)


def _one_volt_growth(index: dict[tuple[str, float], dict[str, Any]], topology: str, bias: float) -> dict[str, Any]:
    current, next_row = index.get((topology, bias)), index.get((topology, bias - 1.0))
    if current is None or next_row is None:
        return {"classification": "unavailable", "value": None, "reason": "next exact one-volt checkpoint unavailable"}
    return ratio_record(abs(_current(next_row)), abs(_current(current)))


def _failures(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    rows = manifest.get("failed_transitions", [])
    if not isinstance(rows, list):
        raise ValueError("failed_transitions must be a list")
    for row in rows:
        if row.get("status") != "rejected" or row.get("observables") is not None:
            raise ValueError("failure transition has fabricated observables")
    return rows


def _deck_hashes(manifest: dict[str, Any], key: str) -> list[str]:
    rows = manifest.get(key, [])
    if not isinstance(rows, list):
        raise ValueError(f"{key} must be a list")
    hashes = {str(row["deck_sha256"]) for row in rows if isinstance(row, dict) and isinstance(row.get("deck_sha256"), str)}
    return sorted(hashes)

def _fixed_state_recheck(common: dict[tuple[str, float], tuple[dict[str, Any], dict[str, Any]]], fixed: dict[str, Any]) -> list[dict[str, Any]]:
    original = fixed.get("root_cause_status", "unavailable")
    reason = fixed.get("root_cause_reason", "fixed-state report was not supplied")
    rows: list[dict[str, Any]] = []
    for bias in (0.0, -12.0, -19.0):
        topologies = [topology for topology in ("sketch", "mirror") if (topology, bias) in common]
        if not topologies:
            detail = "no self-consistent exact common checkpoint; quantity ledger was not rerun"
        else:
            detail = f"quantity-ledger raw state export is required to re-evaluate: {reason}"
        rows.append({"bias_V": bias, "status": "unidentifiable", "fixed_state_status": original,
                     "reason": detail, "topologies": topologies})
    return rows


def compare_sweeps(vela_manifest: dict[str, Any], sentaurus_manifest: dict[str, Any], *, fixed_state_report: dict[str, Any]) -> dict[str, Any]:
    """Produce a comparison object from accepted, exact-bias rows only."""
    vela_rows, sentaurus_rows = _accepted(vela_manifest, "vela"), _accepted(sentaurus_manifest, "sentaurus")
    vela_index, sentaurus_index = _index(vela_rows), _index(sentaurus_rows)
    keys = sorted(set(vela_index) & set(sentaurus_index), key=lambda key: (key[0], -key[1]))
    common = {key: (vela_index[key], sentaurus_index[key]) for key in keys}
    checkpoints: list[dict[str, Any]] = []
    for topology, bias in keys:
        vela, sentaurus = common[(topology, bias)]
        checkpoints.append({
            "topology": topology, "bias_V": bias, "classification": "common_exact", "vela": vela, "sentaurus": sentaurus,
            "terminal_current_sign_alignment": "aligned" if math.copysign(1.0, _current(vela)) == math.copysign(1.0, _current(sentaurus)) else "opposed",
            "terminal_current_ratio": _ratio_rows(vela, sentaurus, "anode_current_A_per_um", absolute=True),
            "maximum_field_ratio": _ratio_rows(vela, sentaurus, "max_field_V_per_m"),
            "native_source_ratio": _ratio_rows(vela, sentaurus, "native_source_integral_s_inv_per_cm"),
            "reconstructed_source_ratio": _ratio_rows(vela, sentaurus, "reconstructed_source_integral_s_inv_per_cm"),
            "vela_one_volt_current_growth": _one_volt_growth(vela_index, topology, bias),
            "sentaurus_one_volt_current_growth": _one_volt_growth(sentaurus_index, topology, bias),
            "gap_closure": {"named_contributions": ["terminal_current", "maximum_field", "native_source", "reconstructed_source"],
                            "residual": {"classification": "unidentifiable", "value": None, "reason": "no cross-solver counterfactual substitution is available"}},
        })
    side_only: list[dict[str, Any]] = []
    missing_tails: list[dict[str, Any]] = []
    for solver, index, other in (("vela", vela_index, sentaurus_index), ("sentaurus", sentaurus_index, vela_index)):
        by_topology: dict[str, list[float]] = {}
        for (topology, bias), row in index.items():
            if (topology, bias) not in other:
                side_only.append({"solver": solver, "topology": topology, "bias_V": bias, "classification": "side_only", "checkpoint": row})
                by_topology.setdefault(topology, []).append(bias)
        for topology, biases in by_topology.items():
            missing_tails.append({"solver": solver, "topology": topology, "biases_V": sorted(biases, reverse=True), "reason": "no accepted exact checkpoint from the other solver"})
    topology_sensitivity: list[dict[str, Any]] = []
    for solver, index in (("vela", vela_index), ("sentaurus", sentaurus_index)):
        biases = sorted({bias for topology, bias in index if topology == "sketch" and ("mirror", bias) in index}, reverse=True)
        for bias in biases:
            sketch, mirror = index[("sketch", bias)], index[("mirror", bias)]
            topology_sensitivity.append({"solver": solver, "bias_V": bias,
                "terminal_current_sketch_over_mirror": ratio_record(abs(_current(sketch)), abs(_current(mirror))),
                "maximum_field_sketch_over_mirror": _ratio_rows(sketch, mirror, "max_field_V_per_m"),
                "native_source_sketch_over_mirror": _ratio_rows(sketch, mirror, "native_source_integral_s_inv_per_cm")})
    failures = _failures(vela_manifest) + _failures(sentaurus_manifest)
    deepest = min((bias for _, bias in keys), default=None)
    report = {
        "schema": SCHEMA, "diagnostic_disclaimer": DISCLAIMER, "interpolation": "forbidden",
        "solver_configurations": {
            "vela": {"template": vela_manifest.get("template"), "topology_input_sha256": vela_manifest.get("topology_input_sha256", {}), "deck_sha256": _deck_hashes(vela_manifest, "segments")},
            "sentaurus": {"template": sentaurus_manifest.get("template"), "topology_input_sha256": sentaurus_manifest.get("topology_input_sha256", {}), "deck_sha256": _deck_hashes(sentaurus_manifest, "sentaurus_segments")}},
        "accepted_transitions": {"vela": vela_rows, "sentaurus": sentaurus_rows},
        "failed_transitions": failures, "failure_transitions": failures,
        "checkpoints": checkpoints, "records": checkpoints, "terminal_currents": checkpoints, "maximum_fields": checkpoints, "source_integrals": checkpoints,
        "convergence_metadata": {"vela_accepted": len(vela_rows), "sentaurus_accepted": len(sentaurus_rows), "common_exact": len(checkpoints)},
        "curve_artifact_hashes": {"vela_manifest": _canonical_sha(vela_manifest), "sentaurus_manifest": _canonical_sha(sentaurus_manifest)},
        "deepest_common_bias_V": {"classification": "available", "value": deepest} if deepest is not None else {"classification": "unavailable", "value": None, "reason": "no accepted exact common checkpoint"},
        "missing_tails": sorted(missing_tails, key=lambda row: (row["solver"], row["topology"])),
        "side_only_checkpoints": sorted(side_only, key=lambda row: (row["solver"], row["topology"], -row["bias_V"])),
        "topology_sensitivity": topology_sensitivity,
        "fixed_state_recheck": _fixed_state_recheck(common, fixed_state_report),
        "artifact_hashes": {},
        "input_artifacts": {},
        "closure": {"status": "closed", "eligible_gaps": len(checkpoints), "rule": "each eligible gap records named contributions and a typed residual; non-common points are side-only"},
    }
    return report


def _finish(fig: plt.Figure, ax: plt.Axes, title: str, ylabel: str) -> None:
    ax.set_title(title); ax.set_ylabel(ylabel); ax.set_xlabel("exact applied bias (V)"); ax.grid(True, alpha=0.25)
    fig.text(0.01, 0.01, DISCLAIMER + "; solver termination is marked; no BV extrapolation", fontsize=7)
    fig.tight_layout(rect=(0, 0.04, 1, 1))


def _side_plot(ax: plt.Axes, report: dict[str, Any], observable: str, title: str, ylabel: str, absolute: bool = False) -> None:
    any_data = False
    for solver, rows in report["accepted_transitions"].items():
        for topology in ("sketch", "mirror"):
            selected = sorted((row for row in rows if row["topology"] == topology), key=lambda row: row["target_bias_V"], reverse=True)
            if selected:
                x = [row["target_bias_V"] for row in selected]; y = [float(row["observables"][observable]) for row in selected]
                ax.plot(x, [abs(value) for value in y] if absolute else y, marker="o", label=f"{solver} {topology}"); any_data = True
    for failure in report["failure_transitions"]:
        ax.axvline(float(failure["target_bias_V"]), color="tab:red", linestyle="--", alpha=0.55)
    if any_data: ax.legend(loc="best")
    else: ax.text(0.5, 0.5, "No accepted checkpoint", ha="center", va="center", transform=ax.transAxes)
    _finish(ax.figure, ax, title, ylabel)


def _render_figures(out_dir: Path, report: dict[str, Any]) -> list[Path]:
    figures: list[Path] = []
    for stem, observable, title, unit, absolute in (
        ("terminal_current", "anode_current_A_per_um", "Terminal current at accepted exact checkpoints", "A/um", False),
        ("maximum_field", "max_field_V_per_m", "Maximum electric field at accepted exact checkpoints", "V/m", False),
        ("source_integrals", "native_source_integral_s_inv_per_cm", "Native avalanche source at accepted exact checkpoints", "s^-1 per 1 cm depth", True)):
        fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=120); _side_plot(ax, report, observable, title, unit, absolute)
        path = out_dir / f"{stem}.png"; fig.savefig(path, dpi=120); plt.close(fig); figures.append(path)
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=120)
    growth = [row for row in report["checkpoints"] if row["vela_one_volt_current_growth"]["classification"] == "available"]
    if growth:
        ax.plot([row["bias_V"] for row in growth], [row["vela_one_volt_current_growth"]["value"] for row in growth], marker="o", label="Vela")
        ax.plot([row["bias_V"] for row in growth], [row["sentaurus_one_volt_current_growth"]["value"] for row in growth], marker="s", label="Sentaurus"); ax.legend(loc="best")
    else: ax.text(0.5, 0.5, "No exact common one-volt pair", ha="center", va="center", transform=ax.transAxes)
    _finish(fig, ax, "One-volt terminal-current growth", "growth ratio")
    path = out_dir / "one_volt_growth.png"; fig.savefig(path, dpi=120); plt.close(fig); figures.append(path)
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=120)
    topology = [row for row in report["topology_sensitivity"] if row["terminal_current_sketch_over_mirror"]["classification"] == "available"]
    if topology:
        for solver in ("vela", "sentaurus"):
            selected = [row for row in topology if row["solver"] == solver]
            ax.plot([row["bias_V"] for row in selected], [row["terminal_current_sketch_over_mirror"]["value"] for row in selected], marker="o", label=solver)
        ax.legend(loc="best")
    else: ax.text(0.5, 0.5, "No exact sketch/mirror pair", ha="center", va="center", transform=ax.transAxes)
    _finish(fig, ax, "Sketch/mirror terminal-current sensitivity", "sketch / mirror ratio")
    path = out_dir / "topology.png"; fig.savefig(path, dpi=120); plt.close(fig); figures.append(path)
    return figures


def _write_csv(path: Path, report: dict[str, Any]) -> None:
    rows: list[dict[str, Any]] = []
    for row in report["checkpoints"]:
        rows.append({"classification": row["classification"], "solver": "both", "topology": row["topology"], "bias_V": row["bias_V"], "terminal_current_ratio": row["terminal_current_ratio"]["value"], "maximum_field_ratio": row["maximum_field_ratio"]["value"], "native_source_ratio": row["native_source_ratio"]["value"], "reason": "exact common checkpoint"})
    for row in report["side_only_checkpoints"]:
        rows.append({"classification": "side_only", "solver": row["solver"], "topology": row["topology"], "bias_V": row["bias_V"], "terminal_current_ratio": None, "maximum_field_ratio": None, "native_source_ratio": None, "reason": "no accepted exact checkpoint from the other solver"})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["classification", "solver", "topology", "bias_V", "terminal_current_ratio", "maximum_field_ratio", "native_source_ratio", "reason"])
        writer.writeheader(); writer.writerows(rows)


def _markdown(report: dict[str, Any]) -> str:
    deepest = report["deepest_common_bias_V"]
    deepest_text = str(deepest["value"]) if deepest["classification"] == "available" else f"unavailable ({deepest['reason']})"
    lines = ["# PN2D minimal6 diagnostic sweep comparison", "", DISCLAIMER, "", "## Exact-checkpoint result", "", f"- Deepest common accepted bias: {deepest_text}.", f"- Common exact checkpoints: {len(report['checkpoints'])}.", f"- Recorded rejected transitions: {len(report['failure_transitions'])}.", "- Interpolation is forbidden; solver tails and physical breakdown voltage are not extrapolated.", "", "## Fixed-state recheck", ""]
    lines.extend(f"- {row['bias_V']:.0f} V: {row['status']} — {row['reason']}" for row in report["fixed_state_recheck"])
    lines.extend(["", "## Termination", ""])
    lines.extend(f"- {row.get('solver')} {row.get('topology')} {row.get('start_bias_V')} V to {row.get('target_bias_V')} V: {row.get('incomplete_reason', 'rejected transition')}" for row in report["failure_transitions"])
    return "\n".join(lines) + "\n"


def verify_comparison_artifacts(report_path: Path) -> bool:
    """Verify every hash-addressed comparison artifact and declared input path."""
    report = json.loads(report_path.read_text(encoding="utf-8"))
    validate_bv_comparison_v1(report)
    hashes = report["artifact_hashes"]
    for name, item in report["input_artifacts"].items():
        if not isinstance(item, dict) or set(item) != {"path", "sha256"}:
            raise ValueError(f"input artifact {name} has an invalid contract")
        path = Path(item["path"])
        if not path.is_file() or _sha_bytes(path.read_bytes()) != item["sha256"]:
            raise ValueError(f"input artifact hash mismatch: {name}")
        if hashes.get(f"input:{name}") != item["sha256"]:
            raise ValueError(f"input artifact is not hash-addressed: {name}")
    for name, digest in hashes.items():
        if name.startswith("input:"):
            continue
        path = report_path.parent / name
        if not path.is_file() or _sha_bytes(path.read_bytes()) != digest:
            raise ValueError(f"generated artifact hash mismatch: {name}")
    return True


def write_comparison_package(out_dir: Path, vela_manifest: dict[str, Any], sentaurus_manifest: dict[str, Any], *, fixed_state_report: dict[str, Any], input_artifacts: dict[str, Path] | None = None) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    report = compare_sweeps(vela_manifest, sentaurus_manifest, fixed_state_report=fixed_state_report)
    csv_path = out_dir / "sweep_comparison.csv"; _write_csv(csv_path, report)
    md_path = out_dir / "sweep_comparison.md"; md_path.write_text(_markdown(report), encoding="utf-8")
    figures = _render_figures(out_dir, report)
    input_records = {name: {"path": str(path.resolve()), "sha256": _sha_bytes(path.read_bytes())} for name, path in (input_artifacts or {}).items()}
    hashes = {"sweep_comparison.csv": _sha_bytes(csv_path.read_bytes()), "sweep_comparison.md": _sha_bytes(md_path.read_bytes())}
    hashes.update({path.name: _sha_bytes(path.read_bytes()) for path in figures})
    hashes.update({f"input:{name}": item["sha256"] for name, item in input_records.items()})
    report["artifact_hashes"] = hashes
    report["input_artifacts"] = input_records
    validate_bv_comparison_v1(report)
    report_path = out_dir / "sweep_comparison.json"
    report_path.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    verify_comparison_artifacts(report_path)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vela-manifest", type=Path, required=True); parser.add_argument("--sentaurus-manifest", type=Path, required=True)
    parser.add_argument("--fixed-state-report", type=Path, required=True); parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    vela = json.loads(args.vela_manifest.read_text(encoding="utf-8")); sentaurus = json.loads(args.sentaurus_manifest.read_text(encoding="utf-8")); fixed = json.loads(args.fixed_state_report.read_text(encoding="utf-8"))
    write_comparison_package(args.out_dir.resolve(), vela, sentaurus, fixed_state_report=fixed, input_artifacts={"vela_manifest": args.vela_manifest, "sentaurus_manifest": args.sentaurus_manifest, "fixed_state_report": args.fixed_state_report})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())