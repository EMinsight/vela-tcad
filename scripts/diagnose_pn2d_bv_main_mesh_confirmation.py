#!/usr/bin/env python3
"""Confirm the compensated-SG mechanism on the PN2D main mesh at five BV anchors."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
import diagnose_pn2d_bv_compensated_source_proxy as compensated

DEFAULT_BIASES = (-10.0, -13.2, -18.0, -19.0, -20.0)
ROW_FIELDS = [
    "bias_V", "edge_id", "node0", "node1", "x0_um", "y0_um", "x1_um", "y1_um",
    "vela_vtk_path", "sentaurus_export_dir", "edge_type",
    "vela_source_physical_m_inv_s",
    "sentaurus_same_area_source_proxy_physical_m_inv_s", "sentaurus_source_basis",
    "vela_source_p99_m_inv_s", "sentaurus_same_area_source_proxy_p99_m_inv_s",
    "vela_active", "sentaurus_active",
    "support_class", "in_active_union", "vela_flux_abs_m2_s",
    "sentaurus_vector_flux_abs_m2_s", "sentaurus_replay_flux_abs_m2_s",
    "vela_over_sentaurus_vector_abs_ratio", "sentaurus_replay_over_vector_abs_ratio",
    "original_flux_gap_dex", "replay_flux_gap_dex", "gap_recovery",
    "vela_alpha_m_inv", "sentaurus_alpha_same_edge_m_inv",
    "vela_over_sentaurus_alpha_abs_ratio", "alpha_gap_dex",
    "vela_source_over_sentaurus_same_area_proxy_abs_ratio", "source_gap_dex",
    "production_highprec_relative_error", "cancellation_condition",
    "any_exponent_clamped", "row_mechanism_classification", "row_mechanism_rule",
]


def parse_biases(raw: str) -> tuple[float, ...]:
    values = tuple(float(item.strip()) for item in raw.split(",") if item.strip())
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError("biases must contain finite values")
    return values


def parse_sentaurus_exports(specs: list[str]) -> dict[float, Path]:
    result: dict[float, Path] = {}
    for spec in specs:
        bias_raw, separator, path_raw = spec.partition("=")
        if not separator or not path_raw:
            raise ValueError("--sentaurus-export must be BIAS=EXPORT_DIR")
        try:
            bias = round(float(bias_raw), 10)
        except ValueError as exc:
            raise ValueError(f"invalid Sentaurus export bias: {bias_raw}") from exc
        if bias in result:
            raise ValueError(f"duplicate Sentaurus export for {bias:g} V")
        result[bias] = Path(path_raw)
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sg-csv", type=Path)
    parser.add_argument("--vtk-root", type=Path)
    parser.add_argument("--imported-doping", type=Path)
    parser.add_argument("--sentaurus-root", type=Path)
    parser.add_argument(
        "--sentaurus-export", action="append", default=[], metavar="BIAS=EXPORT_DIR",
        help="Repeatable exact export mapping; supports anchors stored under different roots.",
    )
    parser.add_argument("--write-diagnostic-deck-from", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default=",".join(f"{v:g}" for v in DEFAULT_BIASES))
    parser.add_argument("--percentile", type=float, default=99.0)
    parser.add_argument("--vtk-prefix", default="dc_sweep")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    args = parser.parse_args(argv)
    args.biases = parse_biases(args.biases)
    try:
        args.sentaurus_exports = parse_sentaurus_exports(args.sentaurus_export)
    except ValueError as exc:
        parser.error(str(exc))
    if args.write_diagnostic_deck_from is None:
        missing = [
            name for name in ("sg_csv", "vtk_root", "imported_doping")
            if getattr(args, name) is None
        ]
        if missing:
            parser.error("analysis mode requires " + ", ".join(f"--{name.replace('_', '-')}" for name in missing))
        if args.sentaurus_root is None and not args.sentaurus_exports:
            parser.error(
                "analysis mode requires --sentaurus-root or repeatable --sentaurus-export"
            )
    if not 0.0 <= args.percentile <= 100.0:
        parser.error("--percentile must be in [0, 100]")
    if not math.isfinite(args.temperature_k) or args.temperature_k <= 0.0:
        parser.error("--temperature-k must be positive and finite")
    return args


def write_diagnostic_deck(template_path: Path, out_dir: Path) -> Path:
    if not template_path.exists():
        raise FileNotFoundError(template_path)
    deck = json.loads(template_path.read_text(encoding="utf-8-sig"))
    if deck.get("simulation_type") != "dc_sweep":
        raise ValueError("diagnostic deck template must use simulation_type=dc_sweep")
    sweep = deck.get("sweep")
    if not isinstance(sweep, dict) or sweep.get("mode") != "bv_reverse":
        raise ValueError("diagnostic deck template must use sweep.mode=bv_reverse")
    solver = deck.get("solver", {})
    impact = solver.get("impact_ionization", {}) if isinstance(solver, dict) else {}
    if (
        impact.get("generation") != "current_density"
        or impact.get("current_approximation") != "density_gradient"
    ):
        raise ValueError("diagnostic deck requires current_density generation and density_gradient current approximation")
    raw_points = sweep.get("bias_points", [])
    if not isinstance(raw_points, list):
        raise ValueError("diagnostic deck sweep.bias_points must be a list")
    points: list[float] = []
    for raw in raw_points:
        value = float(raw)
        if not math.isfinite(value):
            raise ValueError("diagnostic deck bias_points must be finite")
        points.append(value)
    by_key = {round(value, 10): value for value in points}
    for bias in DEFAULT_BIASES:
        by_key[round(bias, 10)] = bias
    sweep["bias_points"] = sorted(by_key.values(), reverse=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    output_root = out_dir.resolve()
    (output_root / "vtk").mkdir(parents=True, exist_ok=True)
    deck["_comment"] = (
        "PN2D current HEAD main-mesh confirmation deck generated from "
        f"{template_path.resolve()}; VTK and assembled SG-edge diagnostics enabled."
    )
    deck["output_csv"] = str(output_root / "main_mesh_confirmation_sweep.csv")
    sweep["stop"] = min(float(sweep.get("stop", -20.0)), min(DEFAULT_BIASES))
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str(output_root / "vtk" / "dc_sweep")
    sweep["write_state_file"] = str(output_root / "main_mesh_confirmation_last_state.csv")
    diagnostics = sweep.get("diagnostics", {})
    if not isinstance(diagnostics, dict):
        raise ValueError("diagnostic deck sweep.diagnostics must be an object")
    diagnostics = dict(diagnostics)
    diagnostics["sg_avalanche_edges"] = {
        "enabled": True,
        "csv_file": str(output_root / "sg_avalanche_edges.csv"),
    }
    sweep["diagnostics"] = diagnostics
    deck_path = output_root / "main_mesh_confirmation_diagnostic_deck.json"
    deck_path.write_text(
        json.dumps(deck, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return deck_path


def percentile(values: list[float], pct: float) -> float:
    clean = sorted(abs(float(v)) for v in values if math.isfinite(float(v)))
    if not clean:
        raise ValueError("percentile requires finite values")
    if pct <= 0.0:
        return clean[0]
    if pct >= 100.0:
        return clean[-1]
    position = (len(clean) - 1) * pct / 100.0
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return clean[low]
    return clean[low] * (high - position) + clean[high] * (position - low)


def support_class(sentaurus_active: bool, vela_active: bool) -> str:
    if sentaurus_active and vela_active:
        return "overlap"
    if sentaurus_active:
        return "false_negative"
    if vela_active:
        return "false_positive"
    return "inactive"


def abs_log10_ratio(numerator: float, denominator: float) -> float:
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        raise ValueError("gap operands must be finite")
    floor = 1.0e-300
    return abs(math.log10(max(abs(numerator), floor) / max(abs(denominator), floor)))


def gap_recovery(original_gap_dex: float, replay_gap_dex: float) -> float:
    if original_gap_dex <= 1.0e-15:
        return 1.0 if replay_gap_dex <= 0.1 else 0.0
    return max(-1.0, min(1.0, 1.0 - replay_gap_dex / original_gap_dex))


def _median(rows: list[dict[str, Any]], field: str) -> float:
    values = [float(row[field]) for row in rows]
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError(f"cannot compute finite median for {field}")
    return float(statistics.median(values))


def _row_evidence(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "double_highprec_relative_error": float(row["production_highprec_relative_error"]),
        "cancellation_condition": float(row["cancellation_condition"]),
        "any_exponent_clamped": bool(row["any_exponent_clamped"]),
        "raw_edge_residual_dex": float(row["original_flux_gap_dex"]),
        "sent_state_vector_residual_dex": float(row["replay_flux_gap_dex"]),
        "sent_state_replay_residual_dex": float(row["replay_flux_gap_dex"]),
        "sent_state_gap_recovery": float(row["gap_recovery"]),
        "source_residual_dex": float(row["source_gap_dex"]),
        "alpha_residual_dex": float(row["alpha_gap_dex"]),
    }


def analyze_anchor_rows(
    bias: float,
    raw_rows: list[dict[str, Any]],
    *,
    percentile_value: float = 99.0,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not raw_rows:
        raise ValueError(f"no main-mesh edge rows for {bias:g} V")
    required = (
        "edge_id", "node0", "node1", "vela_source_physical_m_inv_s",
        "sentaurus_same_area_source_proxy_physical_m_inv_s",
        "vela_alpha_m_inv", "sentaurus_alpha_same_edge_m_inv", "vela_flux_abs_m2_s",
        "sentaurus_vector_flux_abs_m2_s", "sentaurus_replay_flux_abs_m2_s",
        "production_highprec_relative_error", "cancellation_condition",
    )
    rows: list[dict[str, Any]] = []
    for index, raw in enumerate(raw_rows):
        compensated.require_finite_fields(raw, required, context=f"anchor {bias:g} row {index}")
        row = dict(raw)
        row["bias_V"] = float(bias)
        rows.append(row)

    vela_threshold = percentile(
        [float(row["vela_source_physical_m_inv_s"]) for row in rows], percentile_value
    )
    sentaurus_threshold = percentile(
        [float(row["sentaurus_same_area_source_proxy_physical_m_inv_s"]) for row in rows],
        percentile_value,
    )
    for row in rows:
        vela_source = abs(float(row["vela_source_physical_m_inv_s"]))
        sentaurus_source = abs(float(row["sentaurus_same_area_source_proxy_physical_m_inv_s"]))
        vela_active = vela_source >= vela_threshold
        sentaurus_active = sentaurus_source >= sentaurus_threshold
        row.update({
            "vela_source_p99_m_inv_s": vela_threshold,
            "sentaurus_same_area_source_proxy_p99_m_inv_s": sentaurus_threshold,
            "vela_active": int(vela_active),
            "sentaurus_active": int(sentaurus_active),
            "support_class": support_class(sentaurus_active, vela_active),
            "in_active_union": int(sentaurus_active or vela_active),
        })
        row["vela_over_sentaurus_vector_abs_ratio"] = compensated._finite_abs_ratio(
            float(row["vela_flux_abs_m2_s"]), float(row["sentaurus_vector_flux_abs_m2_s"])
        )
        row["sentaurus_replay_over_vector_abs_ratio"] = compensated._finite_abs_ratio(
            float(row["sentaurus_replay_flux_abs_m2_s"]),
            float(row["sentaurus_vector_flux_abs_m2_s"]),
        )
        row["original_flux_gap_dex"] = abs_log10_ratio(
            float(row["vela_flux_abs_m2_s"]), float(row["sentaurus_vector_flux_abs_m2_s"])
        )
        row["replay_flux_gap_dex"] = abs_log10_ratio(
            float(row["sentaurus_replay_flux_abs_m2_s"]),
            float(row["sentaurus_vector_flux_abs_m2_s"]),
        )
        row["gap_recovery"] = gap_recovery(
            float(row["original_flux_gap_dex"]), float(row["replay_flux_gap_dex"])
        )
        row["vela_over_sentaurus_alpha_abs_ratio"] = compensated._finite_abs_ratio(
            float(row["vela_alpha_m_inv"]),
            float(row["sentaurus_alpha_same_edge_m_inv"]),
        )
        row["alpha_gap_dex"] = abs_log10_ratio(
            float(row["vela_alpha_m_inv"]),
            float(row["sentaurus_alpha_same_edge_m_inv"]),
        )
        row["vela_source_over_sentaurus_same_area_proxy_abs_ratio"] = compensated._finite_abs_ratio(
            float(row["vela_source_physical_m_inv_s"]),
            float(row["sentaurus_same_area_source_proxy_physical_m_inv_s"]),
        )
        row["source_gap_dex"] = abs_log10_ratio(
            float(row["vela_source_physical_m_inv_s"]),
            float(row["sentaurus_same_area_source_proxy_physical_m_inv_s"]),
        )
        result = compensated.classify_root_cause(_row_evidence(row))
        row["row_mechanism_classification"] = result["classification"]
        row["row_mechanism_rule"] = result["rule"]

    union = [row for row in rows if row["in_active_union"]]
    if not union:
        raise ValueError(f"p{percentile_value:g} active-support union is empty")
    evidence = {
        "double_highprec_relative_error": _median(union, "production_highprec_relative_error"),
        "cancellation_condition": _median(union, "cancellation_condition"),
        "any_exponent_clamped": any(bool(row.get("any_exponent_clamped")) for row in union),
        "raw_edge_residual_dex": _median(union, "original_flux_gap_dex"),
        "sent_state_vector_residual_dex": _median(union, "replay_flux_gap_dex"),
        "sent_state_replay_residual_dex": _median(union, "replay_flux_gap_dex"),
        "sent_state_gap_recovery": _median(union, "gap_recovery"),
        "source_residual_dex": _median(union, "source_gap_dex"),
        "alpha_residual_dex": _median(union, "alpha_gap_dex"),
    }
    row_mechanism_counts = Counter(
        str(row["row_mechanism_classification"]) for row in union
    )
    eligible = {
        name: count for name, count in row_mechanism_counts.items()
        if name not in {"inconclusive", "coarse_artifact"}
    }
    winner = (
        sorted(eligible, key=lambda name: (-eligible[name], name))[0]
        if eligible else "inconclusive"
    )
    winner_count = eligible.get(winner, 0)
    if winner_count * 2 <= len(union):
        winner = "inconclusive"
        winner_count = row_mechanism_counts.get("inconclusive", 0)
    mechanism = {
        "classification": winner,
        "rule": "strict majority of per-edge classifier results on the p99 active-support union",
        "evidence": evidence,
        "row_counts": dict(sorted(row_mechanism_counts.items())),
        "winning_row_count": winner_count,
        "union_row_count": len(union),
    }
    counts = Counter(row["support_class"] for row in union)
    return rows, {
        "bias_V": float(bias),
        "percentile": float(percentile_value),
        "edge_count": len(rows),
        "vela_source_p99_m_inv_s": vela_threshold,
        "sentaurus_same_area_source_proxy_p99_m_inv_s": sentaurus_threshold,
        "union_count": len(union),
        "overlap_count": counts.get("overlap", 0),
        "false_positive_count": counts.get("false_positive", 0),
        "false_negative_count": counts.get("false_negative", 0),
        "jaccard": counts.get("overlap", 0) / len(union),
        "median_original_flux_gap_dex": evidence["raw_edge_residual_dex"],
        "median_replay_flux_gap_dex": evidence["sent_state_replay_residual_dex"],
        "median_gap_recovery": evidence["sent_state_gap_recovery"],
        "median_source_gap_dex": evidence["source_residual_dex"],
        "median_alpha_gap_dex": evidence["alpha_residual_dex"],
        "mechanism": mechanism,
    }


def evaluate_confirmation_gate(anchors: list[dict[str, Any]]) -> dict[str, Any]:
    expected = {round(value, 10) for value in DEFAULT_BIASES}
    actual = {round(float(anchor["bias_V"]), 10) for anchor in anchors}
    names = [
        str(anchor.get("mechanism", {}).get("classification", "inconclusive"))
        for anchor in anchors
    ]
    counts = Counter(name for name in names if name not in {"inconclusive", "coarse_artifact"})
    dominant = (
        sorted(counts, key=lambda name: (-counts[name], name))[0]
        if counts else "inconclusive"
    )
    same_count = counts.get(dominant, 0)
    coverage_pass = actual == expected and len(anchors) == len(DEFAULT_BIASES)
    same_pass = coverage_pass and same_count >= 4
    false_positive_total = sum(int(anchor.get("false_positive_count", 0)) for anchor in anchors)
    false_negative_total = sum(int(anchor.get("false_negative_count", 0)) for anchor in anchors)
    support_rule = (
        "Across the complete five-anchor set, both false-positive and false-negative "
        "p99 support counts must be nonzero."
    )
    support_pass = coverage_pass and false_positive_total > 0 and false_negative_total > 0
    by_bias = {round(float(anchor["bias_V"]), 10): anchor for anchor in anchors}
    high: dict[str, float | None] = {}
    high_pass = True
    for bias in (-19.0, -20.0):
        anchor = by_bias.get(round(bias, 10))
        recovery = float(anchor["median_gap_recovery"]) if anchor else None
        high[f"{bias:g}"] = recovery
        high_pass = high_pass and recovery is not None and recovery >= 0.8
    status = (
        "pass" if coverage_pass and support_pass and same_pass and high_pass else "fail"
    )
    if not coverage_pass:
        target = "main_mesh_anchor_coverage"
        minimum_test = "test_pn2d_bv_main_mesh_confirmation_requires_five_anchors"
    elif not support_pass:
        target = "main_mesh_bidirectional_support_explanation"
        minimum_test = "test_pn2d_bv_main_mesh_confirmation_requires_bidirectional_support"
    elif not same_pass:
        target = "main_mesh_anchor_mechanism_consistency"
        minimum_test = "test_pn2d_bv_main_mesh_confirmation_same_mechanism_gate"
    elif not high_pass:
        target = "high_bias_sent_state_replay_recovery"
        minimum_test = "test_pn2d_bv_main_mesh_confirmation_high_bias_recovery_gate"
    else:
        target = "main_mesh_continuation_branch_recovery"
        minimum_test = "test_pn2d_bv_main_mesh_continuation_recovers_multiplication_current"
    return {
        "status": status,
        "artifact_contract_pass": True,
        "coverage_pass": coverage_pass,
        "support_bidirectional_pass": support_pass,
        "support_bidirectional_rule": support_rule,
        "support_false_positive_total": false_positive_total,
        "support_false_negative_total": false_negative_total,
        "same_mechanism_pass": same_pass,
        "high_bias_recovery_pass": high_pass,
        "dominant_mechanism": dominant,
        "mechanism_counts": dict(sorted(counts.items())),
        "same_mechanism_count": same_count,
        "required_same_mechanism_count": 4,
        "high_bias_recovery": high,
        "required_high_bias_recovery": 0.8,
        "next_target": target,
        "minimum_failing_test": minimum_test,
    }


def resolve_vtk_for_bias(root: Path, prefix: str, bias: float) -> Path:
    pattern = re.compile(
        rf"^{re.escape(Path(prefix).name)}_\d+_([-+0-9.eEpP]+)V\.vtk$",
        re.IGNORECASE,
    )
    matches: list[Path] = []
    for path in sorted(root.glob(f"{prefix}_*.vtk")):
        match = pattern.match(path.name)
        if match is None:
            continue
        try:
            file_bias = float(match.group(1).replace("p", ".").replace("P", "."))
        except ValueError:
            continue
        if abs(file_bias - bias) <= 1.0e-8:
            matches.append(path)
    if not matches:
        raise FileNotFoundError(f"no VTK filename encodes bias {bias:g} V in {root}")
    if len(matches) != 1:
        raise ValueError(f"ambiguous VTK files for {bias:g} V: {matches}")
    return matches[0]


def resolve_sentaurus_export_dir(args: argparse.Namespace, bias: float) -> Path:
    key = round(bias, 10)
    if key in args.sentaurus_exports:
        return args.sentaurus_exports[key]
    if args.sentaurus_root is None:
        raise FileNotFoundError(f"no exact Sentaurus export mapping for {bias:g} V")
    dotted = args.sentaurus_root / f"sentaurus_{bias:g}v"
    compact = args.sentaurus_root / f"sentaurus_{f'{bias:g}'.replace('.', 'p')}v"
    for candidate in (dotted, compact):
        if candidate.exists():
            return candidate
    return dotted


def load_sentaurus_nodes_from_export(export_dir: Path) -> list[dict[str, Any]]:
    rows = compensated.read_csv(export_dir / "nodes.csv")
    if not rows:
        raise ValueError(f"empty Sentaurus nodes.csv in {export_dir}")
    return [
        {
            "id": int(row["id"]),
            "x_um": compensated.finite_float(row["x_um"]),
            "y_um": compensated.finite_float(row["y_um"]),
        }
        for row in rows
    ]


def inspect_sentaurus_export(export_dir: Path, bias: float) -> dict[str, Any]:
    reasons: list[str] = []
    manifest_path = export_dir / "field_manifest.json"
    manifest: dict[str, Any] | None = None
    if not export_dir.exists():
        reasons.append(f"missing export directory: {export_dir}")
    elif not manifest_path.exists():
        reasons.append(f"missing field_manifest.json: {manifest_path}")
    else:
        try:
            loaded = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
            if not isinstance(loaded, dict):
                raise ValueError("manifest root must be an object")
            manifest = loaded
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            reasons.append(f"invalid field_manifest.json: {exc}")
    if manifest is not None:
        try:
            compensated.validate_manifest_vector_components(
                manifest, "eCurrentDensity", expected_components=2
            )
        except ValueError as exc:
            reasons.append(f"eCurrentDensity components=2 required: {exc}")
        try:
            compensated._require_manifest_field(
                manifest, "eAlphaAvalanche", components=1, unit="cm^-1"
            )
        except ValueError as exc:
            reasons.append(f"eAlphaAvalanche raw export required: {exc}")
    if not reasons:
        try:
            load_sentaurus_nodes_from_export(export_dir)
            compensated.load_sentaurus_electron_state(export_dir)
        except (OSError, ValueError, KeyError) as exc:
            reasons.append(f"strict vector-current/alpha loader failed: {exc}")
    return {
        "bias_V": float(bias),
        "export_dir": str(export_dir),
        "status": "pass" if not reasons else "fail",
        "reasons": reasons,
        "required_raw_fields": ["eCurrentDensity components=2", "eAlphaAvalanche components=1"],
        "fallback_policy": "scalar current magnitude is forbidden as a vector projection substitute",
    }


def preflight_sentaurus_exports(args: argparse.Namespace) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    for bias in args.biases:
        try:
            export_dir = resolve_sentaurus_export_dir(args, bias)
            checks.append(inspect_sentaurus_export(export_dir, bias))
        except (OSError, ValueError) as exc:
            checks.append({
                "bias_V": float(bias),
                "export_dir": "",
                "status": "fail",
                "reasons": [str(exc)],
                "required_raw_fields": [
                    "eCurrentDensity components=2", "eAlphaAvalanche components=1"
                ],
                "fallback_policy": (
                    "scalar current magnitude is forbidden as a vector projection substitute"
                ),
            })
    return checks


def _sentaurus_endpoint(
    nodes: list[dict[str, Any]],
    by_id: dict[int, dict[str, Any]],
    edge: dict[str, str],
    endpoint: int,
) -> dict[str, Any]:
    node_id = int(float(edge[f"node{endpoint}"]))
    direct = by_id.get(node_id)
    x_um = float(edge[f"x{endpoint}_um"])
    y_um = float(edge[f"y{endpoint}_um"])
    if direct is not None and math.hypot(
        float(direct["x_um"]) - x_um, float(direct["y_um"]) - y_um
    ) <= 1.0e-6:
        return direct
    return compensated.nearest_sentaurus_node(nodes, x_um, y_um)


def collect_anchor_edges(args: argparse.Namespace, bias: float) -> list[dict[str, Any]]:
    doping = compensated.load_doping(args.imported_doping)
    vtk_path = resolve_vtk_for_bias(args.vtk_root, args.vtk_prefix, bias)
    vela = compensated.parse_vtk(vtk_path, coordinate_scale_to_um=1.0)
    export_dir = resolve_sentaurus_export_dir(args, bias)
    sent_nodes = load_sentaurus_nodes_from_export(export_dir)
    sent_by_id = {int(node["id"]): node for node in sent_nodes}
    sent_state = compensated.load_sentaurus_electron_state(export_dir)
    selected = [
        row for row in compensated.read_csv(args.sg_csv)
        if abs(float(row.get("bias_V", "nan")) - bias) <= 1.0e-8
    ]
    if not selected:
        raise ValueError(f"SG CSV has no rows for {bias:g} V")
    rows: list[dict[str, Any]] = []
    for edge in selected:
        node0 = int(float(edge["node0"]))
        node1 = int(float(edge["node1"]))
        if node0 not in doping or node1 not in doping:
            raise ValueError(f"edge {edge.get('edge_id')} has missing imported doping")
        if node0 >= len(vela["points"]) or node1 >= len(vela["points"]):
            raise ValueError(f"edge {edge.get('edge_id')} is outside VTK points")
        sent0 = _sentaurus_endpoint(sent_nodes, sent_by_id, edge, 0)
        sent1 = _sentaurus_endpoint(sent_nodes, sent_by_id, edge, 1)
        enriched = compensated.enrich_edge_with_sentaurus_replay(
            edge_row=edge, sentaurus_state=sent_state,
            sentaurus_node0=sent0, sentaurus_node1=sent1,
            temperature_K=args.temperature_k, unit_system="tcad_internal",
        )
        clamp = any(
            bool(int(float(edge[name])))
            for name in (
                "electron_sg_node0_exponent_clamped_low",
                "electron_sg_node0_exponent_clamped_high",
                "electron_sg_node1_exponent_clamped_low",
                "electron_sg_node1_exponent_clamped_high",
            )
        )
        p0 = vela["points"][node0]
        p1 = vela["points"][node1]
        rows.append({
            "edge_id": int(float(edge["edge_id"])), "node0": node0, "node1": node1,
            "x0_um": float(p0[0]), "y0_um": float(p0[1]),
            "x1_um": float(p1[0]), "y1_um": float(p1[1]),
            "vela_vtk_path": str(vtk_path.resolve()),
            "sentaurus_export_dir": str(export_dir.resolve()),
            "edge_type": f"{doping[node0]['type']}-{doping[node1]['type']}",
            "vela_source_physical_m_inv_s": float(
                enriched["vela_e_source_integral_physical_m_inv_s"]
            ),
            "sentaurus_same_area_source_proxy_physical_m_inv_s": float(
                enriched["sentaurus_e_source_on_vela_area_physical_m_inv_s"]
            ),
            "sentaurus_source_basis": "Sentaurus alpha average x vector flux x Vela edge area",
            "vela_alpha_m_inv": float(edge["electron_alpha_m_inv"]),
            "sentaurus_alpha_same_edge_m_inv": float(enriched["sentaurus_e_alpha_edge_average_m_inv"]),
            "vela_flux_abs_m2_s": abs(float(
                enriched["vela_e_sg_production_canonical_signed_flux_m2_s"]
            )),
            "sentaurus_vector_flux_abs_m2_s": abs(float(
                enriched["sentaurus_e_continuity_edge_signed_flux_m2_s"]
            )),
            "sentaurus_replay_flux_abs_m2_s": abs(float(
                enriched["sentaurus_e_sg_vela_mobility_signed_flux_m2_s"]
            )),
            "production_highprec_relative_error": float(
                edge["electron_sg_production_vs_high_precision_reference_relative_error"]
            ),
            "cancellation_condition": float(edge["electron_sg_cancellation_condition"]),
            "any_exponent_clamped": clamp,
        })
    return rows


def run_confirmation(args: argparse.Namespace) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    artifact_checks = preflight_sentaurus_exports(args)
    inputs = {
        "sg_csv": str(args.sg_csv),
        "vtk_root": str(args.vtk_root),
        "imported_doping": str(args.imported_doping),
        "sentaurus_root": str(args.sentaurus_root) if args.sentaurus_root else None,
        "sentaurus_exports": {
            f"{bias:g}": str(path) for bias, path in sorted(args.sentaurus_exports.items())
        },
        "vtk_prefix": args.vtk_prefix,
    }
    failed_checks = [check for check in artifact_checks if check["status"] != "pass"]
    if failed_checks:
        gate = {
            "status": "fail",
            "artifact_contract_pass": False,
            "coverage_pass": False,
            "support_bidirectional_pass": False,
            "support_bidirectional_rule": (
                "Across the complete five-anchor set, both false-positive and false-negative "
                "p99 support counts must be nonzero."
            ),
            "support_false_positive_total": 0,
            "support_false_negative_total": 0,
            "same_mechanism_pass": False,
            "high_bias_recovery_pass": False,
            "dominant_mechanism": "inconclusive",
            "mechanism_counts": {},
            "same_mechanism_count": 0,
            "required_same_mechanism_count": 4,
            "high_bias_recovery": {"-19": None, "-20": None},
            "required_high_bias_recovery": 0.8,
            "next_target": "sentaurus_main_mesh_vector_current_alpha_export",
            "minimum_failing_test": (
                "test_pn2d_bv_main_mesh_confirmation_requires_vector_current_and_alpha_exports"
            ),
            "missing_artifacts": failed_checks,
        }
        return [], {
            "schema": "vela.pn2d_bv_main_mesh_confirmation.v1",
            "inputs": inputs,
            "biases": list(args.biases),
            "percentile": args.percentile,
            "row_count": 0,
            "anchors": [],
            "artifact_checks": artifact_checks,
            "gate": gate,
        }
    all_rows: list[dict[str, Any]] = []
    anchors: list[dict[str, Any]] = []
    for bias in args.biases:
        rows, summary = analyze_anchor_rows(
            bias, collect_anchor_edges(args, bias), percentile_value=args.percentile
        )
        all_rows.extend(rows)
        anchors.append(summary)
    payload = {
        "schema": "vela.pn2d_bv_main_mesh_confirmation.v1",
        "inputs": inputs,
        "biases": list(args.biases), "percentile": args.percentile,
        "row_count": len(all_rows), "anchors": anchors,
        "artifact_checks": artifact_checks,
        "gate": evaluate_confirmation_gate(anchors),
    }
    return all_rows, payload


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in ROW_FIELDS})


def write_report(path: Path, payload: dict[str, Any]) -> None:
    gate = payload["gate"]
    lines = [
        "# PN2D BV Main-Mesh Mechanism Confirmation",
        "",
        (
            "Electron SG evidence on the p99 main-mesh active-support union. The "
            "Sentaurus source is a same-area proxy (Sentaurus alpha average times "
            "projected vector flux times Vela edge area), not a native Sentaurus source field."
        ),
        "No solver behavior is modified.",
        "",
        "## Raw Artifact Contract",
        "",
        f"- Contract pass: {gate['artifact_contract_pass']}",
        "- Required: two-component eCurrentDensity and scalar eAlphaAvalanche at every anchor.",
        "- Scalar current magnitude is never substituted for vector projection.",
    ]
    for check in payload.get("artifact_checks", []):
        details = "; ".join(check.get("reasons", [])) or "strict export contract satisfied"
        lines.append(
            f"- {float(check['bias_V']):g} V: {check['status']} - {details}"
        )
    lines.extend(["", "## Anchor Evidence", ""])
    if payload["anchors"]:
        lines.extend([
            "| bias V | overlap | false positive | false negative | original gap dex | replay gap dex | source gap dex | alpha gap dex | recovery | mechanism |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ])
        for anchor in payload["anchors"]:
            lines.append(
                "| {bias:g} | {overlap} | {fp} | {fn} | {original:.6g} | {replay:.6g} | "
                "{source:.6g} | {alpha:.6g} | {recovery:.6g} | {mechanism} |".format(
                    bias=anchor["bias_V"], overlap=anchor["overlap_count"],
                    fp=anchor["false_positive_count"], fn=anchor["false_negative_count"],
                    original=anchor["median_original_flux_gap_dex"],
                    replay=anchor["median_replay_flux_gap_dex"],
                    source=anchor["median_source_gap_dex"],
                    alpha=anchor["median_alpha_gap_dex"],
                    recovery=anchor["median_gap_recovery"],
                    mechanism=anchor["mechanism"]["classification"],
                )
            )
    else:
        lines.append("Not evaluated because the strict raw-artifact contract failed.")
    lines.extend([
        "", "## Gate", "", f"- Status: {gate['status']}",
        f"- Dominant mechanism: {gate['dominant_mechanism']}",
        f"- Same-mechanism anchors: {gate['same_mechanism_count']}/5 (required 4/5)",
        (
            "- Bidirectional support: "
            f"{gate['support_bidirectional_pass']} "
            f"(false positive={gate['support_false_positive_total']}, "
            f"false negative={gate['support_false_negative_total']})"
        ),
        f"- Bidirectional support rule: {gate['support_bidirectional_rule']}",
        f"- High-bias recovery: {gate['high_bias_recovery']} (required >=80%)",
        "", "## Follow-up", "",
        f"- Unique next target: {gate['next_target']}",
        f"- Minimum failing test: {gate['minimum_failing_test']}", "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.write_diagnostic_deck_from is not None:
        deck_path = write_diagnostic_deck(args.write_diagnostic_deck_from, args.out_dir)
        print(json.dumps({
            "mode": "write_diagnostic_deck", "deck": str(deck_path),
            "required_anchors_V": list(DEFAULT_BIASES),
        }, sort_keys=True))
        return 0
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows, payload = run_confirmation(args)
    csv_path = args.out_dir / "main_mesh_confirmation_edges.csv"
    json_path = args.out_dir / "main_mesh_confirmation_summary.json"
    report_path = args.out_dir / "main_mesh_confirmation_report.md"
    write_csv(csv_path, rows)
    json_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    write_report(report_path, payload)
    print(json.dumps({
        "csv": str(csv_path), "json": str(json_path), "report": str(report_path),
        "status": payload["gate"]["status"],
        "next_target": payload["gate"]["next_target"],
        "minimum_failing_test": payload["gate"]["minimum_failing_test"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

