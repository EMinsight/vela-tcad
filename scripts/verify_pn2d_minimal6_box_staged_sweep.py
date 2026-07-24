#!/usr/bin/env python3
"""Independent verifier for the forty-state Minimal6 box staged replay."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path


FINAL_STAGE = "sentaurus_qfp_density_element_mobility_geometry"
STAGES = (
    "vela_baseline",
    "sentaurus_qfp",
    "sentaurus_qfp_density",
    "sentaurus_qfp_density_element_mobility",
    FINAL_STAGE,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_final_plt(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    info, data = text.split("Data {", 1)
    names = re.findall(
        r'"([^"]+)"', info.split("datasets", 1)[1].split("]", 1)[0]
    )
    values = [
        float(value)
        for value in re.findall(
            r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", data
        )
    ]
    if not names or len(values) % len(names) != 0:
        raise ValueError(f"invalid PLT {path}")
    return dict(zip(names, values[-len(names) :]))


def contact_current(
    edges: dict[tuple[int, int], float], contact: set[int]
) -> float:
    outward = 0.0
    for (start, end), value in edges.items():
        if (start in contact) == (end in contact):
            continue
        outward += value if start in contact else -value
    return -outward


def divergence(edges: dict[tuple[int, int], float], node: int) -> float:
    result = 0.0
    for (start, end), value in edges.items():
        if start == node:
            result += value
        elif end == node:
            result -= value
    return result


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (
        position - lower
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    failures: list[str] = []

    for name, expected in manifest["outputs"].items():
        if sha256(root / name) != expected:
            failures.append(f"hash mismatch: {name}")

    samples = load_csv(root / "stage_edge_samples.csv")
    summaries = load_csv(root / "stage_summary.csv")
    contributions = load_csv(root / "paired_contributions.csv")
    states = load_csv(root / "state_summary.csv")
    density = load_csv(root / "density_recompute_control.csv")
    mobility = load_csv(root / "mobility_comparison.csv")
    geometry = load_csv(root / "geometry_coefficients.csv")
    terminal_recorded = load_csv(root / "terminal_closure.csv")
    kcl_recorded = load_csv(root / "total_current_kcl.csv")
    baseline = load_csv(root / "baseline_operator_crosscheck.csv")
    mapping = load_csv(root / "cell_mapping.csv")
    expected_counts = {
        "samples": (len(samples), 4320),
        "summaries": (len(summaries), 36),
        "contributions": (len(contributions), 24),
        "states": (len(states), 480),
        "density": (len(density), 240),
        "mobility": (len(mobility), 720),
        "geometry": (len(geometry), 24),
        "terminal": (len(terminal_recorded), 160),
        "kcl": (len(kcl_recorded), 80),
        "baseline": (len(baseline), 720),
        "mapping": (len(mapping), 8),
    }
    for name, (actual, expected) in expected_counts.items():
        if actual != expected:
            failures.append(f"{name} count {actual} != {expected}")

    exact_states = {
        (topology, -float(magnitude))
        for topology in ("mirror", "sketch")
        for magnitude in range(1, 21)
    }
    sample_states = {
        (row["topology"], float(row["bias_V"])) for row in samples
    }
    if sample_states != exact_states:
        failures.append("sample state lattice differs from exact 40 states")

    final = [
        row
        for row in samples
        if row["stage"] == FINAL_STAGE and row["status"] == "valid"
    ]
    if len(final) != 400:
        failures.append(f"final valid count {len(final)} != 400")
    final_relative = [
        abs(float(row["candidate_A_per_um"]) - float(row["reference_A_per_um"]))
        / abs(float(row["reference_A_per_um"]))
        for row in final
    ]
    if max(final_relative, default=math.inf) != 0.0:
        failures.append("final stage does not exactly replay its reference")

    mapping_keys = {
        (row["topology"], int(row["vela_triangle_id"])) for row in mapping
    }
    if mapping_keys != {
        (topology, triangle)
        for topology in ("mirror", "sketch")
        for triangle in range(4)
    }:
        failures.append("cell mapping is incomplete")
    mapping_max_residual = max(
        float(row["electric_field_relative_residual"]) for row in mapping
    )
    if mapping_max_residual >= 1.0e-12:
        failures.append(f"cell mapping residual too large: {mapping_max_residual}")

    baseline_max = max(float(row["relative_difference"]) for row in baseline)
    density_max = max(
        max(float(row["n_abs_dex"]), float(row["p_abs_dex"]))
        for row in density
    )
    if baseline_max >= 1.0e-12:
        failures.append(f"baseline operator mismatch: {baseline_max}")
    if density_max >= 5.0e-6:
        failures.append(f"density control mismatch: {density_max}")

    # Independently regenerate pooled summaries.
    for row in summaries:
        scope = row["scope"]
        selected = [
            sample
            for sample in samples
            if sample["stage"] == row["stage"]
            and sample["carrier"] == row["carrier"]
            and sample["status"] == "valid"
            and (scope == "all" or sample["topology"] == scope)
        ]
        errors = [float(sample["absolute_log10_error_dex"]) for sample in selected]
        signs = [float(sample["sign_agreement"]) for sample in selected]
        checks = (
            (int(row["valid_count"]), len(errors), 0.0),
            (
                float(row["median_abs_dex"]),
                statistics.median(errors),
                1.0e-14,
            ),
            (float(row["p95_abs_dex"]), quantile(errors, 0.95), 1.0e-14),
            (float(row["maximum_abs_dex"]), max(errors), 1.0e-14),
            (
                float(row["sign_agreement_fraction"]),
                statistics.mean(signs),
                1.0e-15,
            ),
        )
        if any(
            not math.isclose(float(actual), float(expected), rel_tol=1.0e-14, abs_tol=tolerance)
            for actual, expected, tolerance in checks
        ):
            failures.append(
                f"summary mismatch: {scope} {row['stage']} {row['carrier']}"
            )

    sample_index = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["stage"],
            row["carrier"],
            int(row["edge_id"]),
        ): row
        for row in samples
    }
    for row in contributions:
        paired: list[float] = []
        for topology, bias in sorted(exact_states):
            if row["scope"] != "all" and topology != row["scope"]:
                continue
            for edge in range(9):
                before = sample_index[
                    (
                        topology,
                        bias,
                        row["previous_stage"],
                        row["carrier"],
                        edge,
                    )
                ]
                after = sample_index[
                    (
                        topology,
                        bias,
                        row["current_stage"],
                        row["carrier"],
                        edge,
                    )
                ]
                if before["status"] == "valid" and after["status"] == "valid":
                    paired.append(
                        float(before["absolute_log10_error_dex"])
                        - float(after["absolute_log10_error_dex"])
                    )
        if len(paired) != int(row["paired_count"]) or not math.isclose(
            statistics.median(paired),
            float(row["median_error_reduction_dex"]),
            rel_tol=1.0e-14,
            abs_tol=1.0e-15,
        ):
            failures.append(
                f"contribution mismatch: {row['scope']} {row['carrier']} "
                f"{row['current_stage']}"
            )

    # Recompute terminal and total-current KCL closure from final edge rows.
    mesh_root = Path(manifest["inputs"]["mesh_root"])
    sentaurus_root = Path(manifest["inputs"]["sentaurus_state_root"])
    final_index: dict[
        tuple[str, float, str], dict[tuple[int, int], float]
    ] = {}
    for row in final:
        key = (row["topology"], float(row["bias_V"]), row["carrier"])
        final_index.setdefault(key, {})[
            (int(row["node0"]), int(row["node1"]))
        ] = float(row["reference_A_per_um"])
    terminal_max_carrier = 0.0
    terminal_max_total = 0.0
    kcl_max = 0.0
    for topology, bias in sorted(exact_states):
        magnitude = abs(int(bias))
        label = f"m{magnitude}V"
        mesh = json.loads(
            (mesh_root / topology / "mesh.json").read_text(encoding="utf-8")
        )
        contacts = {
            contact["name"]: {int(node) for node in contact["node_ids"]}
            for contact in mesh["contacts"]
        }
        terminal = parse_final_plt(
            sentaurus_root
            / topology
            / label
            / f"pn2d_minimal6_state_{label}.plt"
        )
        electron = final_index[(topology, bias, "electron")]
        hole = final_index[(topology, bias, "hole")]
        total_edges = {
            pair: electron[pair] + hole[pair] for pair in electron
        }
        scale = max(
            abs(terminal["Anode TotalCurrent"]),
            abs(terminal["Cathode TotalCurrent"]),
            1.0e-300,
        )
        for contact in ("Anode", "Cathode"):
            predicted_total = 0.0
            for carrier, field, edges in (
                ("electron", "eCurrent", electron),
                ("hole", "hCurrent", hole),
            ):
                predicted = contact_current(edges, contacts[contact])
                observed = terminal[f"{contact} {field}"]
                terminal_max_carrier = max(
                    terminal_max_carrier,
                    abs(predicted - observed) / max(abs(observed), 1.0e-300),
                )
                predicted_total += predicted
            terminal_max_total = max(
                terminal_max_total,
                abs(predicted_total - terminal[f"{contact} TotalCurrent"])
                / scale,
            )
        for node in (1, 5):
            kcl_max = max(
                kcl_max, abs(divergence(total_edges, node)) / scale
            )
    if terminal_max_carrier >= 3.0e-3:
        failures.append(f"carrier terminal closure failed: {terminal_max_carrier}")
    if terminal_max_total >= 1.0e-4:
        failures.append(f"total terminal closure failed: {terminal_max_total}")
    if kcl_max >= 1.0e-8:
        failures.append(f"total-current KCL failed: {kcl_max}")

    result = {
        "status": "passed" if not failures else "failed",
        "failure_count": len(failures),
        "failures": failures,
        "checked_states": len(exact_states),
        "checked_stage_samples": len(samples),
        "checked_final_carrier_edges": len(final),
        "cell_mapping_max_relative_residual": mapping_max_residual,
        "baseline_operator_max_relative_difference": baseline_max,
        "recomputed_density_max_abs_dex": density_max,
        "terminal_carrier_max_relative_error": terminal_max_carrier,
        "terminal_total_max_relative_error": terminal_max_total,
        "total_current_kcl_max_relative_error": kcl_max,
    }
    (root / "independent_verification.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
