#!/usr/bin/env python3
"""Compare explicit Sentaurus avalanche driving-force control branches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_minimal6_element_avalanche_replay import (
    currentplot_targets,
    parse_log,
    parse_plt,
)


REFERENCE_VARIANT = "explicit_grad_qf"
CANDIDATE_VARIANTS = (
    "implicit_default",
    "explicit_electric_field",
    "grad_qf_aval_dens_grad_qf",
)
TOPOLOGIES = ("mirror", "sketch")
GROUP_IDENTITIES = {
    "vertices": ("bias_V", "vertex"),
    "elements": ("bias_V", "element"),
    "measures": ("bias_V", "element", "local_vertex", "vertex"),
    "edges": (
        "bias_V",
        "element",
        "local_edge",
        "edge",
        "start",
        "end",
    ),
    "integrals": ("bias_V",),
}
VECTOR_FIELDS = {
    "electric_field": ("efield_x_V_cm", "efield_y_V_cm"),
    "electron_qf_gradient": (
        "grad_qf_n_x_V_cm",
        "grad_qf_n_y_V_cm",
    ),
    "hole_qf_gradient": (
        "grad_qf_p_x_V_cm",
        "grad_qf_p_y_V_cm",
    ),
    "electron_element_current": (
        "current_n_x_A_cm2",
        "current_n_y_A_cm2",
    ),
    "hole_element_current": (
        "current_p_x_A_cm2",
        "current_p_y_A_cm2",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def abs_dex(candidate: float, reference: float) -> float | None:
    candidate_magnitude = abs(candidate)
    reference_magnitude = abs(reference)
    if candidate_magnitude == 0.0 or reference_magnitude == 0.0:
        return None
    return abs(math.log10(candidate_magnitude / reference_magnitude))


def symmetric_relative_error(candidate: float, reference: float) -> float:
    return abs(candidate - reference) / max(
        abs(candidate),
        abs(reference),
        1.0e-300,
    )


def vector_relative_error(
    candidate: tuple[float, float],
    reference: tuple[float, float],
) -> float:
    return math.hypot(
        candidate[0] - reference[0],
        candidate[1] - reference[1],
    ) / max(math.hypot(*candidate), math.hypot(*reference), 1.0e-300)


def vector_angle_error_deg(
    candidate: tuple[float, float],
    reference: tuple[float, float],
) -> float | None:
    if candidate == reference:
        return 0.0
    candidate_norm = math.hypot(*candidate)
    reference_norm = math.hypot(*reference)
    if candidate_norm == 0.0 or reference_norm == 0.0:
        return None
    cosine = (
        candidate[0] * reference[0] + candidate[1] * reference[1]
    ) / (candidate_norm * reference_norm)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def keyed(
    rows: Iterable[dict[str, Any]],
    identity: tuple[str, ...],
) -> dict[tuple[Any, ...], dict[str, Any]]:
    result = {}
    for row in rows:
        key = tuple(row[name] for name in identity)
        if key in result:
            raise ValueError(f"duplicate row identity {identity}={key}")
        result[key] = row
    return result


def compare_group(
    *,
    topology: str,
    candidate_variant: str,
    group: str,
    reference_rows: list[dict[str, Any]],
    candidate_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    identity = GROUP_IDENTITIES[group]
    reference = keyed(reference_rows, identity)
    candidate = keyed(candidate_rows, identity)
    if set(reference) != set(candidate):
        raise ValueError(
            f"{topology}/{candidate_variant}/{group}: identity mismatch"
        )
    rows = []
    for key in sorted(reference):
        reference_row = reference[key]
        candidate_row = candidate[key]
        fields = sorted(set(reference_row) - set(identity))
        if fields != sorted(set(candidate_row) - set(identity)):
            raise ValueError(
                f"{topology}/{candidate_variant}/{group}/{key}: field mismatch"
            )
        identities = dict(zip(identity, key, strict=True))
        for field in fields:
            reference_value = float(reference_row[field])
            candidate_value = float(candidate_row[field])
            rows.append(
                {
                    "topology": topology,
                    "bias_V": float(reference_row["bias_V"]),
                    "candidate_variant": candidate_variant,
                    "reference_variant": REFERENCE_VARIANT,
                    "group": group,
                    "identity": ";".join(
                        f"{name}={identities[name]}" for name in identity
                    ),
                    "field": field,
                    "reference_value": reference_value,
                    "candidate_value": candidate_value,
                    "delta": candidate_value - reference_value,
                    "symmetric_relative_error": symmetric_relative_error(
                        candidate_value,
                        reference_value,
                    ),
                    "absolute_log10_ratio_dex": abs_dex(
                        candidate_value,
                        reference_value,
                    ),
                    "exact": int(candidate_value == reference_value),
                }
            )
    return rows


def currentplot_rows(
    path: Path,
    target_biases: tuple[float, ...],
) -> list[dict[str, float]]:
    names, _ = parse_plt(path)
    ignored = {"time", "runtime_element_avalanche_probe"}
    selected = currentplot_targets(path, target_biases)
    return [
        {
            "bias_V": float(row["bias_V"]),
            **{
                name: float(row[name])
                for name in names
                if name not in ignored
            },
        }
        for row in selected
    ]


def compare_currentplot(
    *,
    topology: str,
    candidate_variant: str,
    reference_rows: list[dict[str, float]],
    candidate_rows: list[dict[str, float]],
) -> list[dict[str, Any]]:
    reference = keyed(reference_rows, ("bias_V",))
    candidate = keyed(candidate_rows, ("bias_V",))
    if set(reference) != set(candidate):
        raise ValueError(
            f"{topology}/{candidate_variant}/currentplot: bias mismatch"
        )
    rows = []
    for key in sorted(reference):
        reference_row = reference[key]
        candidate_row = candidate[key]
        fields = sorted(set(reference_row) - {"bias_V"})
        if fields != sorted(set(candidate_row) - {"bias_V"}):
            raise ValueError(
                f"{topology}/{candidate_variant}/currentplot: field mismatch"
            )
        for field in fields:
            reference_value = reference_row[field]
            candidate_value = candidate_row[field]
            rows.append(
                {
                    "topology": topology,
                    "bias_V": key[0],
                    "candidate_variant": candidate_variant,
                    "reference_variant": REFERENCE_VARIANT,
                    "group": "currentplot",
                    "identity": f"bias_V={key[0]}",
                    "field": field,
                    "reference_value": reference_value,
                    "candidate_value": candidate_value,
                    "delta": candidate_value - reference_value,
                    "symmetric_relative_error": symmetric_relative_error(
                        candidate_value,
                        reference_value,
                    ),
                    "absolute_log10_ratio_dex": abs_dex(
                        candidate_value,
                        reference_value,
                    ),
                    "exact": int(candidate_value == reference_value),
                }
            )
    return rows


def maximum(
    rows: Iterable[dict[str, Any]],
    metric: str,
    *,
    group: str | None = None,
    fields: set[str] | None = None,
) -> float:
    values = []
    for row in rows:
        if group is not None and row["group"] != group:
            continue
        if fields is not None and row["field"] not in fields:
            continue
        value = row[metric]
        if value is not None:
            values.append(abs(float(value)))
    return max(values, default=0.0)


def vector_summary(
    reference_elements: list[dict[str, Any]],
    candidate_elements: list[dict[str, Any]],
) -> dict[str, float]:
    reference = keyed(reference_elements, ("bias_V", "element"))
    candidate = keyed(candidate_elements, ("bias_V", "element"))
    result = {}
    for label, (x_field, y_field) in VECTOR_FIELDS.items():
        relative = []
        angles = []
        for key, reference_row in reference.items():
            candidate_row = candidate[key]
            reference_vector = (
                float(reference_row[x_field]),
                float(reference_row[y_field]),
            )
            candidate_vector = (
                float(candidate_row[x_field]),
                float(candidate_row[y_field]),
            )
            relative.append(
                vector_relative_error(candidate_vector, reference_vector)
            )
            angle = vector_angle_error_deg(candidate_vector, reference_vector)
            if angle is not None:
                angles.append(angle)
        result[f"max_{label}_vector_relative_error"] = max(
            relative,
            default=0.0,
        )
        result[f"max_{label}_angle_error_deg"] = max(angles, default=0.0)
    return result


def summary_row(
    *,
    topology: str,
    bias: float,
    candidate_variant: str,
    long_rows: list[dict[str, Any]],
    reference_elements: list[dict[str, Any]],
    candidate_elements: list[dict[str, Any]],
) -> dict[str, Any]:
    rows = [
        row
        for row in long_rows
        if row["topology"] == topology
        and row["bias_V"] == bias
        and row["candidate_variant"] == candidate_variant
    ]
    terminal_fields = {
        row["field"]
        for row in rows
        if row["group"] == "currentplot"
        and (
            row["field"].endswith("eCurrent")
            or row["field"].endswith("hCurrent")
            or row["field"].endswith("TotalCurrent")
        )
    }
    integral_fields = {
        row["field"]
        for row in rows
        if row["group"] == "currentplot"
        and "AvalancheIntegral" in row["field"]
    }
    result = {
        "topology": topology,
        "bias_V": bias,
        "candidate_variant": candidate_variant,
        "reference_variant": REFERENCE_VARIANT,
        "all_parsed_values_exact": int(all(row["exact"] for row in rows)),
        "max_potential_absolute_delta_V": maximum(
            rows,
            "delta",
            group="vertices",
            fields={"psi_V", "eQFP_V", "hQFP_V"},
        ),
        "max_density_error_dex": maximum(
            rows,
            "absolute_log10_ratio_dex",
            group="vertices",
            fields={"n_cm3", "p_cm3"},
        ),
        "max_mobility_error_dex": maximum(
            rows,
            "absolute_log10_ratio_dex",
            group="elements",
            fields={"mu_n_cm2_Vs", "mu_p_cm2_Vs"},
        ),
        "max_alpha_error_dex": maximum(
            rows,
            "absolute_log10_ratio_dex",
            group="vertices",
            fields={"alpha_n_cm_inv", "alpha_p_cm_inv"},
        ),
        "max_generation_error_dex": maximum(
            rows,
            "absolute_log10_ratio_dex",
            group="vertices",
            fields={
                "generation_n_cm3_s",
                "generation_p_cm3_s",
                "generation_total_cm3_s",
            },
        ),
        "max_edge_sg_current_error_dex": maximum(
            rows,
            "absolute_log10_ratio_dex",
            group="edges",
            fields={"sg_jn_A_cm2", "sg_jp_A_cm2"},
        ),
        "max_element_vertex_qg_error_dex": maximum(
            rows,
            "absolute_log10_ratio_dex",
            group="measures",
            fields={"qg_n_A_um", "qg_p_A_um", "qg_total_A_um"},
        ),
        "max_runtime_integrated_qg_error_dex": maximum(
            rows,
            "absolute_log10_ratio_dex",
            group="integrals",
            fields={"qg_n_A_um", "qg_p_A_um", "qg_total_A_um"},
        ),
        "max_currentplot_integral_error_dex": maximum(
            rows,
            "absolute_log10_ratio_dex",
            group="currentplot",
            fields=integral_fields,
        ),
        "max_terminal_current_relative_error": maximum(
            rows,
            "symmetric_relative_error",
            group="currentplot",
            fields=terminal_fields,
        ),
    }
    reference_at_bias = [
        row for row in reference_elements if float(row["bias_V"]) == bias
    ]
    candidate_at_bias = [
        row for row in candidate_elements if float(row["bias_V"]) == bias
    ]
    result.update(vector_summary(reference_at_bias, candidate_at_bias))
    return result


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        format(value, ".17g")
                        if isinstance(value, float)
                        else "" if value is None else value
                    )
                    for key, value in row.items()
                }
            )


def run(raw_root: Path, output: Path) -> dict[str, Any]:
    raw_root = raw_root.resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_manifest_path = raw_root / "manifest.json"
    source_manifest = json.loads(
        source_manifest_path.read_text(encoding="ascii")
    )
    if source_manifest["status"] != "passed":
        raise ValueError("source experiment manifest is not passed")
    target_biases = tuple(float(value) for value in source_manifest["biases_V"])
    input_hashes = {
        "manifest.json": sha256(source_manifest_path),
    }
    all_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []

    for topology in TOPOLOGIES:
        variant_groups = {}
        variant_currentplot = {}
        for variant in (REFERENCE_VARIANT, *CANDIDATE_VARIANTS):
            fetched = raw_root / topology / variant / "fetched"
            log_path = fetched / f"run_{variant}.out"
            plt_path = (
                fetched
                / f"runtime_element_avalanche_probe_{variant}.plt"
            )
            variant_groups[variant] = parse_log(log_path, target_biases)
            variant_currentplot[variant] = currentplot_rows(
                plt_path,
                target_biases,
            )
            input_hashes[
                log_path.relative_to(raw_root).as_posix()
            ] = sha256(log_path)
            input_hashes[
                plt_path.relative_to(raw_root).as_posix()
            ] = sha256(plt_path)

        for candidate_variant in CANDIDATE_VARIANTS:
            candidate_rows = []
            for group in GROUP_IDENTITIES:
                candidate_rows.extend(
                    compare_group(
                        topology=topology,
                        candidate_variant=candidate_variant,
                        group=group,
                        reference_rows=variant_groups[REFERENCE_VARIANT][group],
                        candidate_rows=variant_groups[candidate_variant][group],
                    )
                )
            candidate_rows.extend(
                compare_currentplot(
                    topology=topology,
                    candidate_variant=candidate_variant,
                    reference_rows=variant_currentplot[REFERENCE_VARIANT],
                    candidate_rows=variant_currentplot[candidate_variant],
                )
            )
            all_rows.extend(candidate_rows)
            for bias in target_biases:
                summaries.append(
                    summary_row(
                        topology=topology,
                        bias=bias,
                        candidate_variant=candidate_variant,
                        long_rows=candidate_rows,
                        reference_elements=variant_groups[
                            REFERENCE_VARIANT
                        ]["elements"],
                        candidate_elements=variant_groups[
                            candidate_variant
                        ]["elements"],
                    )
                )

    long_path = output / "quantity_comparison.csv"
    summary_path = output / "state_summary.csv"
    write_csv(long_path, all_rows)
    write_csv(summary_path, summaries)
    implicit = [
        row
        for row in summaries
        if row["candidate_variant"] == "implicit_default"
    ]
    electric = [
        row
        for row in summaries
        if row["candidate_variant"] == "explicit_electric_field"
    ]
    aval_dens = [
        row
        for row in summaries
        if row["candidate_variant"] == "grad_qf_aval_dens_grad_qf"
    ]
    manifest = {
        "schema_version": 1,
        "status": "valid_sentaurus_avalanche_drive_comparison",
        "experiment": "pn2d_minimal6_sentaurus_avalanche_drive_controls",
        "reference_variant": REFERENCE_VARIANT,
        "candidate_variants": list(CANDIDATE_VARIANTS),
        "topologies": list(TOPOLOGIES),
        "biases_V": list(target_biases),
        "state_count": len(summaries),
        "quantity_comparison_count": len(all_rows),
        "implicit_default_exact_match": all(
            row["all_parsed_values_exact"] == 1 for row in implicit
        ),
        "explicit_electric_field_distinct": any(
            row["max_alpha_error_dex"] > 0.0
            or row["max_generation_error_dex"] > 0.0
            for row in electric
        ),
        "aval_dens_grad_qf_distinct": any(
            row["max_generation_error_dex"] > 0.0
            or row["max_runtime_integrated_qg_error_dex"] > 0.0
            for row in aval_dens
        ),
        "maxima_by_variant": {
            variant: {
                key: max(float(row[key]) for row in summaries if row[
                    "candidate_variant"
                ] == variant)
                for key in (
                    "max_potential_absolute_delta_V",
                    "max_density_error_dex",
                    "max_mobility_error_dex",
                    "max_alpha_error_dex",
                    "max_generation_error_dex",
                    "max_edge_sg_current_error_dex",
                    "max_runtime_integrated_qg_error_dex",
                    "max_terminal_current_relative_error",
                    "max_electric_field_vector_relative_error",
                    "max_electron_qf_gradient_vector_relative_error",
                    "max_hole_qf_gradient_vector_relative_error",
                    "max_electron_element_current_vector_relative_error",
                    "max_hole_element_current_vector_relative_error",
                )
            }
            for variant in CANDIDATE_VARIANTS
        },
        "input_sha256": input_hashes,
        "output_sha256": {
            long_path.name: sha256(long_path),
            summary_path.name: sha256(summary_path),
        },
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    return manifest


def main() -> int:
    args = parse_args()
    manifest = run(args.raw_root, args.output)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
