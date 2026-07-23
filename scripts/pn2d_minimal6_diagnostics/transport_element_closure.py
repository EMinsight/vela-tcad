"""Audit native Sentaurus element transport fields on identical support."""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path

from .edge_flux_experiment import (
    EXPECTED_STATES,
    _load_mesh,
    _load_observations,
    _quantile,
    _sha256,
    _value,
    _write_csv,
)


ELEMENTARY_CHARGE_C = 1.602176634e-19


def effective_density(
    current: tuple[float, float],
    mobility: float,
    grad_qf: tuple[float, float],
) -> dict[str, float | str | None]:
    grad2 = grad_qf[0] ** 2 + grad_qf[1] ** 2
    current2 = current[0] ** 2 + current[1] ** 2
    if mobility <= 0.0 or grad2 == 0.0 or current2 == 0.0:
        return {
            "status": "degenerate",
            "density_m3": None,
            "angle_deg": None,
            "orthogonal_residual": None,
        }
    dot = current[0] * grad_qf[0] + current[1] * grad_qf[1]
    cosine = max(-1.0, min(1.0, dot / math.sqrt(current2 * grad2)))
    density = dot / (ELEMENTARY_CHARGE_C * mobility * grad2)
    projected = (
        ELEMENTARY_CHARGE_C * mobility * density * grad_qf[0],
        ELEMENTARY_CHARGE_C * mobility * density * grad_qf[1],
    )
    residual = math.hypot(
        current[0] - projected[0], current[1] - projected[1]
    ) / math.sqrt(current2)
    return {
        "status": "valid" if density > 0.0 else "sign_incompatible",
        "density_m3": density if density > 0.0 else None,
        "angle_deg": math.degrees(math.acos(cosine)),
        "orthogonal_residual": residual,
    }


def _node_density_models(values: list[float]) -> dict[str, float]:
    if len(values) != 3 or any(value <= 0.0 for value in values):
        raise ValueError("cell density models require three positive values")
    return {
        "arithmetic": statistics.fmean(values),
        "geometric": math.exp(statistics.fmean(math.log(v) for v in values)),
        "harmonic": 3.0 / sum(1.0 / value for value in values),
    }


def _load_transport_rows(
    path: Path,
) -> dict[tuple[str, float, int], dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    result = {
        (row["topology"], float(row["bias_V"]), int(row["cell_id"])): row
        for row in rows
    }
    expected = {
        (topology, bias, cell)
        for topology, bias in EXPECTED_STATES
        for cell in range(4)
    }
    if len(rows) != len(result) or set(result) != expected:
        raise ValueError("transport rows differ from the exact 40 x 4 contract")
    return result


def _summary(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for topology in ("all", "mirror", "sketch"):
        for carrier in ("electron", "hole"):
            selected = [
                row
                for row in samples
                if row["carrier"] == carrier
                and (topology == "all" or row["topology"] == topology)
            ]
            valid = [row for row in selected if row["status"] == "valid"]
            angles = [float(row["current_grad_angle_deg"]) for row in valid]
            residuals = [
                float(row["orthogonal_current_residual"]) for row in valid
            ]
            record: dict[str, object] = {
                "topology": topology,
                "carrier": carrier,
                "sample_count": len(selected),
                "valid_count": len(valid),
                "sign_incompatible_count": sum(
                    row["status"] == "sign_incompatible" for row in selected
                ),
                "degenerate_count": sum(
                    row["status"] == "degenerate" for row in selected
                ),
                "median_current_grad_angle_deg": statistics.median(angles),
                "p95_current_grad_angle_deg": _quantile(angles, 0.95),
                "median_orthogonal_current_residual": statistics.median(
                    residuals
                ),
                "p95_orthogonal_current_residual": _quantile(residuals, 0.95),
            }
            for model in ("arithmetic", "geometric", "harmonic"):
                gaps = [
                    abs(float(row[f"effective_over_{model}_density_dex"]))
                    for row in valid
                ]
                record[f"median_abs_{model}_density_gap_dex"] = (
                    statistics.median(gaps)
                )
                record[f"p95_abs_{model}_density_gap_dex"] = _quantile(
                    gaps, 0.95
                )
            rows.append(record)
    return rows


def _report(summary: list[dict[str, object]]) -> str:
    selected = [row for row in summary if row["topology"] == "all"]
    lines = [
        "# PN2D Minimal6 native element transport closure",
        "",
        "Sentaurus exposes current, mobility, and quasi-Fermi-gradient fields "
        "on the same four-element support. It does not expose element carrier "
        "density for this model/deck, so density is independently inferred "
        "from J = q mu n g and compared with three node-to-cell controls.",
        "",
        "| Carrier | Valid | Angle median/p95 (deg) | Orthogonal residual "
        "median/p95 | Arithmetic density gap median/p95 (dex) | Geometric "
        "density gap median/p95 (dex) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in selected:
        lines.append(
            "| {carrier} | {valid_count}/{sample_count} | "
            "{median_current_grad_angle_deg:.6g} / "
            "{p95_current_grad_angle_deg:.6g} | "
            "{median_orthogonal_current_residual:.6g} / "
            "{p95_orthogonal_current_residual:.6g} | "
            "{median_abs_arithmetic_density_gap_dex:.6g} / "
            "{p95_abs_arithmetic_density_gap_dex:.6g} | "
            "{median_abs_geometric_density_gap_dex:.6g} / "
            "{p95_abs_geometric_density_gap_dex:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "The collinearity test is independent of carrier density. Density "
            "gap controls are diagnostic only because no native element "
            "density field is observable.",
            "",
            "This evidence does not provide a native directed-edge flux and "
            "does not authorize a production transport-formula change.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_transport_element_closure(
    *,
    transport_csv: str | Path,
    transport_manifest: str | Path,
    observations_csv: str | Path,
    inverse_inputs_root: str | Path,
    output_root: str | Path,
) -> dict[str, object]:
    transport_path = Path(transport_csv).resolve()
    transport_manifest_path = Path(transport_manifest).resolve()
    observations_path = Path(observations_csv).resolve()
    inverse_root = Path(inverse_inputs_root).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)

    manifest = json.loads(
        transport_manifest_path.read_text(encoding="utf-8")
    )
    if (
        manifest.get("status") != "valid"
        or manifest.get("state_count") != 40
        or manifest.get("sample_count") != 160
    ):
        raise ValueError("transport element manifest is not valid")
    transport = _load_transport_rows(transport_path)
    observations = _load_observations(observations_path)

    samples: list[dict[str, object]] = []
    for topology, bias in EXPECTED_STATES:
        _, triangles = _load_mesh(inverse_root, topology)
        for cell, triangle in enumerate(triangles):
            row = transport[(topology, bias, cell)]
            for carrier in ("electron", "hole"):
                prefix = "electron" if carrier == "electron" else "hole"
                density_quantity = "eDensity" if carrier == "electron" else "hDensity"
                current = (
                    float(row[f"{prefix}_current_x_A_per_m2"]),
                    float(row[f"{prefix}_current_y_A_per_m2"]),
                )
                grad = (
                    float(row[f"{prefix}_grad_qf_x_V_per_m"]),
                    float(row[f"{prefix}_grad_qf_y_V_per_m"]),
                )
                mobility = float(row[f"{prefix}_mobility_m2_per_Vs"])
                inversion = effective_density(current, mobility, grad)
                density_models = _node_density_models(
                    [
                        _value(
                            observations,
                            "sentaurus",
                            topology,
                            bias,
                            node,
                            density_quantity,
                        )
                        for node in triangle
                    ]
                )
                density = inversion["density_m3"]
                sample: dict[str, object] = {
                    "topology": topology,
                    "bias_V": bias,
                    "cell_id": cell,
                    "carrier": carrier,
                    "status": inversion["status"],
                    "current_x_A_per_m2": current[0],
                    "current_y_A_per_m2": current[1],
                    "grad_qf_x_V_per_m": grad[0],
                    "grad_qf_y_V_per_m": grad[1],
                    "mobility_m2_per_Vs": mobility,
                    "current_grad_angle_deg": (
                        "" if inversion["angle_deg"] is None
                        else inversion["angle_deg"]
                    ),
                    "orthogonal_current_residual": (
                        "" if inversion["orthogonal_residual"] is None
                        else inversion["orthogonal_residual"]
                    ),
                    "effective_density_m3": "" if density is None else density,
                }
                for model, value in density_models.items():
                    sample[f"{model}_node_density_m3"] = value
                    sample[f"effective_over_{model}_density_dex"] = (
                        "" if density is None else math.log10(float(density) / value)
                    )
                samples.append(sample)

    if len(samples) != 320:
        raise ValueError("transport closure count differs from 40 x 4 x 2")
    summary = _summary(samples)
    samples_path = output / "transport_element_closure_samples.csv"
    summary_path = output / "transport_element_closure_summary.csv"
    report_path = output / "report.md"
    _write_csv(samples_path, samples)
    _write_csv(summary_path, summary)
    report_path.write_text(_report(summary), encoding="utf-8")
    outputs = {
        samples_path.name: _sha256(samples_path),
        summary_path.name: _sha256(summary_path),
        report_path.name: _sha256(report_path),
    }
    result: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "minimal6_native_element_transport_closure",
        "state_count": 40,
        "sample_count": len(samples),
        "support_audit": {
            "current_mobility_grad_qf_same_native_element_support": True,
            "native_element_density_available": False,
            "native_directed_edge_flux_available": False,
        },
        "acceptance_policy": {
            "formula_change_authorized": False,
            "reason": (
                "The same-cell identity is observable, but native element "
                "density and directed-edge flux remain unavailable."
            ),
        },
        "inputs": {
            "transport_csv": str(transport_path),
            "transport_csv_sha256": _sha256(transport_path),
            "transport_manifest": str(transport_manifest_path),
            "transport_manifest_sha256": _sha256(transport_manifest_path),
            "observations_csv": str(observations_path),
            "observations_sha256": _sha256(observations_path),
            "inverse_inputs_root": str(inverse_root),
        },
        "outputs": outputs,
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result
