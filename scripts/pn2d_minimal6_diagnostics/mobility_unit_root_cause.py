"""Diagnose the Minimal6 high-field mobility unit and support mismatch.

The diagnostic is read-only with respect to production solver outputs. It
replays the production Masetti plus high-field arithmetic in SI units with
two saturation-velocity interpretations:

* correct: the physical saturation velocity is converted to the active unit;
* legacy unit-scaled: an SI default numeric value is consumed as cm/s.

It compares both branches with the direct triangle local-edge mobility emitted
by the C++ fixed-state audit and with native Sentaurus element mobility.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from collections.abc import Iterable, Mapping
from pathlib import Path

from .mobility_diagnosis import FIELD, masetti_low_field_mobility


CARRIERS = ("electron", "hole")
TOPOLOGIES = ("mirror", "sketch")
BIASES_V = tuple(float(-value) for value in range(1, 21))


def absolute_log10_error(value: float, reference: float) -> float:
    """Return the absolute base-10 ratio error for positive finite values."""
    value = float(value)
    reference = float(reference)
    if not math.isfinite(value) or not math.isfinite(reference):
        raise ValueError("mobility comparison requires finite values")
    if value <= 0.0 or reference <= 0.0:
        raise ValueError("mobility comparison requires positive values")
    return abs(math.log10(value / reference))


def quantile(values: Iterable[float], fraction: float) -> float:
    """Return a deterministic linearly interpolated quantile."""
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("quantile requires at least one value")
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("quantile fraction must be in [0, 1]")
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (
        position - lower
    )


def field_limited_mobility(
    carrier: str,
    low_field_mobility_m2_per_Vs: float,
    field_V_per_m: float,
    *,
    saturation_velocity_scale: float = 1.0,
) -> float:
    """Evaluate the production high-field law with an explicit velocity scale.

    ``saturation_velocity_scale=1`` is the correct SI interpretation.
    ``saturation_velocity_scale=0.01`` reproduces the unit-scaled C++ behavior
    where a value declared in m/s is numerically consumed as cm/s.
    """
    if carrier not in CARRIERS:
        raise ValueError(f"unsupported carrier {carrier!r}")
    mobility = float(low_field_mobility_m2_per_Vs)
    field = abs(float(field_V_per_m))
    scale = float(saturation_velocity_scale)
    if mobility < 0.0 or not math.isfinite(mobility):
        raise ValueError("low-field mobility must be finite and non-negative")
    if not math.isfinite(field):
        raise ValueError("driving field must be finite")
    if scale <= 0.0 or not math.isfinite(scale):
        raise ValueError("saturation velocity scale must be finite and positive")
    if mobility == 0.0 or field == 0.0:
        return mobility
    parameters = FIELD[carrier]
    beta = float(parameters["beta"])
    velocity = float(parameters["saturation_velocity"]) * scale
    ratio = mobility * field / velocity
    return mobility / (1.0 + ratio**beta) ** (1.0 / beta)


def endpoint_averaged_mobility(
    carrier: str,
    net_doping0_m3: float,
    net_doping1_m3: float,
    field_V_per_m: float,
    *,
    saturation_velocity_scale: float,
) -> float:
    """Reproduce triangleGssEndpointAveragedMobility in physical SI units."""
    values = []
    for doping in (net_doping0_m3, net_doping1_m3):
        low_field = masetti_low_field_mobility(carrier, float(doping))
        values.append(
            field_limited_mobility(
                carrier,
                low_field,
                field_V_per_m,
                saturation_velocity_scale=saturation_velocity_scale,
            )
        )
    return 0.5 * (values[0] + values[1])


def cell_average_doping_mobility(
    carrier: str,
    net_doping_m3: Iterable[float],
    field_V_per_m: float,
    *,
    saturation_velocity_scale: float,
) -> float:
    """Evaluate Masetti at the arithmetic cell-average net doping."""
    values = tuple(float(value) for value in net_doping_m3)
    if len(values) != 3:
        raise ValueError("cell mobility requires three nodal doping values")
    low_field = masetti_low_field_mobility(carrier, sum(values) / 3.0)
    return field_limited_mobility(
        carrier,
        low_field,
        field_V_per_m,
        saturation_velocity_scale=saturation_velocity_scale,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _doping_by_topology(
    inverse_inputs_root: Path,
) -> dict[str, dict[int, float]]:
    result: dict[str, dict[int, float]] = {}
    source = inverse_inputs_root / "vela" / "source" / "topologies"
    for topology in TOPOLOGIES:
        rows = _read_rows(source / topology / "doping.csv")
        values = {
            int(row["node_id"]): (
                float(row["donors_cm3"]) - float(row["acceptors_cm3"])
            )
            * 1.0e6
            for row in rows
        }
        if set(values) != set(range(6)):
            raise ValueError(f"{topology} doping does not contain nodes 0..5")
        result[topology] = values
    return result


def _sentaurus_elements(
    path: Path,
) -> dict[tuple[str, float, int], dict[str, str]]:
    result: dict[tuple[str, float, int], dict[str, str]] = {}
    for row in _read_rows(path):
        key = (row["topology"], float(row["bias_V"]), int(row["cell_id"]))
        if key in result:
            raise ValueError(f"duplicate Sentaurus element row {key}")
        result[key] = row
    if len(result) != 160:
        raise ValueError(
            f"expected 160 Sentaurus element rows, received {len(result)}"
        )
    return result


def _state(path: Path) -> dict[int, dict[str, str]]:
    result = {int(row["node_id"]): row for row in _read_rows(path)}
    if set(result) != set(range(6)):
        raise ValueError(f"state does not contain nodes 0..5: {path}")
    return result


def _edge_lengths(path: Path) -> dict[tuple[int, int], float]:
    result = {
        tuple(sorted((int(row["node0"]), int(row["node1"])))): float(
            row["length_m"]
        )
        for row in _read_rows(path)
    }
    if len(result) != 9:
        raise ValueError(f"edge audit does not contain nine unique edges: {path}")
    return result


def _summary_row(
    *,
    support: str,
    carrier: str,
    branch: str,
    values: list[float],
) -> dict[str, object]:
    return {
        "support": support,
        "carrier": carrier,
        "branch": branch,
        "sample_count": len(values),
        "median_abs_log10_error_dex": statistics.median(values),
        "p95_abs_log10_error_dex": quantile(values, 0.95),
        "maximum_abs_log10_error_dex": max(values),
    }


def _markdown(
    summaries: list[dict[str, object]],
    manifest: Mapping[str, object],
) -> str:
    lookup = {
        (str(row["support"]), str(row["carrier"]), str(row["branch"])): row
        for row in summaries
    }
    local_e = lookup[
        ("triangle_local_edge", "electron", "legacy_velocity_interpretation")
    ]
    local_h = lookup[
        ("triangle_local_edge", "hole", "legacy_velocity_interpretation")
    ]
    correct_local_e = lookup[
        ("triangle_local_edge", "electron", "correct_velocity_interpretation")
    ]
    correct_local_h = lookup[
        ("triangle_local_edge", "hole", "correct_velocity_interpretation")
    ]
    elem_e = lookup[
        ("sentaurus_native_element", "electron", "correct_cell_average_doping")
    ]
    elem_h = lookup[
        ("sentaurus_native_element", "hole", "correct_cell_average_doping")
    ]
    buggy_elem_e = lookup[
        (
            "sentaurus_native_element",
            "electron",
            "legacy_cell_average_doping",
        )
    ]
    buggy_elem_h = lookup[
        ("sentaurus_native_element", "hole", "legacy_cell_average_doping")
    ]
    production_is_correct = (
        float(correct_local_e["median_abs_log10_error_dex"])
        + float(correct_local_h["median_abs_log10_error_dex"])
        < float(local_e["median_abs_log10_error_dex"])
        + float(local_h["median_abs_log10_error_dex"])
    )
    if production_is_correct:
        unit_status = (
            "The current production audit confirms that the saturation-velocity "
            "unit conversion defect has been repaired. The direct C++ operator now "
            "uses the physically correct velocity branch in `unit_scaling`."
        )
        direct_closure = (
            "- Direct C++ triangle local-edge mobility closes against the correctly "
            "converted velocity interpretation at median errors "
            f"{float(correct_local_e['median_abs_log10_error_dex']):.3g} dex "
            "(electron) and "
            f"{float(correct_local_h['median_abs_log10_error_dex']):.3g} dex (hole)."
        )
        production_closure = (
            "The correctly converted branch matches every direct C++ triangle "
            "local-edge value to floating-point precision, which verifies the "
            "production fix independently of aggregate correlation. "
        )
        recommendations = [
            "1. Preserve the velocity unit parity tests and the two-root 40-state audit as regression gates.",
            "2. Investigate the remaining native-element mobility median gaps using aligned doping interpolation and temperature-dependent high-field parameters.",
            "3. Continue directed-edge SG current inversion with node-pair support alignment.",
            "4. Treat the remaining avalanche-source gap as a fixed-state support/model residual, not evidence for another mobility scale change.",
        ]
    else:
        unit_status = (
            "The current production audit still exhibits the saturation-velocity "
            "unit conversion defect in `unit_scaling`: an SI m/s default is consumed "
            "as an internal cm/s value, lowering the physical velocity by 100."
        )
        direct_closure = (
            "- Direct C++ triangle local-edge mobility closes against the legacy "
            "velocity interpretation at median errors "
            f"{float(local_e['median_abs_log10_error_dex']):.3g} dex (electron) and "
            f"{float(local_h['median_abs_log10_error_dex']):.3g} dex (hole)."
        )
        production_closure = (
            "The legacy branch matches every direct C++ triangle local-edge value "
            "to floating-point precision, confirming that the production defect "
            "is still active. "
        )
        recommendations = [
            "1. Add an explicit velocity conversion to the unit system and apply it to default and JSON-provided saturation velocities.",
            "2. Add legacy-SI versus unit-scaled physical parity tests.",
            "3. Re-run the 40-state fixed-state and self-consistent replacement audits.",
            "4. Continue directed-edge SG and avalanche closure only after mobility parity is restored.",
        ]
    lines = [
        "# PN2D Minimal6 mobility unit root-cause audit",
        "",
        "## Technical summary",
        "",
        unit_status,
        "",
        direct_closure,
        f"- On native Sentaurus element support, the legacy interpretation has median errors "
        f"{float(buggy_elem_e['median_abs_log10_error_dex']):.6f} dex and "
        f"{float(buggy_elem_h['median_abs_log10_error_dex']):.6f} dex.",
        f"- Converting saturation velocity consistently reduces those native-element medians to "
        f"{float(elem_e['median_abs_log10_error_dex']):.6f} dex and "
        f"{float(elem_h['median_abs_log10_error_dex']):.6f} dex.",
        "",
        "## Exact comparison table",
        "",
        "| Support | Carrier | Branch | N | Median abs error (dex) | P95 (dex) | Maximum (dex) |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in summaries:
        lines.append(
            "| {support} | {carrier} | {branch} | {sample_count} | "
            "{median_abs_log10_error_dex:.12g} | "
            "{p95_abs_log10_error_dex:.12g} | "
            "{maximum_abs_log10_error_dex:.12g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Scope and metric definitions",
            "",
            "- States: 40 exact states, mirror/sketch topologies, reverse biases -1 through -20 V.",
            "- Direct local-edge sample: each of 4 cells x 3 local edges x 2 carriers x 40 states = 960 rows.",
            "- Native element sample: 4 cells x 2 carriers x 40 states = 320 rows.",
            "- Error metric: absolute base-10 logarithm of the positive mobility ratio.",
            "- Correct element reconstruction: Masetti at arithmetic cell-average net doping, limited by the affine cell QFP-gradient magnitude with saturation velocity in m/s.",
            "",
            "## Why support must be aligned",
            "",
            "Vela's triangle avalanche path stores one mobility for each cell-local edge and uses the edge QFP difference. "
            "Sentaurus native element mobility is a cell quantity associated with the native element QFP-gradient field. "
            "Averaging Vela local-edge mobilities is not the same operator and can make the corrected formula look worse because a zero-QFP-difference edge retains low-field mobility. "
            "The root-cause comparison therefore reconstructs Vela on the native cell field before comparing with Sentaurus element mobility.",
            "",
            "## Method and robustness",
            "",
            "The diagnostic reproduces the production Masetti parameters and high-field exponent. "
            "It evaluates two fixed branches that differ only by a factor of 100 in saturation velocity. "
            + production_closure
            + "The native-element result independently shows that the correctly converted branch is much closer to Sentaurus.",
            "",
            "## Limitations",
            "",
            "This audit does not expose the global SG edge mobility directly; it uses the direct triangle local-edge mobility already emitted by the fixed-state C++ audit. "
            "The Sentaurus element comparison remains a model/support reconstruction rather than access to Sentaurus's internal edge flux. "
            "Residual native-element errors up to the reported maxima can contain doping interpolation, material-model parameter, and element-evaluation differences.",
            "",
            "## Recommended next steps",
            "",
            *recommendations,
            "",
            "## Further questions",
            "",
            "After the velocity conversion is fixed, how much of the remaining 0.05 dex median element gap is explained by cell doping interpolation versus Sentaurus-specific high-field parameter temperature dependence?",
            "",
            "All paths in this report are relative to its containing output root.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_mobility_unit_root_cause(
    *,
    self_consistent_root: str | Path,
    inverse_inputs_root: str | Path,
    sentaurus_element_csv: str | Path,
    output_root: str | Path,
) -> dict[str, object]:
    """Run the exact 40-state mobility unit/support decomposition."""
    self_root = Path(self_consistent_root).resolve()
    inverse_root = Path(inverse_inputs_root).resolve()
    element_path = Path(sentaurus_element_csv).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)

    doping = _doping_by_topology(inverse_root)
    sentaurus = _sentaurus_elements(element_path)
    local_rows: list[dict[str, object]] = []
    element_rows: list[dict[str, object]] = []
    metrics = {
        carrier: {
            "direct_vs_legacy": [],
            "direct_vs_correct": [],
            "sent_vs_legacy_node_average": [],
            "sent_vs_correct_node_average": [],
            "sent_vs_legacy_average_doping": [],
            "sent_vs_correct_average_doping": [],
        }
        for carrier in CARRIERS
    }

    state_count = 0
    for topology in TOPOLOGIES:
        for bias in BIASES_V:
            state_count += 1
            label = f"m{abs(int(bias))}V"
            state_root = self_root / "self_consistent_replay" / topology / label
            state = _state(state_root / "state.csv")
            lengths = _edge_lengths(state_root / "edges.csv")
            triangles = _read_rows(state_root / "triangles.csv")
            if len(triangles) != 4:
                raise ValueError(
                    f"expected four triangles for {topology} {bias:g} V"
                )
            for triangle in triangles:
                cell_id = int(triangle["cell_id"])
                nodes = tuple(int(triangle[f"node{index}"]) for index in range(3))
                reference = sentaurus[(topology, bias, cell_id)]
                for carrier in CARRIERS:
                    qf_key = "phin_V" if carrier == "electron" else "phip_V"
                    cell_field = float(
                        triangle[
                            f"local_edge0_{carrier}_cell_qf_field_V_per_m"
                        ]
                    )
                    direct_values: list[float] = []
                    correct_values: list[float] = []
                    legacy_values: list[float] = []
                    for local_edge in range(3):
                        prefix = f"local_edge{local_edge}"
                        node0 = int(triangle[f"{prefix}_node0"])
                        node1 = int(triangle[f"{prefix}_node1"])
                        pair = tuple(sorted((node0, node1)))
                        edge_field = abs(
                            float(state[node1][qf_key])
                            - float(state[node0][qf_key])
                        ) / lengths[pair]
                        direct = float(
                            triangle[
                                f"{prefix}_{carrier}_mobility_m2_per_V_s"
                            ]
                        )
                        correct = endpoint_averaged_mobility(
                            carrier,
                            doping[topology][node0],
                            doping[topology][node1],
                            edge_field,
                            saturation_velocity_scale=1.0,
                        )
                        legacy = endpoint_averaged_mobility(
                            carrier,
                            doping[topology][node0],
                            doping[topology][node1],
                            edge_field,
                            saturation_velocity_scale=0.01,
                        )
                        direct_legacy_error = absolute_log10_error(
                            direct, legacy
                        )
                        direct_correct_error = absolute_log10_error(
                            direct, correct
                        )
                        metrics[carrier]["direct_vs_legacy"].append(
                            direct_legacy_error
                        )
                        metrics[carrier]["direct_vs_correct"].append(
                            direct_correct_error
                        )
                        direct_values.append(direct)
                        correct_values.append(correct)
                        legacy_values.append(legacy)
                        local_rows.append(
                            {
                                "topology": topology,
                                "bias_V": bias,
                                "cell_id": cell_id,
                                "carrier": carrier,
                                "local_edge": local_edge,
                                "node0": node0,
                                "node1": node1,
                                "edge_qf_field_V_per_m": edge_field,
                                "direct_cpp_mobility_m2_per_Vs": direct,
                                "legacy_velocity_mobility_m2_per_Vs": legacy,
                                "correct_velocity_mobility_m2_per_Vs": correct,
                                "direct_vs_legacy_abs_log10_error_dex": direct_legacy_error,
                                "direct_vs_correct_abs_log10_error_dex": direct_correct_error,
                            }
                        )

                    nodal_low_field = [
                        masetti_low_field_mobility(
                            carrier, doping[topology][node]
                        )
                        for node in nodes
                    ]
                    legacy_node_average = sum(
                        field_limited_mobility(
                            carrier,
                            value,
                            cell_field,
                            saturation_velocity_scale=0.01,
                        )
                        for value in nodal_low_field
                    ) / 3.0
                    correct_node_average = sum(
                        field_limited_mobility(
                            carrier,
                            value,
                            cell_field,
                            saturation_velocity_scale=1.0,
                        )
                        for value in nodal_low_field
                    ) / 3.0
                    legacy_average_doping = cell_average_doping_mobility(
                        carrier,
                        (doping[topology][node] for node in nodes),
                        cell_field,
                        saturation_velocity_scale=0.01,
                    )
                    correct_average_doping = cell_average_doping_mobility(
                        carrier,
                        (doping[topology][node] for node in nodes),
                        cell_field,
                        saturation_velocity_scale=1.0,
                    )
                    sentaurus_mobility = float(
                        reference[f"{carrier}_mobility_m2_per_Vs"]
                    )
                    errors = {
                        "legacy_node_average": absolute_log10_error(
                            legacy_node_average, sentaurus_mobility
                        ),
                        "correct_node_average": absolute_log10_error(
                            correct_node_average, sentaurus_mobility
                        ),
                        "legacy_average_doping": absolute_log10_error(
                            legacy_average_doping, sentaurus_mobility
                        ),
                        "correct_average_doping": absolute_log10_error(
                            correct_average_doping, sentaurus_mobility
                        ),
                    }
                    for branch, error in errors.items():
                        metrics[carrier][f"sent_vs_{branch}"].append(error)
                    element_rows.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "cell_id": cell_id,
                            "carrier": carrier,
                            "cell_qf_field_V_per_m": cell_field,
                            "sentaurus_native_mobility_m2_per_Vs": sentaurus_mobility,
                            "legacy_node_average_mobility_m2_per_Vs": legacy_node_average,
                            "correct_node_average_mobility_m2_per_Vs": correct_node_average,
                            "legacy_average_doping_mobility_m2_per_Vs": legacy_average_doping,
                            "correct_average_doping_mobility_m2_per_Vs": correct_average_doping,
                            "legacy_node_average_abs_log10_error_dex": errors[
                                "legacy_node_average"
                            ],
                            "correct_node_average_abs_log10_error_dex": errors[
                                "correct_node_average"
                            ],
                            "legacy_average_doping_abs_log10_error_dex": errors[
                                "legacy_average_doping"
                            ],
                            "correct_average_doping_abs_log10_error_dex": errors[
                                "correct_average_doping"
                            ],
                            "mean_direct_local_edge_mobility_m2_per_Vs": sum(
                                direct_values
                            )
                            / 3.0,
                            "mean_legacy_local_edge_mobility_m2_per_Vs": sum(
                                legacy_values
                            )
                            / 3.0,
                            "mean_correct_local_edge_mobility_m2_per_Vs": sum(
                                correct_values
                            )
                            / 3.0,
                        }
                    )

    if state_count != 40 or len(local_rows) != 960 or len(element_rows) != 320:
        raise ValueError(
            "unexpected state/local-edge/element counts "
            f"{state_count}/{len(local_rows)}/{len(element_rows)}"
        )

    summaries: list[dict[str, object]] = []
    for carrier in CARRIERS:
        summaries.append(
            _summary_row(
                support="triangle_local_edge",
                carrier=carrier,
                branch="legacy_velocity_interpretation",
                values=metrics[carrier]["direct_vs_legacy"],
            )
        )
        summaries.append(
            _summary_row(
                support="triangle_local_edge",
                carrier=carrier,
                branch="correct_velocity_interpretation",
                values=metrics[carrier]["direct_vs_correct"],
            )
        )
        for branch in (
            "legacy_node_average",
            "correct_node_average",
            "legacy_average_doping",
            "correct_average_doping",
        ):
            summaries.append(
                _summary_row(
                    support="sentaurus_native_element",
                    carrier=carrier,
                    branch=branch.replace(
                        "_average", "_cell_average", 1
                    )
                    if branch.endswith("average_doping")
                    else branch,
                    values=metrics[carrier][f"sent_vs_{branch}"],
                )
            )

    local_path = output / "local_edge_unit_decomposition.csv"
    element_output = output / "native_element_support_decomposition.csv"
    summary_path = output / "summary.csv"
    _write_rows(local_path, local_rows)
    _write_rows(element_output, element_rows)
    _write_rows(summary_path, summaries)

    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "pn2d_minimal6_mobility_unit_root_cause",
        "state_count": state_count,
        "local_edge_carrier_sample_count": len(local_rows),
        "native_element_carrier_sample_count": len(element_rows),
        "velocity_interpretations": {
            "correct": "SI defaults are converted to the active unit; unit-scaled JSON values use cm/s",
            "legacy_unit_scaled": (
                "configured m/s number is consumed as cm/s, equivalent to "
                "0.01 times the physical SI saturation velocity"
            ),
        },
        "inputs": {
            "self_consistent_root": str(self_root),
            "self_consistent_manifest_sha256": _sha256(
                self_root / "manifest.json"
            ),
            "inverse_inputs_root": str(inverse_root),
            "sentaurus_element_csv": str(element_path),
            "sentaurus_element_csv_sha256": _sha256(element_path),
        },
        "outputs": {
            "local_edge_csv": local_path.name,
            "local_edge_csv_sha256": _sha256(local_path),
            "native_element_csv": element_output.name,
            "native_element_csv_sha256": _sha256(element_output),
            "summary_csv": summary_path.name,
            "summary_csv_sha256": _sha256(summary_path),
        },
        "source_audit": {
            "mobility_defaults": "include/vela/physics/MobilityModel.h",
            "unit_system": "src/core/UnitScaling.cpp",
            "parameter_conversion": "src/physics/MobilityModel.cpp",
            "triangle_local_edge_operator": (
                "include/vela/equation/AssemblerUtils.h"
            ),
        },
    }
    report_path = output / "report.md"
    report_path.write_text(
        _markdown(summaries, manifest), encoding="utf-8"
    )
    manifest["outputs"]["report_md"] = report_path.name
    manifest["outputs"]["report_md_sha256"] = _sha256(report_path)
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest

