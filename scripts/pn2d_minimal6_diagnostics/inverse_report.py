"""Deterministic writers and report assembly for the Minimal6 inverse audit."""

from __future__ import annotations

import csv
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import statistics
from typing import Iterable

from .inverse_avalanche import impact_generation, invert_van_overstraeten_alpha
from .inverse_contracts import (
    AcceptanceThresholds, Identifiability, Observation, SampleStatus,
    SupportKind, validate_inverse_report_v1,
)
from .inverse_fields import triangle_gradient
from .inverse_inputs import (
    InputBundle, canonical_observations, field_inventory, write_input_manifest,
)
from .inverse_plots import FIGURE_NAMES, render_inverse_figures, write_figure_manifest
from .inverse_replacements import INVERSE_DEPENDENCIES, run_replacement_matrix
from .inverse_transport import current_inverted_qf_gradient


OBSERVATION_COLUMNS = (
    "solver", "topology", "bias_V", "support_kind", "support_id", "quantity",
    "component", "raw_value", "raw_unit", "value_si", "unit_si",
    "coordinate_frame", "orientation", "conversion", "status", "source_path",
    "source_sha256",
)
CANDIDATE_COLUMNS = (
    "candidate", "quantity", "carrier", "split", "topology", "bias_V",
    "support_kind", "valid_count", "median_abs_error", "p95_abs_error",
    "median_angle_deg", "classification",
)
REPLACEMENT_COLUMNS = (
    "sequence", "step", "factor", "value", "incremental_dex", "closure_abs_dex",
)


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path: str | Path, payload: object) -> None:
    Path(path).write_text(
        json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8", newline="\n",
    )


def _format_csv(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("CSV values must be finite")
        return format(value, ".17g")
    if hasattr(value, "value"):
        return str(value.value)
    return str(value)


def write_csv(path: str | Path, columns: tuple[str, ...], rows: Iterable[dict]) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="raise",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _format_csv(row.get(column)) for column in columns})


def _observation_row(row: Observation) -> dict:
    return {
        "solver": row.solver, "topology": row.topology, "bias_V": row.bias_V,
        "support_kind": row.support_kind.value, "support_id": row.support_id,
        "quantity": row.quantity, "component": row.component,
        "raw_value": row.raw_value, "raw_unit": row.raw_unit, "value_si": row.value_si,
        "unit_si": row.unit_si, "coordinate_frame": row.coordinate_frame,
        "orientation": row.orientation, "conversion": row.conversion,
        "status": row.status.value, "source_path": row.source_path,
        "source_sha256": row.source_sha256,
    }


def _index(rows: Iterable[Observation]) -> dict[tuple, Observation]:
    return {(row.solver, row.topology, row.bias_V, str(row.support_id), row.quantity,
             row.component): row for row in rows}


def _finite(row: Observation | None) -> float | None:
    if row is None or row.status is not SampleStatus.VALID or row.value_si is None:
        return None
    value = float(row.value_si)
    return value if math.isfinite(value) else None


def _percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _paired_errors(rows: tuple[Observation, ...], quantities: set[str], split_keys: set[tuple[str, float]]) -> list[float]:
    index = _index(rows)
    errors = []
    for key, vela in sorted(index.items(), key=lambda item: tuple(map(str, item[0]))):
        solver, topology, bias, node, quantity, component = key
        if solver != "vela" or quantity not in quantities or (topology, bias) not in split_keys:
            continue
        sentaurus = index.get(("sentaurus", topology, bias, node, quantity, component))
        first, second = _finite(vela), _finite(sentaurus)
        if first is None or second is None or abs(second) <= 1.0e-300:
            continue
        errors.append(abs(math.log10(max(abs(first), 1.0e-300) / max(abs(second), 1.0e-300))))
    return errors


def _candidate_metrics(bundle: InputBundle, rows: tuple[Observation, ...], limits: AcceptanceThresholds) -> tuple[list[dict], list[dict]]:
    specifications = (
        ("potential_field_direct", "potential_and_field", "", {"ElectrostaticPotential", "ElectricField"}, limits.gradient_median_abs_dex),
        ("current_density_direct", "current_density", "both", {"eCurrentDensity", "hCurrentDensity"}, limits.gradient_median_abs_dex),
        ("alpha_generation_direct", "alpha_and_generation", "both", {"eAlphaAvalanche", "hAlphaAvalanche", "ImpactIonization"}, limits.local_generation_abs_dex),
    )
    result = []
    classifications = []
    split_sets = {
        "discovery": set(bundle.discovery_keys), "holdout": set(bundle.holdout_keys),
        "combined": set(bundle.common_keys),
    }
    candidate_classification = {}
    for candidate, quantity, carrier, fields, gate in specifications:
        split_rows = {}
        for split in ("discovery", "holdout", "combined"):
            errors_by_quantity = [
                _paired_errors(rows, {field}, split_sets[split])
                for field in sorted(fields)
            ]
            errors = [error for field_errors in errors_by_quantity for error in field_errors]
            complete_coverage = all(errors_by_quantity)
            median = statistics.median(errors) if errors else None
            p95 = _percentile(errors, 0.95)
            classification = (Identifiability.INSUFFICIENT_DATA if not complete_coverage
                              else Identifiability.IDENTIFIED if median is not None and p95 is not None
                              and median <= gate and p95 <= limits.gradient_p95_abs_dex
                              else Identifiability.REJECTED)
            split_rows[split] = classification
            result.append({
                "candidate": candidate, "quantity": quantity, "carrier": carrier,
                "split": split, "topology": "all", "bias_V": None,
                "support_kind": SupportKind.NODE.value, "valid_count": len(errors),
                "median_abs_error": median, "p95_abs_error": p95,
                "median_angle_deg": None, "classification": classification.value,
            })
        if any(value is Identifiability.INSUFFICIENT_DATA for value in split_rows.values()):
            final = Identifiability.INSUFFICIENT_DATA
        elif all(value is Identifiability.IDENTIFIED for value in split_rows.values()):
            final = Identifiability.IDENTIFIED
        else:
            final = Identifiability.REJECTED
        candidate_classification[candidate] = final

    # Mobility is independently available only for the Sentaurus rows.  The
    # Vela-to-Sentaurus gradient inference therefore remains a combined product.
    qf_class = Identifiability.CONFOUNDED
    result.append({
        "candidate": "current_inverted_qf_gradient", "quantity": "qf_gradient",
        "carrier": "both", "split": "combined", "topology": "all", "bias_V": None,
        "support_kind": SupportKind.NODE.value, "valid_count": 0,
        "median_abs_error": None, "p95_abs_error": None, "median_angle_deg": None,
        "classification": qf_class.value,
    })
    candidate_classification["current_inverted_qf_gradient"] = qf_class
    for candidate in sorted(candidate_classification):
        classification = candidate_classification[candidate]
        classifications.append({
            "candidate": candidate, "classification": classification.value,
            "claim_type": "identifiability",
            "reason": (
                "discovery, holdout, and combined numerical gates passed without local fitting"
                if classification is Identifiability.IDENTIFIED else
                "mobility and gradient are not independently available for both solvers"
                if classification is Identifiability.CONFOUNDED else
                "one or more declared quantities lacked compatible finite paired support in a required split"
                if classification is Identifiability.INSUFFICIENT_DATA else
                "one or more discovery, holdout, or combined gates failed"
            ),
        })
    return result, classifications


def _replacement() -> dict:
    def operand(factor: str, value: float) -> dict:
        return {
            "factor": factor, "value": value, "status": SampleStatus.VALID.value,
            "support_kind": SupportKind.INTEGRATED.value, "support_id": "global",
            "unit_si": "dimensionless", "carrier": None, "topology": "all", "bias_V": -20.0,
        }
    baseline = {factor: operand(factor, 1.0) for factor in INVERSE_DEPENDENCIES}
    replacement = {factor: operand(factor, 1.0 + (index + 1) / 100.0)
                   for index, factor in enumerate(INVERSE_DEPENDENCIES)}
    target = math.prod(row["value"] for row in replacement.values())
    return run_replacement_matrix(baseline, replacement, direct_target=target)


def _semantic_replay(rows: tuple[Observation, ...], replacement: dict) -> dict:
    index = _index(rows)
    sentaurus_states = sorted({(row.topology, row.bias_V) for row in rows if row.solver == "sentaurus"})
    state = sentaurus_states[0]
    nodes = sorted({str(row.support_id) for row in rows if row.solver == "sentaurus"
                    and (row.topology, row.bias_V) == state and row.quantity == "coordinate"})
    replay: dict[str, object] = {"state": [state[0], state[1]]}
    if len(nodes) >= 3:
        points, potentials = [], []
        for node in nodes[:3]:
            x = _finite(index.get(("sentaurus", *state, node, "coordinate", "x")))
            y = _finite(index.get(("sentaurus", *state, node, "coordinate", "y")))
            psi = _finite(index.get(("sentaurus", *state, node, "ElectrostaticPotential", "component0")))
            points.append([x, y]); potentials.append(psi)
        try:
            gradient = triangle_gradient(points, potentials)
            replay["triangle_gradient"] = {"status": "valid", "points_m": points,
                                           "values_V": potentials, "value_V_per_m": list(gradient)}
        except (TypeError, ValueError):
            replay["triangle_gradient"] = {"status": "incompatible_support"}
    else:
        replay["triangle_gradient"] = {"status": "insufficient_data"}

    node = nodes[0]
    density = _finite(index.get(("sentaurus", *state, node, "eDensity", "component0")))
    mobility = _finite(index.get(("sentaurus", *state, node, "eMobility", "component0")))
    jx = _finite(index.get(("sentaurus", *state, node, "eCurrentDensity", "component0")))
    jy = _finite(index.get(("sentaurus", *state, node, "eCurrentDensity", "component1")))
    if None not in (density, mobility, jx, jy) and density > 0.0 and mobility > 0.0:
        gradient = current_inverted_qf_gradient("electron", density, mobility, (jx, jy))
        replay["current_inverted_gradient"] = {
            "status": "valid", "carrier": "electron", "density_m3": density,
            "mobility_m2_per_Vs": mobility, "current_A_per_m2": [jx, jy],
            "value_V_per_m": list(gradient),
        }
    else:
        replay["current_inverted_gradient"] = {"status": "insufficient_data"}

    alpha_n = _finite(index.get(("sentaurus", *state, node, "eAlphaAvalanche", "component0")))
    alpha_p = _finite(index.get(("sentaurus", *state, node, "hAlphaAvalanche", "component0")))
    hjx = _finite(index.get(("sentaurus", *state, node, "hCurrentDensity", "component0")))
    hjy = _finite(index.get(("sentaurus", *state, node, "hCurrentDensity", "component1")))
    if alpha_n is not None and alpha_n > 0.0:
        prefactor = max(2.0 * alpha_n, 1.0)
        field, status = invert_van_overstraeten_alpha(alpha_n, prefactor=prefactor,
                                                       critical_field=1.0, gamma=1.0)
        replay["inverse_alpha"] = {
            "status": status, "alpha_m_inv": alpha_n, "prefactor_m_inv": prefactor,
            "critical_field_V_per_m": 1.0, "gamma": 1.0, "field_V_per_m": field,
        }
    else:
        replay["inverse_alpha"] = {"status": "insufficient_data"}
    if None not in (alpha_n, alpha_p, jx, jy, hjx, hjy):
        generation = impact_generation(alpha_n, (jx, jy), alpha_p, (hjx, hjy))
        replay["generation"] = {
            "status": "valid", "alpha_n_m_inv": alpha_n, "alpha_p_m_inv": alpha_p,
            "jn_A_per_m2": [jx, jy], "jp_A_per_m2": [hjx, hjy],
            "value_m3_s_inv": generation,
        }
    else:
        replay["generation"] = {"status": "insufficient_data"}
    generation_rows = [
        _finite(index.get(("sentaurus", *state, current, "ImpactIonization", "component0")))
        for current in nodes[:3]
    ]
    triangle = replay["triangle_gradient"]
    if (triangle["status"] == "valid" and len(generation_rows) == 3
            and all(value is not None for value in generation_rows)):
        (x0, y0), (x1, y1), (x2, y2) = triangle["points_m"]
        area_m2 = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        depth_m = 0.01
        replay["generation_support_integral"] = {
            "status": "valid", "support_kind": "triangle_from_raw_nodes",
            "support_ids": nodes[:3], "native_generation_m3_s_inv": generation_rows,
            "area_m2": area_m2, "depth_m": depth_m,
            "integral_s_inv": statistics.mean(generation_rows) * area_m2 * depth_m,
        }
    else:
        replay["generation_support_integral"] = {"status": "insufficient_data"}
    replay["replacement_closure"] = replacement["closure"]
    return replay


def _series(rows: tuple[Observation, ...], quantity: str, *, magnitude: bool = False) -> tuple[list[float], list[dict]]:
    biases = sorted({abs(row.bias_V) for row in rows if row.quantity == quantity})
    output = []
    for solver in ("vela", "sentaurus"):
        values = []
        for bias in biases:
            selected = [row for row in rows if row.solver == solver and abs(row.bias_V) == bias
                        and row.quantity == quantity and row.status is SampleStatus.VALID
                        and row.value_si is not None]
            if magnitude:
                grouped = {}
                for row in selected:
                    grouped.setdefault((row.topology, row.support_id), {})[row.component] = row.value_si
                complete = [value for value in grouped.values()
                            if value.get("component0") is not None
                            and value.get("component1") is not None]
                values.append(statistics.mean(math.hypot(float(value["component0"]),
                                                         float(value["component1"]))
                                              for value in complete) if complete else None)
            else:
                values.append(statistics.mean(abs(float(row.value_si)) for row in selected) if selected else None)
        output.append({"label": solver.capitalize(), "values": values})
    return biases, output


def _qf_series(rows: tuple[Observation, ...]) -> tuple[list[float], list[dict]]:
    index = _index(rows)
    biases = sorted({abs(row.bias_V) for row in rows if row.solver == "sentaurus"})
    result = []
    for carrier, prefix in (("electron", "e"), ("hole", "h")):
        values = []
        for bias in biases:
            samples = []
            keys = sorted({(row.topology, str(row.support_id)) for row in rows
                           if row.solver == "sentaurus" and abs(row.bias_V) == bias})
            for topology, node in keys:
                state = (topology, -bias)
                density = _finite(index.get(("sentaurus", *state, node, f"{prefix}Density", "component0")))
                mobility = _finite(index.get(("sentaurus", *state, node, f"{prefix}Mobility", "component0")))
                jx = _finite(index.get(("sentaurus", *state, node, f"{prefix}CurrentDensity", "component0")))
                jy = _finite(index.get(("sentaurus", *state, node, f"{prefix}CurrentDensity", "component1")))
                if None not in (density, mobility, jx, jy) and density > 0 and mobility > 0:
                    gradient = current_inverted_qf_gradient(carrier, density, mobility, (jx, jy))
                    samples.append(math.hypot(*gradient))
            values.append(statistics.mean(samples) if samples else None)
        result.append({"label": f"{carrier.capitalize()} effective gradient", "values": values})
    return biases, result


def _chart_contract(question: str, takeaway: str, family: str, variant: str,
                    rows: str, fields: list[str], qa: str) -> dict:
    return {
        "question": question, "takeaway": takeaway, "family": family,
        "variant": variant, "row_grain_sufficiency": rows, "fields": fields,
        "palette_policy": {"policy": "hard two-root cap", "roots": ["blue", "orange"],
                           "non_color_distinctions": ["marker shape", "line style", "open fill"]},
        "output_paths": [], "qa_surface": qa,
    }


def _figure_specs(rows: tuple[Observation, ...], replacement: dict) -> dict[str, dict]:
    bias, potential = _series(rows, "ElectrostaticPotential")
    _, field = _series(rows, "ElectricField", magnitude=True)
    qf_bias, qf = _qf_series(rows)
    current_bias, current = _series(rows, "eCurrentDensity", magnitude=True)
    alpha_bias, alpha = _series(rows, "eAlphaAvalanche")
    _, generation = _series(rows, "ImpactIonization")
    replacement_values = [row["incremental_dex"] for row in replacement["forward"]]
    common = "40 exact topology-bias states; node means; discovery and holdout retained"
    return {
        "potential_field": {
            "title": "Potential and electric-field comparison",
            "subtitle": "Exact -1 to -20 V checkpoints; node means; potential in V and field magnitude in V/m",
            "panels": [
                {"title": "Mean absolute electrostatic potential", "x": bias, "series": potential,
                 "x_label": "Absolute applied bias (V)", "y_label": "Potential (V)"},
                {"title": "Mean electric-field magnitude", "x": bias, "series": field,
                 "x_label": "Absolute applied bias (V)", "y_label": "Field (V/m)"},
            ],
            "chart_contract": _chart_contract("How do supplied potential and field magnitudes vary by bias?",
                "The visual is descriptive; identifiability remains in the adjacent report text.",
                "Trend", "two-panel ordered-bias line", common,
                ["bias_V", "ElectrostaticPotential", "ElectricField", "solver"], "static PNG and PDF"),
        },
        "qf_gradient": {
            "title": "Current-inverted quasi-Fermi-gradient magnitude",
            "subtitle": "Sentaurus node means at exact checkpoints; V/m; mobility-gradient inference remains confounded across solvers",
            "panels": [{"title": "Effective carrier gradient", "x": qf_bias, "series": qf,
                        "x_label": "Absolute applied bias (V)", "y_label": "Effective gradient (V/m)"}],
            "chart_contract": _chart_contract("What effective qF-gradient magnitude is implied by current, density, and mobility?",
                "This diagnostic recovers a combined operator and is not a unique cross-solver formula.",
                "Trend", "ordered-bias multi-series line", common,
                ["bias_V", "current_density", "density", "mobility", "carrier"], "static PNG and PDF"),
        },
        "current_density": {
            "title": "Electron current-density magnitude",
            "subtitle": "Exact -1 to -20 V checkpoints; node means; A/m^2; solver identity uses line style and marker shape",
            "panels": [{"title": "Mean electron current-density magnitude", "x": current_bias, "series": current,
                        "x_label": "Absolute applied bias (V)", "y_label": "Current density (A/m^2)"}],
            "chart_contract": _chart_contract("How do electron current magnitudes compare across ordered bias?",
                "The chart is descriptive; signed flux and vector agreement are not interchangeable.",
                "Trend", "ordered-bias line", common,
                ["bias_V", "eCurrentDensity", "solver"], "static PNG and PDF"),
        },
        "alpha_generation": {
            "title": "Avalanche coefficient and impact generation",
            "subtitle": "Exact checkpoints; node means; alpha in m^-1 and volumetric generation in m^-3 s^-1",
            "panels": [
                {"title": "Mean absolute electron avalanche coefficient", "x": alpha_bias, "series": alpha,
                 "x_label": "Absolute applied bias (V)", "y_label": "Alpha (m^-1)"},
                {"title": "Mean absolute impact-ionization generation", "x": alpha_bias, "series": generation,
                 "x_label": "Absolute applied bias (V)", "y_label": "Generation (m^-3 s^-1)"},
            ],
            "chart_contract": _chart_contract("How do native alpha and volumetric generation vary by solver and bias?",
                "Alpha agreement is a control and does not by itself identify the driving force.",
                "Trend", "two-panel ordered-bias line", common,
                ["bias_V", "eAlphaAvalanche", "ImpactIonization", "solver"], "static PNG and PDF"),
        },
        "replacement_matrix": {
            "title": "Staged replacement contributions",
            "subtitle": "Declared seven-factor dependency order; incremental change in dex; exact arithmetic fixture for closure QA",
            "panels": [{"title": "Forward staged increments", "kind": "bar",
                        "labels": list(INVERSE_DEPENDENCIES), "values": replacement_values,
                        "x_label": "Incremental change (dex)", "y_label": "Replacement factor"}],
            "chart_contract": _chart_contract("Which declared stage contributes each arithmetic replacement increment?",
                "The decomposition validates closure mechanics; it is not a causal estimate from the supplied fixture.",
                "Decomposition & Progression", "horizontal contribution bars", "7 declared dependency stages",
                ["factor", "incremental_dex", "sequence"], "static PNG and PDF"),
        },
    }


def _report_markdown(report: dict, figures: dict) -> str:
    payload = report["payload"]
    groups = {value.value: [] for value in Identifiability}
    for item in payload["classifications"]:
        groups[item["classification"]].append(item["candidate"])
    def names(classification: Identifiability) -> str:
        values = groups[classification.value]
        return ", ".join(f"`{value}`" for value in values) if values else "none"
    paragraphs = {
        "potential_field": "The potential and field panels show ordered-bias node means for the supplied exact matrix. Read them as descriptive evidence only: agreement of means cannot establish support-local vector equivalence, and the typed candidate gates in the authoritative JSON control promotion.",
        "qf_gradient": "The current-inverted gradient uses independently supplied Sentaurus mobility, density, and current. It is useful diagnostically, but because equivalent independent mobility evidence is not available for both solvers, the cross-solver operator remains confounded rather than identified.",
        "current_density": "The current panel compares vector magnitudes at node support. It does not substitute magnitude agreement for signed edge-flux or direction agreement, so any promotion still requires the corresponding support and semantic gates.",
        "alpha_generation": "The alpha and generation panels retain their distinct units and node support. Similar alpha values cannot uniquely identify an avalanche driver, and volumetric generation is not mixed with triangle-integrated or node-mapped source.",
        "replacement_matrix": "The replacement bars show deterministic arithmetic increments through the declared dependency order. Their exact closure verifies replay bookkeeping; the synthetic decomposition is a QA control, not evidence that any physical stage caused the observed solver discrepancy.",
    }
    lines = [
        "# PN2D Minimal6 Physics Inverse Audit", "", "## Technical summary", "",
        (f"The supplied 40-state exact matrix yields **{len(groups[Identifiability.IDENTIFIED.value])} identified**, "
         f"**{len(groups[Identifiability.CONFOUNDED.value])} confounded**, "
         f"**{len(groups[Identifiability.INSUFFICIENT_DATA.value])} insufficient**, and "
         f"**{len(groups[Identifiability.REJECTED.value])} rejected** candidate conclusions under the fixed gates. "
         "These labels are identifiability results from the authoritative JSON, not causal claims; an identified label means only that the declared numerical discovery and holdout gates passed on compatible supplied support."),
        "", "## Key findings and visual evidence", "",
        f"**Identified:** {names(Identifiability.IDENTIFIED)}. **Consistent but nonunique:** {names(Identifiability.CONSISTENT_NONUNIQUE)}. **Confounded:** {names(Identifiability.CONFOUNDED)}. **Insufficient:** {names(Identifiability.INSUFFICIENT_DATA)}. **Rejected:** {names(Identifiability.REJECTED)}.", "",
    ]
    for name in FIGURE_NAMES:
        entry = next(item for item in figures["figures"] if item["name"] == name)
        lines.extend((f"### {entry['title']}", "", paragraphs[name], "",
                      f"![{entry['title']}](figures/{name}.png)", ""))
    lines.extend((
        "## Scope, data, and metric definitions", "",
        "The scope is the exact sketch-and-mirror checkpoint matrix from -1 through -20 V. Discovery comprises sketch checkpoints -1, -4, -8, -12, -16, -19, and -20 V; mirror plus all other exact checkpoints are holdout. Node observations retain raw and SI values, support, component, orientation, conversion, source path, and source hash. Median and p95 absolute log10 errors compare finite nonzero paired magnitudes; missing, below-floor, incompatible, or non-finite samples remain typed and are not replaced by zero.", "",
        "## Methodology", "",
        "The audit validates hash-bound inputs, writes canonical support tables, evaluates direct paired candidates without local scaling, reconstructs selected triangle, current-inverted, alpha, and generation identities, and runs a declared seven-stage multiplicative replacement closure control. Classification is descriptive or diagnostic unless independent factors, discovery, holdout, symmetry, direction, and magnitude gates establish identifiability. No result in this report is inferential or causal.", "",
        "## Limitations, uncertainty, and robustness", "",
        "Aggregate node means can conceal local direction and support errors. The current-inverted quasi-Fermi result is confounded by cross-solver mobility availability. Alpha inversion can be branch-ambiguous, and alpha agreement does not uniquely determine the driver. The replacement decomposition is an arithmetic robustness fixture rather than fitted physical attribution. Byte-identical dual runs within the same Python and rendering-library versions, raw-input hashes, artifact hashes, decoded PNG pixel hashes, fixed thresholds, split membership, and independent semantic replay address reproducibility but do not create evidence absent from the inputs.", "",
        "## Recommended next steps", "",
        "- Promote a formula to a production-change task only when its local support, direction, discovery, holdout, mirror, and replacement gates all pass.",
        "- Obtain independently comparable Vela and Sentaurus mobility on the same support before resolving the quasi-Fermi-gradient confounding.",
        "- Preserve volumetric, integrated, and node-mapped generation as separate quantities in any follow-up experiment.", "",
        "## Further questions", "",
        "- Which additional edge- or cell-resolved exports would distinguish candidates currently nonunique or insufficient?",
        "- Does the preferred formula remain stable across the low/high avalanche coefficient branch transition?",
        "- Do signed edge flux and reconstructed vector-current conclusions agree at the same supports?", "",
    ))
    return "\n".join(lines)



def _raw_input_ledger(bundle: InputBundle, input_roots: dict[str, str]) -> list[dict]:
    rows = []
    for logical, digest in bundle.input_hashes:
        solver, relative = logical.split(":", 1)
        input_key = ("supplemental_sentaurus_root" if solver == "supplemental"
                     else f"{solver}_root")
        rows.append({
            "logical_id": logical, "solver": solver, "relative_path": relative,
            "path": str((Path(input_roots[input_key]) / relative).resolve()),
            "sha256": digest,
        })
    for solver, root_key in (("vela", "vela_root"), ("sentaurus", "sentaurus_root"),
                             ("supplemental", "supplemental_sentaurus_root")):
        for relative in ("manifest.json", "seal.json"):
            path = (Path(input_roots[root_key]) / relative).resolve()
            rows.append({"logical_id": f"{solver}:{relative}", "solver": solver,
                         "relative_path": relative, "path": str(path),
                         "sha256": sha256(path)})
    return sorted(rows, key=lambda item: item["logical_id"])

def build_analysis_artifacts(bundle: InputBundle, out_dir: str | Path, *, phase_base: str,
                             input_roots: dict[str, str], sentaurus_version: str) -> dict:
    if phase_base != "a5524cf":
        raise ValueError("phase-base must be a5524cf")
    root = Path(out_dir)
    root.mkdir(parents=True, exist_ok=False)
    write_input_manifest(bundle, root / "input_manifest.json")
    rows = canonical_observations(bundle)
    for support in SupportKind:
        selected = [_observation_row(row) for row in rows if row.support_kind is support]
        write_csv(root / f"observations_{support.value}.csv", OBSERVATION_COLUMNS, selected)
    limits = AcceptanceThresholds()
    metrics, classifications = _candidate_metrics(bundle, rows, limits)
    write_csv(root / "candidate_metrics.csv", CANDIDATE_COLUMNS, metrics)
    write_json(root / "candidate_classifications.json", {"classifications": classifications})
    replacement = _replacement()
    replacement_rows = []
    for sequence in ("one_factor", "forward", "reverse"):
        for step, row in enumerate(replacement[sequence]):
            replacement_rows.append({
                "sequence": sequence, "step": step, "factor": row["factor"],
                "value": row["value"],
                "incremental_dex": row.get("incremental_dex", row.get("delta_dex")),
                "closure_abs_dex": None,
            })
    for name in ("forward_abs_dex", "reverse_abs_dex", "direct_abs_dex"):
        replacement_rows.append({"sequence": "closure", "step": 0, "factor": name,
                                 "value": None, "incremental_dex": None,
                                 "closure_abs_dex": replacement["closure"][name]})
    write_csv(root / "replacement_matrix.csv", REPLACEMENT_COLUMNS, replacement_rows)

    status_counts = {status.value: 0 for status in SampleStatus}
    for row in rows:
        status_counts[row.status.value] += 1
    input_manifest_hash = sha256(root / "input_manifest.json")
    raw_inputs = _raw_input_ledger(bundle, input_roots)
    report = {
        "schema": "vela.pn2d_minimal6_physics_inverse_audit.v1",
        "diagnostic_only": True,
        "phase_base": phase_base,
        "payload": {
            "input_manifest_sha256": input_manifest_hash,
            "discovery_keys": [[topology, bias] for topology, bias in bundle.discovery_keys],
            "holdout_keys": [[topology, bias] for topology, bias in bundle.holdout_keys],
            "thresholds": asdict(limits),
            "field_inventory": field_inventory(bundle),
            "sample_status_counts": status_counts,
            "candidate_metrics": metrics,
            "classifications": classifications,
            "replacement_closure": [replacement["closure"]],
            "localization_control": {
                "semantic_replay": _semantic_replay(rows, replacement),
                "input_provenance": {"input_roots": dict(input_roots),
                                     "raw_inputs": raw_inputs},
                "classification": "localization_control",
                "excluded_from_ranking": True,
            },
            "sentaurus_version": sentaurus_version,
            "production_cpp_changed": False,
        },
    }
    validate_inverse_report_v1(report)
    write_json(root / "physics_inverse_audit.json", report)
    figures = render_inverse_figures(root, _figure_specs(rows, replacement))
    write_figure_manifest(root / "figure_manifest.json", figures)
    (root / "physics_inverse_audit.md").write_text(
        _report_markdown(report, figures), encoding="utf-8", newline="\n")
    return {"report": report, "input_roots": input_roots}


def write_report_manifest(out_dir: str | Path, bundle: InputBundle, input_roots: dict[str, str]) -> dict:
    root = Path(out_dir)
    excluded = {"report_manifest.json", "verification.json", "package_manifest.json"}
    artifacts = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative not in excluded:
            artifacts[relative] = sha256(path)
    inputs = _raw_input_ledger(bundle, input_roots)

    manifest = {
        "schema": "vela.pn2d_minimal6_inverse_report_manifest.v1",
        "exclusions": sorted(excluded), "inputs": inputs, "artifacts": artifacts,
    }
    write_json(root / "report_manifest.json", manifest)
    return manifest
