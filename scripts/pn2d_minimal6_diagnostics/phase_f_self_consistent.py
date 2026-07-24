"""Phase F self-consistent Minimal6 dependency-chain comparison."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
import subprocess
from pathlib import Path
from typing import Any, Iterable


TOPOLOGIES = ("mirror", "sketch")
BIASES = tuple(range(1, 21))
INTERNAL_NODES = (1, 5)
ELEMENTARY_CHARGE_C = 1.602176634e-19
DEPENDENCY_ORDER = (
    "psi",
    "electron_qfp",
    "hole_qfp",
    "electron_density",
    "hole_density",
    "electron_mobility",
    "hole_mobility",
    "electron_directed_current",
    "hole_directed_current",
    "terminal_current",
    "impact_source",
)
THRESHOLDS = {
    "psi_max_V": 1.0e-6,
    "qfp_median_V": 0.01,
    "qfp_p95_V": 0.025,
    "density_median_dex": 0.10,
    "density_p95_dex": 0.25,
    "mobility_median_dex": 0.05,
    "mobility_p95_dex": 0.20,
    "directed_current_median_dex": 0.10,
    "directed_current_p95_dex": 0.25,
    "terminal_current_median_dex": 0.10,
    "impact_source_median_dex": 0.30,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _triangle_source_per_cm_s(rows: list[dict[str, str]]) -> float:
    """Integrate the production triangle-GSS source on physical SI support."""
    total_per_m_s = 0.0
    for row in rows:
        for local in range(3):
            volume = float(
                row[f"local_edge{local}_truncated_partial_volume_m2"]
            )
            electron = float(
                row[f"local_edge{local}_electron_source_integral_per_m_s"]
            )
            hole = float(
                row[f"local_edge{local}_hole_source_integral_per_m_s"]
            )
            if volume == 0.0 and (electron != 0.0 or hole != 0.0):
                raise ValueError("geometric-zero triangle edge has nonzero source")
            total_per_m_s += electron + hole
    return total_per_m_s * 1.0e-2


def percentile(values: Iterable[float], fraction: float) -> float:
    """Return the frozen linear-interpolation percentile."""
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("percentile requires at least one value")
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("percentile fraction must be within [0, 1]")
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def classify_log_error(
    reference: float, candidate: float
) -> tuple[str, float | None, float | None]:
    """Return a typed magnitude error without inventing values at zero."""
    reference = float(reference)
    candidate = float(candidate)
    if not math.isfinite(reference) or not math.isfinite(candidate):
        return "nonfinite", None, None
    if reference == 0.0 and candidate == 0.0:
        return "exact_zero", None, None
    if reference == 0.0 or candidate == 0.0:
        return "zero_reference_mismatch", None, None
    error = abs(math.log10(abs(candidate) / abs(reference)))
    sign = 1.0 if math.copysign(1.0, candidate) == math.copysign(1.0, reference) else 0.0
    return "valid", error, sign


def directed_current_A_per_um(
    *,
    current_density_A_per_m2: float,
    dual_length_m: float,
    candidate_node0: int,
    candidate_node1: int,
    reference_node0: int,
    reference_node1: int,
) -> float:
    candidate_pair = (int(candidate_node0), int(candidate_node1))
    reference_pair = (int(reference_node0), int(reference_node1))
    if set(candidate_pair) != set(reference_pair):
        raise ValueError("candidate and reference edges do not share one node pair")
    orientation = 1.0 if candidate_pair == reference_pair else -1.0
    return (
        orientation
        * float(current_density_A_per_m2)
        * float(dual_length_m)
        * 1.0e-6
    )


def first_failed_metric(gates: dict[str, bool]) -> str | None:
    for metric in DEPENDENCY_ORDER:
        if metric in gates and not gates[metric]:
            return metric
    return None


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _bias_tag(bias: int) -> str:
    return f"m{bias}V"


def _accepted_vela_states(
    sweep_root: Path,
) -> tuple[dict[tuple[str, int], dict[str, Any]], dict[str, Any]]:
    manifest_path = sweep_root / "sweep_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    failed = manifest.get("failed_transitions", [])
    if failed:
        first = failed[0]
        raise RuntimeError(
            "candidate sweep has a preserved failed transition: "
            f"{first.get('topology')} {first.get('target_bias_V')}"
        )
    states: dict[tuple[str, int], dict[str, Any]] = {}
    for row in manifest.get("accepted_checkpoints", []):
        if row.get("solver") != "vela":
            continue
        target = float(row["target_bias_V"])
        magnitude = int(round(abs(target)))
        if magnitude not in BIASES or abs(target + magnitude) > 1.0e-12:
            continue
        key = (str(row["topology"]), magnitude)
        if key in states:
            raise ValueError(f"duplicate Vela checkpoint {key}")
        state_path = Path(str(row["state_path"]))
        if not state_path.is_absolute():
            state_path = sweep_root / state_path
        if not state_path.is_file() or _sha256(state_path) != row["state_sha256"]:
            raise ValueError(f"candidate checkpoint is missing or hash-tampered: {key}")
        states[key] = {**row, "resolved_state_path": state_path.resolve()}
    expected = {(topology, bias) for topology in TOPOLOGIES for bias in BIASES}
    if set(states) != expected:
        missing = sorted(expected - set(states))
        raise ValueError(f"candidate sweep lacks exact checkpoints: {missing}")
    return states, manifest


def _sentaurus_manifest_rows(
    root: Path,
) -> dict[tuple[str, int], dict[str, Any]]:
    manifest = json.loads((root / "sweep_manifest.json").read_text(encoding="utf-8"))
    rows: dict[tuple[str, int], dict[str, Any]] = {}
    for row in manifest.get("accepted_checkpoints", []):
        if row.get("solver") != "sentaurus":
            continue
        target = float(row["target_bias_V"])
        bias = int(round(abs(target)))
        if bias not in BIASES or abs(target + bias) > 1.0e-12:
            continue
        rows[(str(row["topology"]), bias)] = row
    expected = {(topology, bias) for topology in TOPOLOGIES for bias in BIASES}
    if set(rows) != expected:
        raise ValueError("Sentaurus sweep reference lacks the exact 40-state lattice")
    return rows


def _candidate_state(rows: list[dict[str, str]]) -> dict[int, dict[str, float]]:
    required = {"node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("candidate state CSV lacks the Phase F state contract")
    return {
        int(row["node_id"]): {
            "psi_V": float(row["psi"]),
            "phin_V": float(row["phin"]),
            "phip_V": float(row["phip"]),
            "n_m3": float(row["electrons_m3"]),
            "p_m3": float(row["holes_m3"]),
        }
        for row in rows
    }


def _sentaurus_state(path: Path) -> dict[int, dict[str, float]]:
    rows = _read_csv(path)
    return {
        int(row["canonical_node_id"]): {
            "psi_V": float(row["ElectrostaticPotential_component0"]),
            "phin_V": float(row["eQuasiFermiPotential_component0"]),
            "phip_V": float(row["hQuasiFermiPotential_component0"]),
            "n_m3": float(row["eDensity_component0"]) * 1.0e6,
            "p_m3": float(row["hDensity_component0"]) * 1.0e6,
        }
        for row in rows
    }


def _run_operator_audit(
    *,
    executable: Path,
    sweep_root: Path,
    topology: str,
    bias: int,
    state: Path,
    raw_root: Path,
) -> Path:
    segment = bias - 1
    base = sweep_root / "vela" / topology
    work = raw_root / topology / _bias_tag(bias)
    work.mkdir(parents=True, exist_ok=True)
    canonical_state = work / "candidate_state.csv"
    state_values = _candidate_state(_read_csv(state))
    _write_csv(
        canonical_state,
        [
            {
                "node_id": node,
                "psi_V": values["psi_V"],
                "phin_V": values["phin_V"],
                "phip_V": values["phip_V"],
                "n_m3": values["n_m3"],
                "p_m3": values["p_m3"],
            }
            for node, values in sorted(state_values.items())
        ],
    )
    triangle = work / "triangle.csv"
    command = [
        str(executable),
        "--mesh",
        str(sweep_root / "inputs" / topology / "mesh.json"),
        "--doping",
        str(sweep_root / "inputs" / topology / "doping.csv"),
        "--state",
        str(canonical_state),
        "--config",
        str(base / "decks" / f"segment_{segment:02d}.json"),
        "--node-out",
        str(work / "node.csv"),
        "--edge-out",
        str(work / "edge.csv"),
        "--triangle-out",
        str(triangle),
    ]
    completed = subprocess.run(
        command, cwd=sweep_root, text=True, capture_output=True
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"operator audit failed for {topology} -{bias} V: "
            f"{completed.stderr.strip()}"
        )
    return triangle


def _candidate_element_mobility(
    triangle_rows: list[dict[str, str]],
) -> dict[tuple[int, str], float]:
    result: dict[tuple[int, str], float] = {}
    for row in triangle_rows:
        cell = int(row["cell_id"])
        for carrier in ("electron", "hole"):
            values = [
                float(row[f"local_edge{local}_{carrier}_mobility_m2_per_V_s"])
                for local in range(3)
            ]
            result[(cell, carrier)] = statistics.fmean(values)
    return result


def _endpoint_edge_rows(
    sweep_root: Path, topology: str, bias: int
) -> dict[tuple[int, int], dict[str, str]]:
    path = (
        sweep_root
        / "vela"
        / topology
        / "diagnostics"
        / f"segment_{bias - 1:02d}_sg_avalanche_edges.csv"
    )
    selected = [
        row
        for row in _read_csv(path)
        if abs(float(row["bias_V"]) + bias) <= 1.0e-12
    ]
    if not selected:
        raise ValueError(f"missing exact endpoint edge rows for {topology} -{bias} V")
    result: dict[tuple[int, int], dict[str, str]] = {}
    for row in selected:
        pair = tuple(sorted((int(row["node0"]), int(row["node1"]))))
        if pair in result:
            raise ValueError(f"duplicate candidate global edge {topology} -{bias} {pair}")
        result[pair] = row
    return result


def _stat(values: list[float]) -> dict[str, float | int]:
    if not values:
        return {"count": 0, "median": math.nan, "p95": math.nan, "maximum": math.nan}
    return {
        "count": len(values),
        "median": statistics.median(values),
        "p95": percentile(values, 0.95),
        "maximum": max(values),
    }


def _normalized_sweep_digest(
    states: dict[tuple[str, int], dict[str, Any]],
) -> str:
    payload = [
        {
            "topology": topology,
            "bias_V": -bias,
            "state_sha256": states[(topology, bias)]["state_sha256"],
            "observables": states[(topology, bias)]["observables"],
            "convergence_metadata": states[(topology, bias)].get(
                "convergence_metadata", {}
            ),
        }
        for topology in TOPOLOGIES
        for bias in BIASES
    ]
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _report(
    outcome: str,
    first_failure: str | None,
    summaries: dict[str, dict[str, float | int]],
    gates: dict[str, bool],
) -> str:
    lines = [
        "# PN2D Minimal6 Phase F self-consistent comparison",
        "",
        f"Status: `{outcome}`.",
        "",
        "The comparison covers the exact `mirror/sketch x -1..-20 V` lattice.",
        "Sentaurus directed edge current remains a box-operator reconstruction,",
        "not a native directed-edge observation.",
        "",
        f"First failed dependency metric: `{first_failure or 'none'}`.",
        "",
        "| Metric | Count | Median | P95 | Maximum | Gate |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for name in DEPENDENCY_ORDER:
        summary = summaries.get(name, {})
        lines.append(
            f"| {name} | {summary.get('count', 0)} | "
            f"{summary.get('median', math.nan):.9g} | "
            f"{summary.get('p95', math.nan):.9g} | "
            f"{summary.get('maximum', math.nan):.9g} | "
            f"{'pass' if gates.get(name, False) else 'fail'} |"
        )
    lines.extend(
        [
            "",
            "Thresholds were frozen in the revised Phase B-G plan before this run.",
            "Zero/reference-missing samples remain typed and are never assigned an",
            "invented dex value.",
            "",
        ]
    )
    return "\n".join(lines)


def run_phase_f(
    *,
    candidate_sweep_root: str | Path,
    sentaurus_sweep_root: str | Path,
    inverse_inputs_root: str | Path,
    phase_c_root: str | Path,
    phase_d_root: str | Path,
    operator_audit: str | Path,
    output_root: str | Path,
) -> dict[str, Any]:
    candidate_root = Path(candidate_sweep_root).resolve()
    sentaurus_root = Path(sentaurus_sweep_root).resolve()
    inverse_root = Path(inverse_inputs_root).resolve()
    phase_c = Path(phase_c_root).resolve()
    phase_d = Path(phase_d_root).resolve()
    audit = Path(operator_audit).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)

    candidate_rows, _ = _accepted_vela_states(candidate_root)
    sentaurus_rows = _sentaurus_manifest_rows(sentaurus_root)
    phase_d_rows = {
        (
            row["topology"],
            int(round(abs(float(row["bias_V"])))),
            int(row["cell_id"]),
            row["carrier"],
        ): row
        for row in _read_csv(phase_d / "native_element_decomposition.csv")
        if row["status"] == "valid"
    }
    edge_reference = [
        row
        for row in _read_csv(phase_c / "stage_edge_samples.csv")
        if row["stage"] == "vela_baseline"
        and row["status"] == "valid"
        and float(row["reference_A_per_um"]) != 0.0
    ]
    edges_by_state: dict[tuple[str, int], list[dict[str, str]]] = {}
    for row in edge_reference:
        key = (row["topology"], int(round(abs(float(row["bias_V"])))))
        edges_by_state.setdefault(key, []).append(row)

    state_output: list[dict[str, Any]] = []
    mobility_output: list[dict[str, Any]] = []
    current_output: list[dict[str, Any]] = []
    terminal_output: list[dict[str, Any]] = []
    raw_root = output / "raw"

    for topology in TOPOLOGIES:
        for bias in BIASES:
            key = (topology, bias)
            candidate_record = candidate_rows[key]
            candidate = _candidate_state(
                _read_csv(candidate_record["resolved_state_path"])
            )
            sentaurus = _sentaurus_state(
                inverse_root / "sentaurus" / "states" / topology / f"m{bias}V.csv"
            )
            for node in INTERNAL_NODES:
                row: dict[str, Any] = {
                    "topology": topology,
                    "bias_V": -bias,
                    "node_id": node,
                }
                for quantity in ("psi_V", "phin_V", "phip_V"):
                    row[f"sentaurus_{quantity}"] = sentaurus[node][quantity]
                    row[f"vela_{quantity}"] = candidate[node][quantity]
                    row[f"{quantity}_abs_error_V"] = abs(
                        candidate[node][quantity] - sentaurus[node][quantity]
                    )
                for quantity in ("n_m3", "p_m3"):
                    row[f"sentaurus_{quantity}"] = sentaurus[node][quantity]
                    row[f"vela_{quantity}"] = candidate[node][quantity]
                    status, error, _ = classify_log_error(
                        sentaurus[node][quantity], candidate[node][quantity]
                    )
                    row[f"{quantity}_status"] = status
                    row[f"{quantity}_abs_error_dex"] = "" if error is None else error
                state_output.append(row)

            triangle = _run_operator_audit(
                executable=audit,
                sweep_root=candidate_root,
                topology=topology,
                bias=bias,
                state=candidate_record["resolved_state_path"],
                raw_root=raw_root,
            )
            candidate_mobility = _candidate_element_mobility(_read_csv(triangle))
            for cell in range(4):
                for carrier in ("electron", "hole"):
                    reference = phase_d_rows[(topology, bias, cell, carrier)]
                    reference_value = float(reference["sentaurus_native_final_m2_per_Vs"])
                    candidate_value = candidate_mobility[(cell, carrier)]
                    status, error, _ = classify_log_error(
                        reference_value, candidate_value
                    )
                    mobility_output.append(
                        {
                            "topology": topology,
                            "bias_V": -bias,
                            "cell_id": cell,
                            "carrier": carrier,
                            "sentaurus_native_mobility_m2_per_Vs": reference_value,
                            "vela_candidate_cell_average_mobility_m2_per_Vs": candidate_value,
                            "absolute_log10_error_dex": "" if error is None else error,
                            "status": status,
                        }
                    )

            candidate_edges = _endpoint_edge_rows(candidate_root, topology, bias)
            for reference in edges_by_state[key]:
                node0 = int(reference["node0"])
                node1 = int(reference["node1"])
                pair = tuple(sorted((node0, node1)))
                edge = candidate_edges.get(pair)
                if edge is None:
                    raise ValueError(
                        f"candidate lacks active reconstructed box edge {key} {pair}"
                    )
                carrier = reference["carrier"]
                raw_particle_flux = float(
                    edge[f"{carrier}_raw_signed_flux_proxy"]
                )
                carrier_orientation = 1.0 if carrier == "electron" else -1.0
                current_density = (
                    carrier_orientation * ELEMENTARY_CHARGE_C
                    * raw_particle_flux * 1.0e4
                )
                candidate_current = directed_current_A_per_um(
                    current_density_A_per_m2=current_density,
                    dual_length_m=float(edge["edge_couple_m"]),
                    candidate_node0=int(edge["node0"]),
                    candidate_node1=int(edge["node1"]),
                    reference_node0=node0,
                    reference_node1=node1,
                )
                reference_current = float(reference["reference_A_per_um"])
                status, error, sign = classify_log_error(
                    reference_current, candidate_current
                )
                current_output.append(
                    {
                        "topology": topology,
                        "bias_V": -bias,
                        "carrier": carrier,
                        "node0": node0,
                        "node1": node1,
                        "sentaurus_box_reconstructed_A_per_um": reference_current,
                        "vela_directed_A_per_um": candidate_current,
                        "absolute_log10_error_dex": "" if error is None else error,
                        "sign_agreement": "" if sign is None else sign,
                        "status": status,
                    }
                )

            sentaurus_endpoint = sentaurus_rows[key]["observables"]
            candidate_endpoint = candidate_record["observables"]
            reference_current = float(sentaurus_endpoint["anode_current_A_per_um"])
            candidate_current = float(candidate_endpoint["anode_current_A_per_um"])
            current_status, current_error, current_sign = classify_log_error(
                reference_current, candidate_current
            )
            reference_source = float(
                sentaurus_endpoint["native_source_integral_s_inv_per_cm"]
            )
            candidate_source = _triangle_source_per_cm_s(
                _read_csv(triangle)
            )
            source_status, source_error, _ = classify_log_error(
                reference_source, candidate_source
            )
            terminal_output.append(
                {
                    "topology": topology,
                    "bias_V": -bias,
                    "sentaurus_terminal_current_A_per_um": reference_current,
                    "vela_terminal_current_A_per_um": candidate_current,
                    "terminal_absolute_log10_error_dex": (
                        "" if current_error is None else current_error
                    ),
                    "terminal_sign_agreement": (
                        "" if current_sign is None else current_sign
                    ),
                    "terminal_status": current_status,
                    "sentaurus_impact_source_s_inv_per_cm": reference_source,
                    "vela_impact_source_s_inv_per_cm": candidate_source,
                    "impact_absolute_log10_error_dex": (
                        "" if source_error is None else source_error
                    ),
                    "impact_status": source_status,
                }
            )

    if len(state_output) != 80:
        raise ValueError("Phase F requires 80 internal-node comparison rows")
    if len(mobility_output) != 320:
        raise ValueError("Phase F requires 320 carrier-element mobility rows")
    if len(current_output) != 400:
        raise ValueError("Phase F requires 400 active directed-current rows")
    if len(terminal_output) != 40:
        raise ValueError("Phase F requires 40 terminal/source rows")

    _write_csv(output / "state_node_comparison.csv", state_output)
    _write_csv(output / "mobility_element_comparison.csv", mobility_output)
    _write_csv(output / "directed_edge_current_comparison.csv", current_output)
    _write_csv(output / "terminal_source_comparison.csv", terminal_output)

    summaries = {
        "psi": _stat([float(row["psi_V_abs_error_V"]) for row in state_output]),
        "electron_qfp": _stat(
            [float(row["phin_V_abs_error_V"]) for row in state_output]
        ),
        "hole_qfp": _stat(
            [float(row["phip_V_abs_error_V"]) for row in state_output]
        ),
        "electron_density": _stat(
            [
                float(row["n_m3_abs_error_dex"])
                for row in state_output
                if row["n_m3_status"] == "valid"
            ]
        ),
        "hole_density": _stat(
            [
                float(row["p_m3_abs_error_dex"])
                for row in state_output
                if row["p_m3_status"] == "valid"
            ]
        ),
    }
    for carrier in ("electron", "hole"):
        summaries[f"{carrier}_mobility"] = _stat(
            [
                float(row["absolute_log10_error_dex"])
                for row in mobility_output
                if row["carrier"] == carrier and row["status"] == "valid"
            ]
        )
        summaries[f"{carrier}_directed_current"] = _stat(
            [
                float(row["absolute_log10_error_dex"])
                for row in current_output
                if row["carrier"] == carrier and row["status"] == "valid"
            ]
        )
    summaries["terminal_current"] = _stat(
        [
            float(row["terminal_absolute_log10_error_dex"])
            for row in terminal_output
            if row["terminal_status"] == "valid"
        ]
    )
    summaries["impact_source"] = _stat(
        [
            float(row["impact_absolute_log10_error_dex"])
            for row in terminal_output
            if row["impact_status"] == "valid"
        ]
    )

    all_current_signs = {
        carrier: [
            float(row["sign_agreement"])
            for row in current_output
            if row["carrier"] == carrier and row["status"] == "valid"
        ]
        for carrier in ("electron", "hole")
    }
    terminal_signs = [
        float(row["terminal_sign_agreement"])
        for row in terminal_output
        if row["terminal_status"] == "valid"
    ]
    order_failures: list[dict[str, Any]] = []
    for topology in TOPOLOGIES:
        rows = sorted(
            (row for row in terminal_output if row["topology"] == topology),
            key=lambda row: abs(float(row["bias_V"])),
        )
        for previous, current in zip(rows, rows[1:]):
            previous_abs = abs(float(previous["vela_terminal_current_A_per_um"]))
            current_abs = abs(float(current["vela_terminal_current_A_per_um"]))
            if current_abs < previous_abs * (1.0 - 1.0e-12):
                order_failures.append(
                    {
                        "topology": topology,
                        "previous_bias_V": previous["bias_V"],
                        "bias_V": current["bias_V"],
                        "previous_abs_current_A_per_um": previous_abs,
                        "abs_current_A_per_um": current_abs,
                    }
                )

    gates = {
        "psi": summaries["psi"]["maximum"] <= THRESHOLDS["psi_max_V"],
        "electron_qfp": (
            summaries["electron_qfp"]["median"] <= THRESHOLDS["qfp_median_V"]
            and summaries["electron_qfp"]["p95"] <= THRESHOLDS["qfp_p95_V"]
        ),
        "hole_qfp": (
            summaries["hole_qfp"]["median"] <= THRESHOLDS["qfp_median_V"]
            and summaries["hole_qfp"]["p95"] <= THRESHOLDS["qfp_p95_V"]
        ),
        "electron_density": (
            summaries["electron_density"]["median"]
            <= THRESHOLDS["density_median_dex"]
            and summaries["electron_density"]["p95"]
            <= THRESHOLDS["density_p95_dex"]
        ),
        "hole_density": (
            summaries["hole_density"]["median"]
            <= THRESHOLDS["density_median_dex"]
            and summaries["hole_density"]["p95"]
            <= THRESHOLDS["density_p95_dex"]
        ),
        "electron_mobility": (
            summaries["electron_mobility"]["median"]
            <= THRESHOLDS["mobility_median_dex"]
            and summaries["electron_mobility"]["p95"]
            <= THRESHOLDS["mobility_p95_dex"]
        ),
        "hole_mobility": (
            summaries["hole_mobility"]["median"]
            <= THRESHOLDS["mobility_median_dex"]
            and summaries["hole_mobility"]["p95"]
            <= THRESHOLDS["mobility_p95_dex"]
        ),
        "electron_directed_current": (
            summaries["electron_directed_current"]["median"]
            <= THRESHOLDS["directed_current_median_dex"]
            and summaries["electron_directed_current"]["p95"]
            <= THRESHOLDS["directed_current_p95_dex"]
            and len(all_current_signs["electron"]) == 200
            and all(value == 1.0 for value in all_current_signs["electron"])
        ),
        "hole_directed_current": (
            summaries["hole_directed_current"]["median"]
            <= THRESHOLDS["directed_current_median_dex"]
            and summaries["hole_directed_current"]["p95"]
            <= THRESHOLDS["directed_current_p95_dex"]
            and len(all_current_signs["hole"]) == 200
            and all(value == 1.0 for value in all_current_signs["hole"])
        ),
        "terminal_current": (
            summaries["terminal_current"]["median"]
            <= THRESHOLDS["terminal_current_median_dex"]
            and len(terminal_signs) == 40
            and all(value == 1.0 for value in terminal_signs)
            and not order_failures
        ),
        "impact_source": (
            summaries["impact_source"]["median"]
            <= THRESHOLDS["impact_source_median_dex"]
            and summaries["impact_source"]["count"] == 40
        ),
    }
    failure = first_failed_metric(gates)
    outcome = "parity_passed" if failure is None else "model_difference"

    summary_rows = [
        {
            "metric": metric,
            **summaries[metric],
            "gate_passed": int(gates[metric]),
        }
        for metric in DEPENDENCY_ORDER
    ]
    _write_csv(output / "summary.csv", summary_rows)
    if order_failures:
        _write_csv(output / "branch_order_failures.csv", order_failures)

    report_text = _report(outcome, failure, summaries, gates)
    (output / "report.md").write_text(
        report_text, encoding="ascii", errors="strict", newline="\n"
    )
    output_names = [
        "state_node_comparison.csv",
        "mobility_element_comparison.csv",
        "directed_edge_current_comparison.csv",
        "terminal_source_comparison.csv",
        "summary.csv",
        "report.md",
    ]
    if order_failures:
        output_names.append("branch_order_failures.csv")
    manifest = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "pn2d_minimal6_phase_f_self_consistent",
        "outcome": {
            "status": outcome,
            "first_failed_metric": failure,
            "gates": gates,
            "branch_order_failure_count": len(order_failures),
        },
        "thresholds": THRESHOLDS,
        "contracts": {
            "state_count": 40,
            "internal_node_row_count": len(state_output),
            "mobility_row_count": len(mobility_output),
            "directed_current_row_count": len(current_output),
            "terminal_source_row_count": len(terminal_output),
            "sentaurus_current_support": "box_operator_reconstruction",
        },
        "summaries": summaries,
        "inputs": {
            "candidate_sweep_normalized_sha256": _normalized_sweep_digest(
                candidate_rows
            ),
            "operator_audit_sha256": _sha256(audit),
            "sentaurus_sweep_manifest_sha256": _sha256(
                sentaurus_root / "sweep_manifest.json"
            ),
            "inverse_sentaurus_manifest_sha256": _sha256(
                inverse_root / "sentaurus" / "manifest.json"
            ),
            "phase_c_manifest_sha256": _sha256(phase_c / "manifest.json"),
            "phase_d_manifest_sha256": _sha256(phase_d / "manifest.json"),
        },
        "outputs": {name: _sha256(output / name) for name in output_names},
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    return manifest
