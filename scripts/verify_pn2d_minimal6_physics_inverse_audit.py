#!/usr/bin/env python3
"""Raw-input semantic verifier for the Minimal6 inverse-audit package.

This module does not import the Task 8 report, plot, or diagnose modules.  It
reconstructs expected artifacts from the sealed Task 2 input bundle and uses
Task 4-7 diagnostic primitives only for independently specified calculations.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from PIL import Image

from scripts.pn2d_minimal6_diagnostics.inverse_avalanche import (
    impact_generation, invert_van_overstraeten_alpha,
)
from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (
    AcceptanceThresholds, Identifiability, Observation, SampleStatus, SupportKind,
)
from scripts.pn2d_minimal6_diagnostics.inverse_fields import triangle_gradient
from scripts.pn2d_minimal6_diagnostics.inverse_inputs import (
    InputBundle, canonical_observations, field_inventory, load_input_bundle,
)
from scripts.pn2d_minimal6_diagnostics.inverse_replacements import (
    INVERSE_DEPENDENCIES, run_replacement_matrix,
)
from scripts.pn2d_minimal6_diagnostics.inverse_transport import current_inverted_qf_gradient


REPORT_EXCLUSIONS = {"report_manifest.json", "verification.json", "package_manifest.json"}
PACKAGE_EXCLUSIONS = {"package_manifest.json"}
FIGURES = ("potential_field", "qf_gradient", "current_density",
           "alpha_generation", "replacement_matrix")
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


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pixel_sha256(path: Path) -> str:
    with Image.open(path) as image:
        return hashlib.sha256(image.convert("RGBA").tobytes()).hexdigest()


def _load_json(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"),
                          parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid JSON: {path.name}") from error


def _write_json_bytes(payload: object) -> bytes:
    return (json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def _format_csv(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("expected CSV value is non-finite")
        return format(value, ".17g")
    if hasattr(value, "value"):
        return str(value.value)
    return str(value)


def _read_csv(path: Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != columns:
            raise ValueError(f"{path.name} column contract mismatch")
        return list(reader)


def _raw_ledger(bundle: InputBundle, roots: dict[str, str]) -> list[dict]:
    ledger = []
    for logical, digest in bundle.input_hashes:
        solver, relative = logical.split(":", 1)
        root_key = "supplemental_sentaurus_root" if solver == "supplemental" else f"{solver}_root"
        ledger.append({
            "logical_id": logical, "solver": solver, "relative_path": relative,
            "path": str((Path(roots[root_key]) / relative).resolve()), "sha256": digest,
        })
    for solver, root_key in (("vela", "vela_root"), ("sentaurus", "sentaurus_root"),
                             ("supplemental", "supplemental_sentaurus_root")):
        for relative in ("manifest.json", "seal.json"):
            path = (Path(roots[root_key]) / relative).resolve()
            ledger.append({"logical_id": f"{solver}:{relative}", "solver": solver,
                           "relative_path": relative, "path": str(path),
                           "sha256": _sha256(path)})
    return sorted(ledger, key=lambda item: item["logical_id"])


def _load_raw_bundle(report: dict) -> tuple[InputBundle, dict[str, str], list[dict]]:
    try:
        provenance = report["payload"]["localization_control"]["input_provenance"]
        roots = provenance["input_roots"]
        persisted_ledger = provenance["raw_inputs"]
    except (KeyError, TypeError) as error:
        raise ValueError("authoritative raw-input provenance is missing") from error
    required = {"vela_root", "sentaurus_root", "supplemental_sentaurus_root"}
    if not isinstance(roots, dict) or set(roots) != required:
        raise ValueError("authoritative input roots mismatch")
    canonical_roots = {key: str(Path(value).resolve()) for key, value in roots.items()}
    if canonical_roots != roots:
        raise ValueError("authoritative input roots are not exact resolved paths")
    bundle = load_input_bundle(roots["vela_root"], roots["sentaurus_root"],
                               roots["supplemental_sentaurus_root"])
    expected_ledger = _raw_ledger(bundle, roots)
    if persisted_ledger != expected_ledger:
        raise ValueError("authoritative raw-input ledger semantic mismatch")
    for item in expected_ledger:
        path = Path(item["path"])
        if not path.is_file() or _sha256(path) != item["sha256"]:
            raise ValueError(f"raw input hash mismatch: {item['logical_id']}")
    return bundle, roots, expected_ledger


def _verify_report_manifest(root: Path, raw_ledger: list[dict]) -> dict:
    manifest = _load_json(root / "report_manifest.json")
    if not isinstance(manifest, dict) or manifest.get("schema") != "vela.pn2d_minimal6_inverse_report_manifest.v1":
        raise ValueError("report manifest schema mismatch")
    if set(manifest.get("exclusions", ())) != REPORT_EXCLUSIONS:
        raise ValueError("report manifest exclusions mismatch")
    if manifest.get("inputs") != raw_ledger:
        raise ValueError("report manifest raw-input ledger mismatch")
    actual = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative not in REPORT_EXCLUSIONS:
            actual[relative] = _sha256(path)
    if manifest.get("artifacts") != actual:
        raise ValueError("report artifact hash mismatch")
    return manifest


def _verify_package_manifest(root: Path, raw_ledger: list[dict] | None = None) -> None:
    path = root / "package_manifest.json"
    if not path.exists():
        return
    manifest = _load_json(path)
    if not isinstance(manifest, dict) or manifest.get("schema") != "vela.pn2d_minimal6_inverse_package_manifest.v1":
        raise ValueError("package manifest schema mismatch")
    if set(manifest.get("exclusions", ())) != PACKAGE_EXCLUSIONS:
        raise ValueError("package manifest exclusions mismatch")
    for item in manifest.get("raw_inputs", []):
        input_path = Path(item["path"])
        if not input_path.is_file() or _sha256(input_path) != item.get("sha256"):
            raise ValueError(f"package raw input hash mismatch: {item.get('logical_id')}")
    if raw_ledger is not None and manifest.get("raw_inputs") != raw_ledger:
        raise ValueError("package raw-input ledger mismatch")
    actual = {}
    for member in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = member.relative_to(root).as_posix()
        if relative not in PACKAGE_EXCLUSIONS:
            actual[relative] = _sha256(member)
    if manifest.get("artifacts") != actual:
        raise ValueError("package manifest hash mismatch")


def _observation_expected(row: Observation) -> dict[str, str]:
    values = {
        "solver": row.solver, "topology": row.topology, "bias_V": row.bias_V,
        "support_kind": row.support_kind.value, "support_id": row.support_id,
        "quantity": row.quantity, "component": row.component, "raw_value": row.raw_value,
        "raw_unit": row.raw_unit, "value_si": row.value_si, "unit_si": row.unit_si,
        "coordinate_frame": row.coordinate_frame, "orientation": row.orientation,
        "conversion": row.conversion, "status": row.status.value,
        "source_path": row.source_path, "source_sha256": row.source_sha256,
    }
    return {column: _format_csv(values[column]) for column in OBSERVATION_COLUMNS}


def _verify_observations(root: Path, rows: tuple[Observation, ...]) -> None:
    for support in SupportKind:
        expected = [_observation_expected(row) for row in rows if row.support_kind is support]
        actual = _read_csv(root / f"observations_{support.value}.csv", OBSERVATION_COLUMNS)
        if actual != expected:
            raise ValueError(f"observations_{support.value}.csv raw semantic mismatch")


def _index(rows: tuple[Observation, ...]) -> dict[tuple, Observation]:
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
    lower, upper = math.floor(position), math.ceil(position)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _paired_errors(rows: tuple[Observation, ...], quantities: set[str],
                   split_keys: set[tuple[str, float]]) -> list[float]:
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
        errors.append(abs(math.log10(max(abs(first), 1.0e-300) /
                                     max(abs(second), 1.0e-300))))
    return errors


def _expected_candidates(bundle: InputBundle, rows: tuple[Observation, ...]) -> tuple[list[dict], list[dict]]:
    limits = AcceptanceThresholds()
    specifications = (
        ("potential_field_direct", "potential_and_field", "", {"ElectrostaticPotential", "ElectricField"}, limits.gradient_median_abs_dex),
        ("current_density_direct", "current_density", "both", {"eCurrentDensity", "hCurrentDensity"}, limits.gradient_median_abs_dex),
        ("alpha_generation_direct", "alpha_and_generation", "both", {"eAlphaAvalanche", "hAlphaAvalanche", "ImpactIonization"}, limits.local_generation_abs_dex),
    )
    split_sets = {"discovery": set(bundle.discovery_keys), "holdout": set(bundle.holdout_keys),
                  "combined": set(bundle.common_keys)}
    metrics, final_classes = [], {}
    for candidate, quantity, carrier, fields, gate in specifications:
        split_classes = {}
        for split in ("discovery", "holdout", "combined"):
            errors = _paired_errors(rows, fields, split_sets[split])
            median = statistics.median(errors) if errors else None
            p95 = _percentile(errors, 0.95)
            classification = (Identifiability.IDENTIFIED if median is not None and p95 is not None
                              and median <= gate and p95 <= limits.gradient_p95_abs_dex
                              else Identifiability.INSUFFICIENT_DATA if not errors
                              else Identifiability.REJECTED)
            split_classes[split] = classification
            metrics.append({
                "candidate": candidate, "quantity": quantity, "carrier": carrier,
                "split": split, "topology": "all", "bias_V": None,
                "support_kind": SupportKind.NODE.value, "valid_count": len(errors),
                "median_abs_error": median, "p95_abs_error": p95,
                "median_angle_deg": None, "classification": classification.value,
            })
        final = split_classes["combined"]
        if final is Identifiability.IDENTIFIED and split_classes["holdout"] is not Identifiability.IDENTIFIED:
            final = Identifiability.REJECTED
        final_classes[candidate] = final
    metrics.append({
        "candidate": "current_inverted_qf_gradient", "quantity": "qf_gradient",
        "carrier": "both", "split": "combined", "topology": "all", "bias_V": None,
        "support_kind": SupportKind.NODE.value, "valid_count": 0,
        "median_abs_error": None, "p95_abs_error": None, "median_angle_deg": None,
        "classification": Identifiability.CONFOUNDED.value,
    })
    final_classes["current_inverted_qf_gradient"] = Identifiability.CONFOUNDED
    classifications = []
    for candidate in sorted(final_classes):
        classification = final_classes[candidate]
        classifications.append({
            "candidate": candidate, "classification": classification.value,
            "claim_type": "identifiability",
            "reason": ("combined and holdout numerical gates passed without local fitting"
                       if classification is Identifiability.IDENTIFIED else
                       "mobility and gradient are not independently available for both solvers"
                       if classification is Identifiability.CONFOUNDED else
                       "no compatible finite paired support was available"
                       if classification is Identifiability.INSUFFICIENT_DATA else
                       "one or more discovery, holdout, or combined gates failed"),
        })
    return metrics, classifications


def _candidate_csv(metrics: list[dict]) -> list[dict[str, str]]:
    return [{column: _format_csv(row.get(column)) for column in CANDIDATE_COLUMNS}
            for row in metrics]


def _expected_replacement() -> dict:
    def operand(factor: str, value: float) -> dict:
        return {"factor": factor, "value": value, "status": SampleStatus.VALID.value,
                "support_kind": SupportKind.INTEGRATED.value, "support_id": "global",
                "unit_si": "dimensionless", "carrier": None, "topology": "all", "bias_V": -20.0}
    baseline = {factor: operand(factor, 1.0) for factor in INVERSE_DEPENDENCIES}
    replacement = {factor: operand(factor, 1.0 + (index + 1) / 100.0)
                   for index, factor in enumerate(INVERSE_DEPENDENCIES)}
    return run_replacement_matrix(baseline, replacement,
                                  direct_target=math.prod(row["value"] for row in replacement.values()))


def _replacement_csv(replacement: dict) -> list[dict[str, str]]:
    rows = []
    for sequence in ("one_factor", "forward", "reverse"):
        for step, item in enumerate(replacement[sequence]):
            row = {"sequence": sequence, "step": step, "factor": item["factor"],
                   "value": item["value"],
                   "incremental_dex": item.get("incremental_dex", item.get("delta_dex")),
                   "closure_abs_dex": None}
            rows.append({column: _format_csv(row.get(column)) for column in REPLACEMENT_COLUMNS})
    for name in ("forward_abs_dex", "reverse_abs_dex", "direct_abs_dex"):
        row = {"sequence": "closure", "step": 0, "factor": name, "value": None,
               "incremental_dex": None, "closure_abs_dex": replacement["closure"][name]}
        rows.append({column: _format_csv(row.get(column)) for column in REPLACEMENT_COLUMNS})
    return rows


def _expected_semantic(rows: tuple[Observation, ...], replacement: dict) -> dict:
    index = _index(rows)
    states = sorted({(row.topology, row.bias_V) for row in rows if row.solver == "sentaurus"})
    state = states[0]
    nodes = sorted({str(row.support_id) for row in rows if row.solver == "sentaurus"
                    and (row.topology, row.bias_V) == state and row.quantity == "coordinate"})
    replay: dict[str, object] = {"state": [state[0], state[1]]}
    if len(nodes) >= 3:
        points, potentials = [], []
        for node in nodes[:3]:
            points.append([_finite(index.get(("sentaurus", *state, node, "coordinate", "x"))),
                           _finite(index.get(("sentaurus", *state, node, "coordinate", "y")))])
            potentials.append(_finite(index.get(("sentaurus", *state, node,
                                                  "ElectrostaticPotential", "component0"))))
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
            "value_V_per_m": list(gradient)}
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
        replay["inverse_alpha"] = {"status": status, "alpha_m_inv": alpha_n,
                                   "prefactor_m_inv": prefactor,
                                   "critical_field_V_per_m": 1.0, "gamma": 1.0,
                                   "field_V_per_m": field}
    else:
        replay["inverse_alpha"] = {"status": "insufficient_data"}
    if None not in (alpha_n, alpha_p, jx, jy, hjx, hjy):
        generation = impact_generation(alpha_n, (jx, jy), alpha_p, (hjx, hjy))
        replay["generation"] = {"status": "valid", "alpha_n_m_inv": alpha_n,
                                "alpha_p_m_inv": alpha_p, "jn_A_per_m2": [jx, jy],
                                "jp_A_per_m2": [hjx, hjy], "value_m3_s_inv": generation}
    else:
        replay["generation"] = {"status": "insufficient_data"}
    native = [_finite(index.get(("sentaurus", *state, current, "ImpactIonization", "component0")))
              for current in nodes[:3]]
    triangle = replay["triangle_gradient"]
    if triangle["status"] == "valid" and len(native) == 3 and all(value is not None for value in native):
        (x0, y0), (x1, y1), (x2, y2) = triangle["points_m"]
        area = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        replay["generation_support_integral"] = {
            "status": "valid", "support_kind": "triangle_from_raw_nodes",
            "support_ids": nodes[:3], "native_generation_m3_s_inv": native,
            "area_m2": area, "depth_m": 0.01,
            "integral_s_inv": statistics.mean(native) * area * 0.01}
    else:
        replay["generation_support_integral"] = {"status": "insufficient_data"}
    replay["replacement_closure"] = replacement["closure"]
    return replay


def _sentaurus_version(roots: dict[str, str]) -> str:
    versions = set()
    for key in ("sentaurus_root", "supplemental_sentaurus_root"):
        manifest = _load_json(Path(roots[key]) / "manifest.json")
        value = manifest.get("provenance", {}).get("sentaurus_version", manifest.get("sentaurus_version"))
        if isinstance(value, str) and value:
            versions.add(value)
    if len(versions) > 1:
        raise ValueError("raw inputs declare inconsistent Sentaurus versions")
    return next(iter(versions)) if versions else "not_declared_by_input_contract"


def _verify_raw_semantics(root: Path, report: dict, bundle: InputBundle,
                          roots: dict[str, str], rows: tuple[Observation, ...]) -> list[str]:
    payload = report["payload"]
    if report.get("schema") != "vela.pn2d_minimal6_physics_inverse_audit.v1" or report.get("diagnostic_only") is not True:
        raise ValueError("authoritative report schema mismatch")
    if report.get("phase_base") != "a5524cf" or payload.get("production_cpp_changed") is not False:
        raise ValueError("production baseline guard mismatch")
    if payload.get("discovery_keys") != [[a, b] for a, b in bundle.discovery_keys] or payload.get("holdout_keys") != [[a, b] for a, b in bundle.holdout_keys]:
        raise ValueError("raw discovery/holdout reconstruction mismatch")
    if payload.get("thresholds") != asdict(AcceptanceThresholds()):
        raise ValueError("raw threshold reconstruction mismatch")
    if payload.get("field_inventory") != field_inventory(bundle):
        raise ValueError("raw field inventory reconstruction mismatch")
    counts = {status.value: 0 for status in SampleStatus}
    for row in rows:
        counts[row.status.value] += 1
    if payload.get("sample_status_counts") != counts:
        raise ValueError("raw sample-status reconstruction mismatch")
    _verify_observations(root, rows)
    metrics, classifications = _expected_candidates(bundle, rows)
    if payload.get("candidate_metrics") != metrics:
        raise ValueError("authoritative candidate metric reconstruction mismatch")
    if payload.get("classifications") != classifications:
        raise ValueError("authoritative classification reconstruction mismatch")
    if _read_csv(root / "candidate_metrics.csv", CANDIDATE_COLUMNS) != _candidate_csv(metrics):
        raise ValueError("candidate_metrics.csv raw semantic mismatch")
    persisted = _load_json(root / "candidate_classifications.json")
    if persisted != {"classifications": classifications}:
        raise ValueError("candidate_classifications.json raw semantic mismatch")
    replacement = _expected_replacement()
    if payload.get("replacement_closure") != [replacement["closure"]]:
        raise ValueError("authoritative replacement closure mismatch")
    if _read_csv(root / "replacement_matrix.csv", REPLACEMENT_COLUMNS) != _replacement_csv(replacement):
        raise ValueError("replacement_matrix.csv raw semantic mismatch")
    expected_replay = _expected_semantic(rows, replacement)
    if payload["localization_control"].get("semantic_replay") != expected_replay:
        raise ValueError("raw semantic replay mismatch")
    if payload.get("sentaurus_version") != _sentaurus_version(roots):
        raise ValueError("raw Sentaurus version reconstruction mismatch")
    expected_input_manifest = {
        "schema": "vela.pn2d_minimal6_inverse_inputs.v1",
        "common_keys": [[a, b] for a, b in bundle.common_keys],
        "discovery_keys": [[a, b] for a, b in bundle.discovery_keys],
        "holdout_keys": [[a, b] for a, b in bundle.holdout_keys],
        "field_inventory": field_inventory(bundle),
        "executable_sha256": list(bundle.executable_hashes),
        "tracked_source_sha256": dict(bundle.tracked_source_hashes),
        "input_sha256": dict(bundle.input_hashes),
    }
    if _load_json(root / "input_manifest.json") != expected_input_manifest:
        raise ValueError("input_manifest.json raw semantic mismatch")
    if payload.get("input_manifest_sha256") != _sha256(root / "input_manifest.json"):
        raise ValueError("authoritative input-manifest binding mismatch")
    return ["raw_observations", "raw_candidate_metrics", "raw_classifications",
            "raw_replacement_sequences", "raw_triangle_gradient",
            "raw_current_inverted_gradient", "raw_inverse_alpha",
            "raw_local_generation", "raw_generation_support_integral"]


def _verify_figures(root: Path) -> list[str]:
    manifest = _load_json(root / "figure_manifest.json")
    entries = manifest.get("figures") if isinstance(manifest, dict) else None
    if not isinstance(entries, list) or [item.get("name") for item in entries] != list(FIGURES):
        raise ValueError("fixed figure set mismatch")
    renderer = manifest.get("renderer", {})
    if renderer.get("backend") != "Pillow PNG + ReportLab invariant PDF" or "Agg" in json.dumps(renderer):
        raise ValueError("figure renderer metadata mismatch")
    required = {"question", "takeaway", "family", "variant", "row_grain_sufficiency",
                "fields", "palette_policy", "output_paths", "qa_surface"}
    for item in entries:
        name = item["name"]
        png, pdf = root / "figures" / f"{name}.png", root / "figures" / f"{name}.pdf"
        if _sha256(png) != item["png_sha256"] or _sha256(pdf) != item["pdf_sha256"]:
            raise ValueError("figure file hash mismatch")
        if _pixel_sha256(png) != item["png_pixel_sha256"]:
            raise ValueError("figure decoded-pixel hash mismatch")
        if set(item["chart_contract"]) != required:
            raise ValueError("figure chart-contract mismatch")
    return ["five_figure_pairs", "png_pixel_hashes", "truthful_renderer_metadata"]


def _write_or_validate(path: Path, payload: object, label: str) -> None:
    encoded = _write_json_bytes(payload)
    if path.exists():
        if path.read_bytes() != encoded:
            raise ValueError(f"{label} is not byte-stable")
    else:
        path.write_bytes(encoded)


def _write_package(root: Path, raw_ledger: list[dict]) -> None:
    artifacts = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative not in PACKAGE_EXCLUSIONS:
            artifacts[relative] = _sha256(path)
    payload = {"schema": "vela.pn2d_minimal6_inverse_package_manifest.v1",
               "exclusions": sorted(PACKAGE_EXCLUSIONS),
               "raw_inputs": raw_ledger, "artifacts": artifacts}
    _write_or_validate(root / "package_manifest.json", payload, "package manifest")


def verify_report(root: str | Path) -> dict:
    report_root = Path(root).resolve()
    if not report_root.is_dir():
        raise ValueError("report root is not a directory")
    # Validate an existing package envelope before trusting any persisted report.
    _verify_package_manifest(report_root)
    report = _load_json(report_root / "physics_inverse_audit.json")
    if not isinstance(report, dict):
        raise ValueError("authoritative report must be an object")
    bundle, roots, raw_ledger = _load_raw_bundle(report)
    _verify_package_manifest(report_root, raw_ledger)
    report_manifest = _verify_report_manifest(report_root, raw_ledger)
    rows = canonical_observations(bundle)
    checks = _verify_raw_semantics(report_root, report, bundle, roots, rows)
    checks.extend(_verify_figures(report_root))
    checks.extend(("raw_input_hashes", "report_artifact_hashes"))
    result = {
        "schema": "vela.pn2d_minimal6_inverse_verification.v1", "passed": True,
        "checks": sorted(checks),
        "report_manifest_sha256": _sha256(report_root / "report_manifest.json"),
        "verified_input_count": len(raw_ledger),
        "verified_artifact_count": len(report_manifest["artifacts"]),
    }
    _write_or_validate(report_root / "verification.json", result, "verification result")
    _write_package(report_root, raw_ledger)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    result = verify_report(parse_args().report_root)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
