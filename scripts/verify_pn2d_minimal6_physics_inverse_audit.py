#!/usr/bin/env python3
"""Independent semantic and integrity verifier for an inverse-audit package.

This module deliberately does not import report-building or plotting functions.
Selected physics identities are recomputed directly from persisted operands.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from PIL import Image


REPORT_EXCLUSIONS = {"report_manifest.json", "verification.json", "package_manifest.json"}
PACKAGE_EXCLUSIONS = {"package_manifest.json"}
DISCOVERY_KEYS = (("sketch", -1.0), ("sketch", -4.0), ("sketch", -8.0),
                  ("sketch", -12.0), ("sketch", -16.0), ("sketch", -19.0),
                  ("sketch", -20.0))
COMMON_KEYS = tuple(sorted((topology, float(-bias)) for topology in ("sketch", "mirror")
                           for bias in range(1, 21)))
THRESHOLDS = {
    "field_median_relative": 0.02, "field_median_angle_deg": 1.0,
    "gradient_median_abs_dex": 0.1, "gradient_p95_abs_dex": 0.3,
    "gradient_median_angle_deg": 5.0, "integrated_generation_abs_dex": 0.1,
    "local_generation_abs_dex": 0.3, "replacement_closure_abs_dex": 1.0e-10,
}
FIGURES = ("potential_field", "qf_gradient", "current_density",
           "alpha_generation", "replacement_matrix")
Q = 1.602176634e-19


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pixel_sha256(path: Path) -> str:
    with Image.open(path) as image:
        return hashlib.sha256(image.convert("RGBA").tobytes()).hexdigest()


def _load_json(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"), parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid JSON: {path.name}") from error


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def _same(a: float, b: float, *, tolerance: float = 1.0e-12) -> bool:
    return math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def _verify_report_manifest(root: Path) -> dict:
    manifest = _load_json(root / "report_manifest.json")
    if not isinstance(manifest, dict) or manifest.get("schema") != "vela.pn2d_minimal6_inverse_report_manifest.v1":
        raise ValueError("report manifest schema mismatch")
    if set(manifest.get("exclusions", ())) != REPORT_EXCLUSIONS:
        raise ValueError("report manifest exclusions mismatch")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        raise ValueError("report manifest artifact ledger is empty")
    actual = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative not in REPORT_EXCLUSIONS:
            actual[relative] = _sha256(path)
    if actual != artifacts:
        raise ValueError("report artifact hash mismatch")
    inputs = manifest.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError("raw input ledger is empty")
    for item in inputs:
        path = Path(item["path"])
        if not path.is_file() or _sha256(path) != item.get("sha256"):
            raise ValueError(f"raw input hash mismatch: {item.get('logical_id')}")
        if path.name != Path(item["relative_path"]).name:
            raise ValueError("raw input path binding mismatch")
    return manifest


def _verify_package_manifest(root: Path) -> None:
    path = root / "package_manifest.json"
    if not path.exists():
        return
    manifest = _load_json(path)
    if not isinstance(manifest, dict) or manifest.get("schema") != "vela.pn2d_minimal6_inverse_package_manifest.v1":
        raise ValueError("package manifest schema mismatch")
    if set(manifest.get("exclusions", ())) != PACKAGE_EXCLUSIONS:
        raise ValueError("package manifest exclusions mismatch")
    actual = {}
    for member in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = member.relative_to(root).as_posix()
        if relative not in PACKAGE_EXCLUSIONS:
            actual[relative] = _sha256(member)
    if actual != manifest.get("artifacts"):
        raise ValueError("package manifest hash mismatch")


def _triangle_gradient(points: list[list[float]], values: list[float]) -> tuple[float, float]:
    (x0, y0), (x1, y1), (x2, y2) = points
    f0, f1, f2 = values
    determinant = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if not math.isfinite(determinant) or abs(determinant) <= 1.0e-300:
        raise ValueError("semantic triangle is degenerate")
    return (((f1 - f0) * (y2 - y0) - (f2 - f0) * (y1 - y0)) / determinant,
            ((x1 - x0) * (f2 - f0) - (x2 - x0) * (f1 - f0)) / determinant)


def _verify_semantics(report: dict) -> list[str]:
    if set(report) != {"schema", "diagnostic_only", "phase_base", "payload"}:
        raise ValueError("authoritative report top-level contract mismatch")
    if report["schema"] != "vela.pn2d_minimal6_physics_inverse_audit.v1" or report["diagnostic_only"] is not True:
        raise ValueError("authoritative report schema mismatch")
    if report["phase_base"] != "a5524cf":
        raise ValueError("production phase baseline mismatch")
    payload = report["payload"]
    discovery = tuple((str(item[0]), float(item[1])) for item in payload["discovery_keys"])
    holdout = tuple((str(item[0]), float(item[1])) for item in payload["holdout_keys"])
    if discovery != DISCOVERY_KEYS or holdout != tuple(key for key in COMMON_KEYS if key not in DISCOVERY_KEYS):
        raise ValueError("discovery/holdout membership mismatch")
    if payload["thresholds"] != THRESHOLDS:
        raise ValueError("acceptance threshold mismatch")
    if payload.get("production_cpp_changed") is not False:
        raise ValueError("production baseline guard failed")
    replay = payload["localization_control"]["semantic_replay"]
    checks = ["split_membership", "fixed_thresholds", "production_baseline"]
    triangle = replay["triangle_gradient"]
    if triangle["status"] == "valid":
        expected = _triangle_gradient(triangle["points_m"], triangle["values_V"])
        if not all(_same(a, b) for a, b in zip(expected, triangle["value_V_per_m"])):
            raise ValueError("triangle-gradient semantic replay mismatch")
        checks.append("triangle_gradient")
    gradient = replay["current_inverted_gradient"]
    if gradient["status"] == "valid":
        scale = Q * float(gradient["mobility_m2_per_Vs"]) * float(gradient["density_m3"])
        if not math.isfinite(scale) or scale <= 0.0:
            raise ValueError("current-inverted gradient operands invalid")
        expected = tuple(-float(value) / scale for value in gradient["current_A_per_m2"])
        if not all(_same(a, b) for a, b in zip(expected, gradient["value_V_per_m"])):
            raise ValueError("current-inverted gradient semantic replay mismatch")
        checks.append("current_inverted_gradient")
    alpha = replay["inverse_alpha"]
    if alpha["status"] == "valid":
        coefficient = float(alpha["alpha_m_inv"])
        ceiling = float(alpha["gamma"]) * float(alpha["prefactor_m_inv"])
        expected = -float(alpha["gamma"]) * float(alpha["critical_field_V_per_m"]) / math.log(coefficient / ceiling)
        if not _same(expected, float(alpha["field_V_per_m"])):
            raise ValueError("inverse-alpha semantic replay mismatch")
        checks.append("inverse_alpha")
    generation = replay["generation"]
    if generation["status"] == "valid":
        expected = (float(generation["alpha_n_m_inv"]) * math.hypot(*generation["jn_A_per_m2"])
                    + float(generation["alpha_p_m_inv"]) * math.hypot(*generation["jp_A_per_m2"])) / Q
        if not _same(expected, float(generation["value_m3_s_inv"]), tolerance=1.0e-11):
            raise ValueError("generation semantic replay mismatch")
        checks.append("generation_reconstruction")
    closure = replay["replacement_closure"]
    for name in ("forward_abs_dex", "reverse_abs_dex", "direct_abs_dex"):
        if not math.isfinite(float(closure[name])) or float(closure[name]) > 1.0e-10:
            raise ValueError("replacement closure semantic replay mismatch")
    checks.append("replacement_closure")
    return checks


def _verify_tables(root: Path, report: dict) -> list[str]:
    with (root / "candidate_metrics.csv").open(newline="", encoding="utf-8") as handle:
        metrics = list(csv.DictReader(handle))
    if len(metrics) != len(report["payload"]["candidate_metrics"]):
        raise ValueError("candidate metric row count mismatch")
    with (root / "candidate_classifications.json").open(encoding="utf-8") as handle:
        classifications = json.load(handle)["classifications"]
    if classifications != report["payload"]["classifications"]:
        raise ValueError("candidate classification semantic mismatch")
    allowed = {"identified", "consistent_nonunique", "confounded", "insufficient_data", "rejected"}
    if any(item["classification"] not in allowed for item in classifications):
        raise ValueError("candidate classification is not typed")
    for name in ("node", "edge", "cell", "contact", "integrated"):
        with (root / f"observations_{name}.csv").open(newline="", encoding="utf-8") as handle:
            tuple(csv.DictReader(handle))
    return ["canonical_csv_parse", "classification_consistency"]


def _verify_figures(root: Path) -> list[str]:
    manifest = _load_json(root / "figure_manifest.json")
    if not isinstance(manifest, dict):
        raise ValueError("figure manifest schema mismatch")
    entries = manifest.get("figures")
    if [item.get("name") for item in entries] != list(FIGURES):
        raise ValueError("fixed figure set mismatch")
    required_contract = {"question", "takeaway", "family", "variant", "row_grain_sufficiency",
                         "fields", "palette_policy", "output_paths", "qa_surface"}
    for item in entries:
        name = item["name"]
        png, pdf = root / "figures" / f"{name}.png", root / "figures" / f"{name}.pdf"
        if _sha256(png) != item["png_sha256"] or _sha256(pdf) != item["pdf_sha256"]:
            raise ValueError("figure file hash mismatch")
        if _pixel_sha256(png) != item["png_pixel_sha256"]:
            raise ValueError("figure decoded pixel hash mismatch")
        if set(item["chart_contract"]) != required_contract:
            raise ValueError("figure chart contract mismatch")
    return ["five_figure_pairs", "png_pixel_hashes", "chart_contracts"]


def _write_or_validate_verification(root: Path, result: dict) -> None:
    path = root / "verification.json"
    encoded = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")
    if path.exists():
        if path.read_bytes() != encoded:
            raise ValueError("verification result is not byte-stable")
    else:
        path.write_bytes(encoded)


def _write_or_validate_package(root: Path) -> None:
    artifacts = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative not in PACKAGE_EXCLUSIONS:
            artifacts[relative] = _sha256(path)
    payload = {
        "schema": "vela.pn2d_minimal6_inverse_package_manifest.v1",
        "exclusions": sorted(PACKAGE_EXCLUSIONS), "artifacts": artifacts,
    }
    path = root / "package_manifest.json"
    encoded = (json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")
    if path.exists():
        if path.read_bytes() != encoded:
            raise ValueError("package manifest is not byte-stable")
    else:
        path.write_bytes(encoded)


def verify_report(root: str | Path) -> dict:
    report_root = Path(root).resolve()
    if not report_root.is_dir():
        raise ValueError("report root is not a directory")
    _verify_package_manifest(report_root)
    manifest = _verify_report_manifest(report_root)
    report = _load_json(report_root / "physics_inverse_audit.json")
    if not isinstance(report, dict):
        raise ValueError("authoritative report must be an object")
    if _sha256(report_root / "input_manifest.json") != report["payload"]["input_manifest_sha256"]:
        raise ValueError("input manifest binding mismatch")
    checks = _verify_semantics(report)
    checks.extend(_verify_tables(report_root, report))
    checks.extend(_verify_figures(report_root))
    checks.extend(("raw_input_hashes", "report_artifact_hashes"))
    result = {
        "schema": "vela.pn2d_minimal6_inverse_verification.v1",
        "passed": True, "checks": sorted(checks),
        "report_manifest_sha256": _sha256(report_root / "report_manifest.json"),
        "verified_input_count": len(manifest["inputs"]),
        "verified_artifact_count": len(manifest["artifacts"]),
    }
    _write_or_validate_verification(report_root, result)
    _write_or_validate_package(report_root)
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
