#!/usr/bin/env python3
"""Independently verify the Minimal6 element-avalanche replay evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median


Q_LEGACY_C = 1.6021918e-19


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    return parser.parse_args()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def close(a: float, b: float, relative: float = 1.0e-12) -> bool:
    return abs(a - b) <= relative * max(abs(a), abs(b), 1.0e-300)


def f(value: str) -> float:
    return float(value)


def main() -> int:
    root = args.root.resolve()
    analysis = root / "analysis"
    controls = root / "controls"
    failures: list[str] = []

    manifest = json.loads(
        (analysis / "manifest_closed.json").read_text(encoding="ascii")
    )
    for name, digest in manifest["output_sha256"].items():
        path = analysis / name
        if not path.is_file() or sha256(path) != digest:
            failures.append(f"analysis hash mismatch: {name}")
    control_manifest = json.loads(
        (controls / "manifest.json").read_text(encoding="ascii")
    )
    for name, digest in control_manifest["output_sha256"].items():
        path = controls / name
        if not path.is_file() or sha256(path) != digest:
            failures.append(f"control hash mismatch: {name}")

    reconstructions = rows(analysis / "element_reconstructions.csv")
    if len(reconstructions) != 2 * 3 * 4 * 2 * 5:
        failures.append("unexpected element-reconstruction lattice size")
    by_key = {
        (
            row["topology"],
            row["bias_V"],
            row["element"],
            row["carrier"],
            row["candidate"],
        ): row
        for row in reconstructions
    }
    for row in reconstructions:
        expected = (
            f(row["alpha_cm_inv"])
            * f(row["magnitude_A_cm2"])
            / Q_LEGACY_C
        )
        if not close(expected, f(row["generation_cm3_s"]), 2.0e-15):
            failures.append("element alpha-current generation mismatch")
            break

    base_keys = {
        key[:-1]
        for key in by_key
        if key[-1] == "gss_laux_edge_volume_weighted"
    }
    for key in base_keys:
        laux = by_key[(*key, "gss_laux_edge_volume_weighted")]
        active = by_key[(*key, "box_active_edge_exact")]
        for field in (
            "vector_x_A_cm2",
            "vector_y_A_cm2",
            "magnitude_A_cm2",
            "generation_cm3_s",
        ):
            if not close(f(laux[field]), f(active[field]), 2.0e-14):
                failures.append(f"Laux/active mismatch: {key} {field}")

    state_rows = rows(analysis / "state_source_summary_corrected.csv")
    laux_errors = [
        f(row["integral_absolute_error_dex"])
        for row in state_rows
        if row["candidate"] == "gss_laux_edge_volume_weighted"
    ]
    native_errors = [
        f(row["integral_absolute_error_dex"])
        for row in state_rows
        if row["candidate"] == "native_element_vector_control"
    ]
    vela_m20_errors = [
        f(row["integral_absolute_error_dex"])
        for row in state_rows
        if row["candidate"] == "vela_triangle_proxy_existing"
        and f(row["bias_V"]) == -20.0
    ]
    if len(laux_errors) != 12 or max(laux_errors) > 5.0e-4:
        failures.append("GSS/Laux integral closure exceeds 5e-4 dex")
    if median(native_errors) < 6.0:
        failures.append("native element-vector rejection is not reproduced")
    if len(vela_m20_errors) != 2 or min(vela_m20_errors) < 11.9:
        failures.append("corrected Vela -20 V source gap is not reproduced")

    box_rows = rows(analysis / "box_current_closure.csv")
    max_contact = max(
        f(row["relative_error"])
        for row in box_rows
        if row["location"] in {"Anode", "Cathode"}
    )
    max_kcl = max(
        f(row["relative_error"])
        for row in box_rows
        if row["location"].startswith("internal_vertex_")
    )
    if max_contact > 2.4e-3:
        failures.append("carrier contact closure exceeds 2.4e-3")
    if max_kcl > 3.2e-9:
        failures.append("internal total-current KCL exceeds 3.2e-9")

    source_rows = rows(analysis / "source_integral_closure.csv")
    max_source = max(f(row["relative_error"]) for row in source_rows)
    if max_source > 2.0e-15:
        failures.append("ReadMeasure/CurrentPlot closure exceeds 2e-15")

    comparison = rows(controls / "control_vs_default.csv")
    volume_rows = [
        row
        for row in comparison
        if row["variant"] == "element_volume_avalanche"
    ]
    for row in volume_rows:
        for field, value in row.items():
            if field in {"topology", "bias_V", "variant"}:
                continue
            if f(value) != 0.0:
                failures.append("ElementVolumeAvalanche is not exact-zero")
                break
    gradqf_m20 = [
        row
        for row in comparison
        if row["variant"] == "aval_dens_grad_qf"
        and f(row["bias_V"]) == -20.0
    ]
    if len(gradqf_m20) != 2:
        failures.append("missing AvalDensGradQF -20 V controls")
    else:
        values = [f(row["total_qg_log10_ratio_to_default_dex"]) for row in gradqf_m20]
        if max(values) - min(values) > 1.0e-15:
            failures.append("mirror/sketch AvalDensGradQF mismatch")
        if not 0.262 <= values[0] <= 0.263:
            failures.append("unexpected AvalDensGradQF -20 V source shift")

    unsupported_log = (
        root
        / "raw"
        / "mirror"
        / "element_vertex_probe"
        / "run_element_vertex_probe.out"
    )
    unsupported_text = unsupported_log.read_text(
        encoding="ascii", errors="strict"
    )
    unsupported_message = (
        "Tried to read undefined ElementVertex-RegionWise "
        "eAvalancheGeneration"
    )
    if unsupported_message not in unsupported_text:
        failures.append("missing typed unsupported element-vertex evidence")

    result = {
        "status": "passed" if not failures else "failed",
        "failure_count": len(failures),
        "failures": failures,
        "checks": {
            "element_reconstruction_rows": len(reconstructions),
            "gss_laux_max_integral_error_dex": max(laux_errors),
            "native_element_vector_median_error_dex": median(native_errors),
            "vela_m20_min_error_dex": min(vela_m20_errors),
            "max_carrier_contact_relative_error": max_contact,
            "max_internal_kcl_relative_error": max_kcl,
            "max_readmeasure_currentplot_relative_error": max_source,
            "element_vertex_generation_support": "unsupported",
        },
    }
    output = root / "independent_verification.json"
    output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    args = parse_args()
    raise SystemExit(main())
