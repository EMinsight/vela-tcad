#!/usr/bin/env python3
"""Render three source-bound PN2D Minimal6 figures for a daily report."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import tempfile
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, UnidentifiedImageError

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import scripts.compare_pn2d_minimal6_diagnostic_sweeps as comparison


DAILY_FIGURE_NAMES = (
    "terminal_current.png",
    "maximum_field.png",
    "source_integrals.png",
)
SCHEMA = "vela.pn2d_minimal6_daily_report_figures.v1"
AXIS_ORDER = "decreasing_left_to_right"
MANIFEST_NAME = "daily_report_figure_manifest.json"
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
VOLTAGE_TICKS = (-1.0, -5.0, -10.0, -15.0, -20.0)


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _reject_nonfinite(value: str) -> None:
    raise ValueError(f"non-finite JSON constant {value}")


def _load_verified_report(
    report_path: Path,
) -> tuple[dict[str, Any], dict[str, str]]:
    report_path = report_path.resolve()
    comparison.verify_comparison_artifacts(report_path)
    payload = report_path.read_bytes()
    report = json.loads(payload.decode("utf-8"), parse_constant=_reject_nonfinite)
    return report, {"path": str(report_path), "sha256": _sha(payload)}


def _common_biases(report: dict[str, Any]) -> list[float]:
    biases = sorted(
        {float(row["bias_V"]) for row in report["checkpoints"]}, reverse=True
    )
    if not biases:
        raise ValueError("daily-report figures require exact common checkpoints")
    if any(not math.isfinite(value) for value in biases):
        raise ValueError("daily-report voltage axis contains non-finite values")
    limits = [biases[0], biases[-1]]
    if limits != [-1.0, -20.0] or not limits[0] > limits[1]:
        raise ValueError("daily-report voltage axis must run from -1 V to -20 V")
    return biases


def _rows(
    report: dict[str, Any], solver: str, topology: str
) -> list[dict[str, Any]]:
    selected = [
        {"bias_V": float(row["bias_V"]), "state": row[solver]}
        for row in report["checkpoints"]
        if row["topology"] == topology
    ]
    return sorted(selected, key=lambda row: row["bias_V"], reverse=True)


def _finite_values(
    rows: list[dict[str, Any]], observable: str, *, absolute: bool
) -> list[float]:
    values = [float(row["state"]["observables"][observable]) for row in rows]
    if any(not math.isfinite(value) for value in values):
        raise ValueError(f"daily-report observable is non-finite: {observable}")
    return [abs(value) for value in values] if absolute else values


def _finish_descending_axis(
    fig: plt.Figure,
    ax: plt.Axes,
    title: str,
    ylabel: str,
    limits: list[float],
) -> None:
    comparison._finish(fig, ax, title, ylabel)
    ax.set_xlim(*limits)
    ax.set_xticks(VOLTAGE_TICKS)


def _save(
    fig: plt.Figure,
    path: Path,
    source: dict[str, str],
    limits: list[float],
    series: list[dict[str, Any]],
    markers: list[dict[str, Any]],
) -> dict[str, Any]:
    metadata = {
        **comparison.FIGURE_METADATA,
        "VoltageAxisOrder": AXIS_ORDER,
        "VoltageAxisLimitsV": json.dumps(limits, separators=(",", ":")),
        "SeriesIdentities": json.dumps(
            series, sort_keys=True, separators=(",", ":")
        ),
        "FailureTransitionMarkers": json.dumps(
            markers, sort_keys=True, separators=(",", ":")
        ),
    }
    fig.savefig(path, dpi=comparison.FIGURE_DPI, metadata=metadata)
    plt.close(fig)
    payload = path.read_bytes()
    return {
        "source_comparison_path": source["path"],
        "source_comparison_sha256": source["sha256"],
        "sha256": _sha(payload),
        "width_px": comparison.FIGURE_WIDTH_PX,
        "height_px": comparison.FIGURE_HEIGHT_PX,
        "x_quantity": "applied_bias_V",
        "x_axis_order": AXIS_ORDER,
        "x_limits_V": limits,
        "series_identities": series,
    }


def _render(
    report: dict[str, Any], source: dict[str, str], out_dir: Path
) -> dict[str, dict[str, Any]]:
    biases = _common_biases(report)
    limits = [biases[0], biases[-1]]
    markers = comparison._failure_marker_identities(report)
    entries: dict[str, dict[str, Any]] = {}

    for stem, observable, title, unit, absolute, quantity in (
        (
            "terminal_current",
            "anode_current_A_per_um",
            "Terminal current at accepted exact checkpoints",
            "A/um",
            True,
            "terminal_current",
        ),
        (
            "maximum_field",
            "max_field_V_per_m",
            "Maximum electric field at accepted exact checkpoints",
            "V/m",
            False,
            "maximum_field",
        ),
    ):
        fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=comparison.FIGURE_DPI)
        series: list[dict[str, Any]] = []
        for solver in ("vela", "sentaurus"):
            for topology in ("sketch", "mirror"):
                selected = _rows(report, solver, topology)
                values = _finite_values(selected, observable, absolute=absolute)
                ax.plot(
                    [row["bias_V"] for row in selected],
                    values,
                    marker="o",
                    label=f"{solver} {topology}",
                )
                series.append(
                    {"solver": solver, "topology": topology, "quantity": quantity}
                )
        comparison._mark_failure_transitions(ax, report)
        ax.legend(loc="upper left")
        _finish_descending_axis(fig, ax, title, unit, limits)
        name = f"{stem}.png"
        entries[name] = _save(
            fig, out_dir / name, source, limits, series, markers
        )

    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=comparison.FIGURE_DPI)
    source_series: list[dict[str, Any]] = []
    for solver in ("vela", "sentaurus"):
        for topology in ("sketch", "mirror"):
            selected = _rows(report, solver, topology)
            for observable, quantity, label in (
                (
                    "native_source_integral_s_inv_per_cm",
                    "native_source",
                    "native source",
                ),
                (
                    "reconstructed_source_integral_s_inv_per_cm",
                    "reconstructed_source",
                    "reconstructed source",
                ),
            ):
                values = _finite_values(selected, observable, absolute=True)
                ax.plot(
                    [row["bias_V"] for row in selected],
                    values,
                    marker="o",
                    label=f"{solver} {topology} {label}",
                )
                source_series.append(
                    {"solver": solver, "topology": topology, "quantity": quantity}
                )
    comparison._mark_failure_transitions(ax, report)
    ax.legend(loc="upper left")
    _finish_descending_axis(
        fig,
        ax,
        "Native and reconstructed avalanche sources",
        "s^-1 per 1 cm depth",
        limits,
    )
    entries["source_integrals.png"] = _save(
        fig,
        out_dir / "source_integrals.png",
        source,
        limits,
        source_series,
        markers,
    )
    return entries


def render_daily_report_figures(
    report_path: Path, out_dir: Path
) -> dict[str, Any]:
    report, source = _load_verified_report(Path(report_path))
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    entries = _render(report, source, out_dir)
    contract = {
        "schema": SCHEMA,
        "source_comparison": source,
        "figures": {name: entries[name] for name in DAILY_FIGURE_NAMES},
    }
    (out_dir / MANIFEST_NAME).write_text(
        json.dumps(contract, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return contract


def verify_daily_report_figures(manifest_path: Path) -> bool:
    manifest_path = Path(manifest_path).resolve()
    contract = json.loads(
        manifest_path.read_text(encoding="utf-8"), parse_constant=_reject_nonfinite
    )
    if set(contract) != {"schema", "source_comparison", "figures"}:
        raise ValueError("daily-report manifest has invalid top-level keys")
    if (
        contract["schema"] != SCHEMA
        or set(contract["figures"]) != set(DAILY_FIGURE_NAMES)
    ):
        raise ValueError("daily-report figure contract is invalid")

    source = contract["source_comparison"]
    if (
        set(source) != {"path", "sha256"}
        or not isinstance(source["path"], str)
        or not SHA256.fullmatch(source["sha256"])
    ):
        raise ValueError("daily-report source binding is invalid")
    report, actual_source = _load_verified_report(Path(source["path"]))
    if actual_source != source:
        raise ValueError("daily-report source comparison hash mismatch")

    expected_keys = {
        "source_comparison_path",
        "source_comparison_sha256",
        "sha256",
        "width_px",
        "height_px",
        "x_quantity",
        "x_axis_order",
        "x_limits_V",
        "series_identities",
    }
    for name in DAILY_FIGURE_NAMES:
        entry = contract["figures"][name]
        if set(entry) != expected_keys:
            raise ValueError(f"daily-report figure entry is invalid: {name}")
        if (
            entry["source_comparison_path"] != source["path"]
            or entry["source_comparison_sha256"] != source["sha256"]
        ):
            raise ValueError(f"daily-report figure source binding mismatch: {name}")
        limits = entry["x_limits_V"]
        if (
            entry["x_quantity"] != "applied_bias_V"
            or entry["x_axis_order"] != AXIS_ORDER
            or limits != [-1.0, -20.0]
            or not limits[0] > limits[1]
        ):
            raise ValueError(f"daily-report voltage axis contract mismatch: {name}")
        path = manifest_path.parent / name
        if (
            not path.is_file()
            or not SHA256.fullmatch(entry["sha256"])
            or _sha(path.read_bytes()) != entry["sha256"]
        ):
            raise ValueError(f"daily-report PNG hash mismatch: {name}")
        try:
            with Image.open(path) as image:
                if image.size != (
                    comparison.FIGURE_WIDTH_PX,
                    comparison.FIGURE_HEIGHT_PX,
                ):
                    raise ValueError(
                        f"daily-report PNG dimensions mismatch: {name}"
                    )
                if (
                    image.info.get("VoltageAxisOrder") != AXIS_ORDER
                    or json.loads(image.info["VoltageAxisLimitsV"]) != limits
                    or json.loads(image.info["SeriesIdentities"])
                    != entry["series_identities"]
                ):
                    raise ValueError(f"daily-report PNG metadata mismatch: {name}")
        except (
            OSError,
            KeyError,
            json.JSONDecodeError,
            UnidentifiedImageError,
        ) as exc:
            raise ValueError(f"daily-report PNG decode mismatch: {name}") from exc

    with tempfile.TemporaryDirectory(
        prefix="pn2d-minimal6-daily-report-"
    ) as temp:
        rerender = _render(report, source, Path(temp))
        for name in DAILY_FIGURE_NAMES:
            with Image.open(manifest_path.parent / name) as actual, Image.open(
                Path(temp) / name
            ) as expected:
                if (
                    actual.convert("RGB").tobytes()
                    != expected.convert("RGB").tobytes()
                ):
                    raise ValueError(
                        f"daily-report deterministic pixel mismatch: {name}"
                    )
            if (
                rerender[name]["series_identities"]
                != contract["figures"][name]["series_identities"]
            ):
                raise ValueError(f"daily-report rerender identity mismatch: {name}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-report", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    render_daily_report_figures(args.comparison_report, args.out_dir)
    manifest = args.out_dir.resolve() / MANIFEST_NAME
    verify_daily_report_figures(manifest)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
