#!/usr/bin/env python3
"""Invert the Sentaurus Minimal6 high-field mobility driving force."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

from pn2d_minimal6_diagnostics.highfield_box_replay import (
    CELL_MAPPING,
)


FIELD = {
    "electron": {"saturation_velocity": 1.07e5, "beta": 1.109},
    "hole": {"saturation_velocity": 8.37e4, "beta": 1.213},
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def abs_log_error(candidate: float, reference: float) -> float:
    if candidate <= 0.0 or reference <= 0.0:
        raise ValueError("log error requires positive values")
    return abs(math.log10(candidate / reference))


def field_limited_mobility(
    carrier: str, low_field_mobility: float, field_V_per_m: float
) -> float:
    parameters = FIELD[carrier]
    ratio = (
        low_field_mobility
        * abs(field_V_per_m)
        / parameters["saturation_velocity"]
    )
    return low_field_mobility / (
        1.0 + ratio ** parameters["beta"]
    ) ** (1.0 / parameters["beta"])


def inverted_field(
    carrier: str, low_field_mobility: float, final_mobility: float
) -> float:
    if final_mobility > low_field_mobility:
        raise ValueError("final mobility exceeds low-field mobility")
    parameters = FIELD[carrier]
    beta = parameters["beta"]
    power = max(
        0.0, (low_field_mobility / final_mobility) ** beta - 1.0
    )
    return (
        parameters["saturation_velocity"]
        / low_field_mobility
        * power ** (1.0 / beta)
    )


def stats(values: list[float]) -> dict[str, float | int]:
    if not values:
        raise ValueError("cannot summarize empty samples")
    ordered = sorted(values)
    return {
        "count": len(values),
        "median": statistics.median(values),
        "p95": ordered[math.ceil(0.95 * len(ordered)) - 1],
        "maximum": ordered[-1],
    }


def load_low_field(
    topology: str, directory: Path
) -> dict[tuple[str, str, int], float]:
    result: dict[tuple[str, str, int], float] = {}
    for carrier, quantity in (
        ("electron", "eMobility"),
        ("hole", "hMobility"),
    ):
        path = directory / f"{quantity}_region0_cells.csv"
        rows = read_csv(path)
        if [int(row["cell_id"]) for row in rows] != [0, 1, 2, 3]:
            raise ValueError(f"{path} lacks canonical cells 0..3")
        by_region_cell = {
            int(row["cell_id"]): float(row["component0"]) * 1.0e-4
            for row in rows
        }
        for vela_cell, region_cell in enumerate(CELL_MAPPING[topology]):
            result[(topology, carrier, vela_cell)] = by_region_cell[
                region_cell
            ]
    return result


def validate_control_deck(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    required = (
        "DopingDependence",
        "eMobility/Element",
        "hMobility/Element",
    )
    if any(token not in text for token in required):
        raise ValueError(f"low-field deck lacks required contract: {path}")
    if "HighFieldSaturation" in text:
        raise ValueError(f"low-field deck still enables high-field saturation: {path}")


def run(args: argparse.Namespace) -> dict[str, object]:
    paths = {
        "mirror_low_field": args.mirror_low_field.resolve(),
        "sketch_low_field": args.sketch_low_field.resolve(),
        "mirror_deck": args.mirror_deck.resolve(),
        "sketch_deck": args.sketch_deck.resolve(),
        "native": args.native.resolve(),
        "triangle": args.triangle.resolve(),
        "transport": args.transport.resolve(),
    }
    for name in ("mirror_deck", "sketch_deck"):
        validate_control_deck(paths[name])

    low_field = {}
    low_field.update(load_low_field("mirror", paths["mirror_low_field"]))
    low_field.update(load_low_field("sketch", paths["sketch_low_field"]))
    if len(low_field) != 16:
        raise ValueError("expected two topologies x two carriers x four cells")

    native_rows = read_csv(paths["native"])
    triangle_rows = read_csv(paths["triangle"])
    transport_rows = read_csv(paths["transport"])
    triangle = {
        (
            row["topology"],
            float(row["bias_V"]),
            int(row["cell_id"]),
            row["carrier"],
        ): row
        for row in triangle_rows
    }
    transport = {
        (
            row["topology"],
            float(row["bias_V"]),
            int(row["cell_id"]),
        ): row
        for row in transport_rows
    }
    if len(native_rows) != 320 or len(triangle) != 320 or len(transport) != 160:
        raise ValueError("expected 40 states x four cells and two carriers")

    samples: list[dict[str, object]] = []
    for native in native_rows:
        topology = native["topology"]
        bias = float(native["bias_V"])
        cell = int(native["cell_id"])
        carrier = native["carrier"]
        key = (topology, bias, cell, carrier)
        transport_key = (topology, bias, cell)
        low = low_field[(topology, carrier, cell)]
        final = float(native["sentaurus_native_final_m2_per_Vs"])
        inferred = inverted_field(carrier, low, final)
        native_qf = float(native["sentaurus_native_qf_field_V_per_m"])
        triangle_qf = float(triangle[key]["cell_qf_field_V_per_m"])
        transport_row = transport[transport_key]
        electric = math.hypot(
            float(transport_row["electric_field_x_V_per_m"]),
            float(transport_row["electric_field_y_V_per_m"]),
        )
        fields = {
            "native_qf": native_qf,
            "triangle_qf": triangle_qf,
            "electric": electric,
        }
        replay = {
            name: field_limited_mobility(carrier, low, field)
            for name, field in fields.items()
        }
        samples.append(
            {
                "topology": topology,
                "bias_V": bias,
                "cell_id": cell,
                "carrier": carrier,
                "sentaurus_low_field_m2_per_Vs": low,
                "vela_cell_average_low_field_m2_per_Vs": float(
                    native["vela_cell_average_low_field_m2_per_Vs"]
                ),
                "low_field_abs_error_dex": abs_log_error(
                    low,
                    float(native["vela_cell_average_low_field_m2_per_Vs"]),
                ),
                "sentaurus_final_m2_per_Vs": final,
                "inverted_effective_field_V_per_m": inferred,
                "native_qf_field_V_per_m": native_qf,
                "triangle_qf_field_V_per_m": triangle_qf,
                "electric_field_V_per_m": electric,
                "inverted_vs_native_qf_abs_error_dex": abs_log_error(
                    inferred, native_qf
                ),
                "inverted_vs_triangle_qf_abs_error_dex": abs_log_error(
                    inferred, triangle_qf
                ),
                "inverted_vs_electric_abs_error_dex": abs_log_error(
                    inferred, electric
                ),
                "native_qf_replay_m2_per_Vs": replay["native_qf"],
                "triangle_qf_replay_m2_per_Vs": replay["triangle_qf"],
                "electric_replay_m2_per_Vs": replay["electric"],
                "native_qf_replay_abs_error_dex": abs_log_error(
                    replay["native_qf"], final
                ),
                "triangle_qf_replay_abs_error_dex": abs_log_error(
                    replay["triangle_qf"], final
                ),
                "electric_replay_abs_error_dex": abs_log_error(
                    replay["electric"], final
                ),
            }
        )

    summary: list[dict[str, object]] = []
    for carrier in ("electron", "hole"):
        carrier_rows = [row for row in samples if row["carrier"] == carrier]
        row: dict[str, object] = {
            "carrier": carrier,
            "sample_count": len(carrier_rows),
        }
        for output_name, field_name in (
            ("low_field_vs_vela", "low_field_abs_error_dex"),
            ("inverted_vs_native_qf", "inverted_vs_native_qf_abs_error_dex"),
            (
                "inverted_vs_triangle_qf",
                "inverted_vs_triangle_qf_abs_error_dex",
            ),
            ("inverted_vs_electric", "inverted_vs_electric_abs_error_dex"),
            ("native_qf_replay", "native_qf_replay_abs_error_dex"),
            ("triangle_qf_replay", "triangle_qf_replay_abs_error_dex"),
            ("electric_replay", "electric_replay_abs_error_dex"),
        ):
            result = stats([float(item[field_name]) for item in carrier_rows])
            for statistic in ("median", "p95", "maximum"):
                row[f"{output_name}_{statistic}_dex"] = result[statistic]
        summary.append(row)

    report_lines = [
        "# PN2D Minimal6 Sentaurus high-field driving-force control",
        "",
        "Status: `valid`",
        "",
        "Primary outcome: `electric_field_best_supported_candidate`",
        "",
        "The control disables only `HighFieldSaturation` in one -20 V run per "
        "topology and exports native element mobility. Because the retained "
        "`DopingDependence` Masetti branch depends on doping and temperature, "
        "the resulting element mobility is used as the Sentaurus-native "
        "low-field coefficient for all 40 exact states.",
        "",
        "## Full 40-state replay",
        "",
        "| Carrier | Candidate drive | N | Median mobility error (dex) | "
        "P95 (dex) | Maximum (dex) |",
        "|---|---|---:|---:|---:|---:|",
    ]
    labels = (
        ("native_qf", "exported native QFP gradient"),
        ("triangle_qf", "affine triangle QFP gradient"),
        ("electric", "native element electric field"),
    )
    for row in summary:
        for prefix, label in labels:
            report_lines.append(
                f"| {row['carrier']} | {label} | {row['sample_count']} | "
                f"{float(row[f'{prefix}_replay_median_dex']):.6g} | "
                f"{float(row[f'{prefix}_replay_p95_dex']):.6g} | "
                f"{float(row[f'{prefix}_replay_maximum_dex']):.6g} |"
            )
    report_lines.extend(
        [
            "",
            "The native element electric field is the strongest tested "
            "driving-force candidate. This result rejects treating exported "
            "electron `eGradQuasiFermi/Element` as the internal mobility drive.",
            "",
            "## Boundary of inference",
            "",
            "This is a controlled Sentaurus-operator reconstruction, not a "
            "direct observation of a proprietary internal flag. The low-field "
            "control changes the solved state, but its Masetti coefficient is "
            "state-independent at fixed doping and 300 K. No production Vela "
            "formula is modified by this diagnostic.",
            "",
        ]
    )

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "highfield_drive_samples.csv", samples)
    write_csv(output / "highfield_drive_summary.csv", summary)
    (output / "report.md").write_text(
        "\n".join(report_lines), encoding="utf-8"
    )
    outputs = {}
    for name in (
        "highfield_drive_samples.csv",
        "highfield_drive_summary.csv",
        "report.md",
    ):
        outputs[name] = {"sha256": sha256(output / name)}
    manifest = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "pn2d_minimal6_sentaurus_highfield_drive_control",
        "remote_sentaurus_release": args.sentaurus_release,
        "state_count": 40,
        "carrier_element_sample_count": len(samples),
        "control_runs": 2,
        "control_contract": {
            "retained": [
                "DopingDependence",
                "SRH",
                "Avalanche(VanOverstraeten)",
                "EffectiveIntrinsicDensity(OldSlotboom)",
            ],
            "disabled": ["HighFieldSaturation"],
            "temperature_K": 300.0,
            "native_element_mobility": True,
        },
        "primary_outcome": "electric_field_best_supported_candidate",
        "production_formula_modified": False,
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in paths.items()
            if path.is_file()
        },
        "input_directories": {
            name: str(path)
            for name, path in paths.items()
            if path.is_dir()
        },
        "low_field_members": {
            f"{topology}/{quantity}_region0_cells.csv": {
                "path": str(
                    paths[f"{topology}_low_field"]
                    / f"{quantity}_region0_cells.csv"
                ),
                "sha256": sha256(
                    paths[f"{topology}_low_field"]
                    / f"{quantity}_region0_cells.csv"
                ),
            }
            for topology in ("mirror", "sketch")
            for quantity in ("eMobility", "hMobility")
        },
        "outputs": outputs,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mirror-low-field", type=Path, required=True)
    parser.add_argument("--sketch-low-field", type=Path, required=True)
    parser.add_argument("--mirror-deck", type=Path, required=True)
    parser.add_argument("--sketch-deck", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--triangle", type=Path, required=True)
    parser.add_argument("--transport", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--sentaurus-release", default="O-2018.06-SP2"
    )
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
