#!/usr/bin/env python3
"""Factor the BVmethods NMOS avalanche source on identical Vela edges."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


Q_C = 1.602176634e-19
REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_VELA = (
    RUN_ROOT
    / "vela_validation/btbt_e2_iic_reclosure_20260804/postprocess_only/"
      "sg_avalanche_edges.csv"
)
DEFAULT_SENT = RUN_ROOT / "sentaurus_iic_multibias_exact_extended_20260803/imported"
DEFAULT_OUT = (
    RUN_ROOT / "vela_validation/btbt_e2_iic_edge_ledger_20260804"
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def sent_tag(bias: float) -> str:
    return f"iic_v{bias:.6f}".replace(".", "p")


def scalar_field(root: Path, name: str) -> dict[int, float]:
    result: dict[int, float] = {}
    for row in rows(root / "fields" / f"{name}_region3.csv"):
        result[int(row["node_id"])] = float(row["component0"])
    return result


def vector_field(root: Path, name: str) -> dict[int, tuple[float, float]]:
    result: dict[int, tuple[float, float]] = {}
    for row in rows(root / "fields" / f"{name}_region3.csv"):
        result[int(row["node_id"])] = (
            float(row["component0"]), float(row["component1"])
        )
    return result


def edge_scalar(row: dict[str, str], field: dict[int, float], scale: float) -> float:
    return 0.5 * (field[int(row["node0"])] + field[int(row["node1"])]) * scale


def edge_vector_projection(
    row: dict[str, str], field: dict[int, tuple[float, float]], scale: float
) -> float:
    node0, node1 = int(row["node0"]), int(row["node1"])
    dx = float(row["x1_um"]) - float(row["x0_um"])
    dy = float(row["y1_um"]) - float(row["y0_um"])
    length = math.hypot(dx, dy)
    if length <= 0.0:
        return 0.0
    vx = 0.5 * (field[node0][0] + field[node1][0])
    vy = 0.5 * (field[node0][1] + field[node1][1])
    return abs((vx * dx + vy * dy) / length) * scale


def edge_vector_magnitude(
    row: dict[str, str], field: dict[int, tuple[float, float]], scale: float
) -> float:
    node0, node1 = int(row["node0"]), int(row["node1"])
    vx = 0.5 * (field[node0][0] + field[node1][0])
    vy = 0.5 * (field[node0][1] + field[node1][1])
    return math.hypot(vx, vy) * scale


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-edges", type=Path, default=DEFAULT_VELA)
    parser.add_argument("--sentaurus-root", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--bias", type=float, default=6.4)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    sent = args.sentaurus_root / sent_tag(args.bias)
    sent_alpha_n = scalar_field(sent, "eAlphaAvalanche")
    sent_alpha_p = scalar_field(sent, "hAlphaAvalanche")
    sent_generation = scalar_field(sent, "ImpactIonization")
    sent_current_n = vector_field(sent, "eCurrentDensity")
    sent_current_p = vector_field(sent, "hCurrentDensity")
    supported = set(sent_generation)
    selected = [
        row for row in rows(args.vela_edges)
        if math.isclose(float(row["bias_V"]), args.bias, abs_tol=1.0e-10)
        and int(row["node0"]) in supported
        and int(row["node1"]) in supported
        and float(row["edge_area_proxy_m2"]) > 0.0
    ]
    if not selected:
        raise RuntimeError(f"no matched Vela edges at {args.bias:g} V")

    details: list[dict[str, Any]] = []
    totals = {
        "vela_stored_electron_per_m_s": 0.0,
        "vela_stored_hole_per_m_s": 0.0,
        "vela_rebuilt_electron_per_m_s": 0.0,
        "vela_rebuilt_hole_per_m_s": 0.0,
        "sentaurus_alpha_vela_flux_electron_per_m_s": 0.0,
        "sentaurus_alpha_vela_flux_hole_per_m_s": 0.0,
        "vela_alpha_sentaurus_flux_electron_per_m_s": 0.0,
        "vela_alpha_sentaurus_flux_hole_per_m_s": 0.0,
        "sentaurus_alpha_sentaurus_flux_electron_per_m_s": 0.0,
        "sentaurus_alpha_sentaurus_flux_hole_per_m_s": 0.0,
        "sentaurus_alpha_sentaurus_flux_magnitude_electron_per_m_s": 0.0,
        "sentaurus_alpha_sentaurus_flux_magnitude_hole_per_m_s": 0.0,
        "vela_alpha_sentaurus_flux_magnitude_electron_per_m_s": 0.0,
        "vela_alpha_sentaurus_flux_magnitude_hole_per_m_s": 0.0,
        "sentaurus_generation_vela_geometry_per_m_s": 0.0,
        "matched_edge_area_m2": 0.0,
    }
    for row in selected:
        area = float(row["edge_area_proxy_m2"])
        vela_alpha_n = abs(float(row["electron_alpha_m_inv"]))
        vela_alpha_p = abs(float(row["hole_alpha_m_inv"]))
        # The native edge dump stores particle flux in cm^-2 s^-1.
        vela_flux_n = abs(float(row["electron_flux_proxy"])) * 1.0e4
        vela_flux_p = abs(float(row["hole_flux_proxy"])) * 1.0e4
        alpha_n = edge_scalar(row, sent_alpha_n, 100.0)
        alpha_p = edge_scalar(row, sent_alpha_p, 100.0)
        flux_n = edge_vector_projection(row, sent_current_n, 1.0e4) / Q_C
        flux_p = edge_vector_projection(row, sent_current_p, 1.0e4) / Q_C
        flux_n_magnitude = edge_vector_magnitude(row, sent_current_n, 1.0e4) / Q_C
        flux_p_magnitude = edge_vector_magnitude(row, sent_current_p, 1.0e4) / Q_C
        generation = edge_scalar(row, sent_generation, 1.0e6)

        item = {
            "bias_V": args.bias,
            "edge_id": int(row["edge_id"]),
            "node0": int(row["node0"]),
            "node1": int(row["node1"]),
            "x_mid_um": 0.5 * (float(row["x0_um"]) + float(row["x1_um"])),
            "y_mid_um": 0.5 * (float(row["y0_um"]) + float(row["y1_um"])),
            "edge_class": row["edge_class"],
            "edge_area_m2": area,
            "vela_electron_alpha_m_inv": vela_alpha_n,
            "sentaurus_electron_alpha_m_inv": alpha_n,
            "vela_hole_alpha_m_inv": vela_alpha_p,
            "sentaurus_hole_alpha_m_inv": alpha_p,
            "vela_electron_flux_per_m2_s": vela_flux_n,
            "sentaurus_electron_flux_per_m2_s": flux_n,
            "sentaurus_electron_flux_magnitude_per_m2_s": flux_n_magnitude,
            "vela_hole_flux_per_m2_s": vela_flux_p,
            "sentaurus_hole_flux_per_m2_s": flux_p,
            "sentaurus_hole_flux_magnitude_per_m2_s": flux_p_magnitude,
            "vela_stored_source_per_m_s": float(row["edge_source_integral"]) * 1.0e-6,
            "vela_rebuilt_source_per_m_s": (
                vela_alpha_n * vela_flux_n + vela_alpha_p * vela_flux_p
            ) * area,
            "sentaurus_alpha_vela_flux_per_m_s": (
                alpha_n * vela_flux_n + alpha_p * vela_flux_p
            ) * area,
            "vela_alpha_sentaurus_flux_per_m_s": (
                vela_alpha_n * flux_n + vela_alpha_p * flux_p
            ) * area,
            "sentaurus_alpha_sentaurus_flux_per_m_s": (
                alpha_n * flux_n + alpha_p * flux_p
            ) * area,
            "sentaurus_alpha_sentaurus_flux_magnitude_per_m_s": (
                alpha_n * flux_n_magnitude + alpha_p * flux_p_magnitude
            ) * area,
            "sentaurus_generation_vela_geometry_per_m_s": generation * area,
        }
        details.append(item)

        totals["matched_edge_area_m2"] += area
        totals["vela_stored_electron_per_m_s"] += (
            float(row["electron_source_integral"]) * 1.0e-6
        )
        totals["vela_stored_hole_per_m_s"] += (
            float(row["hole_source_integral"]) * 1.0e-6
        )
        totals["vela_rebuilt_electron_per_m_s"] += vela_alpha_n * vela_flux_n * area
        totals["vela_rebuilt_hole_per_m_s"] += vela_alpha_p * vela_flux_p * area
        totals["sentaurus_alpha_vela_flux_electron_per_m_s"] += alpha_n * vela_flux_n * area
        totals["sentaurus_alpha_vela_flux_hole_per_m_s"] += alpha_p * vela_flux_p * area
        totals["vela_alpha_sentaurus_flux_electron_per_m_s"] += vela_alpha_n * flux_n * area
        totals["vela_alpha_sentaurus_flux_hole_per_m_s"] += vela_alpha_p * flux_p * area
        totals["sentaurus_alpha_sentaurus_flux_electron_per_m_s"] += alpha_n * flux_n * area
        totals["sentaurus_alpha_sentaurus_flux_hole_per_m_s"] += alpha_p * flux_p * area
        totals["sentaurus_alpha_sentaurus_flux_magnitude_electron_per_m_s"] += (
            alpha_n * flux_n_magnitude * area
        )
        totals["sentaurus_alpha_sentaurus_flux_magnitude_hole_per_m_s"] += (
            alpha_p * flux_p_magnitude * area
        )
        totals["vela_alpha_sentaurus_flux_magnitude_electron_per_m_s"] += (
            vela_alpha_n * flux_n_magnitude * area
        )
        totals["vela_alpha_sentaurus_flux_magnitude_hole_per_m_s"] += (
            vela_alpha_p * flux_p_magnitude * area
        )
        totals["sentaurus_generation_vela_geometry_per_m_s"] += generation * area

    vela_stored = totals["vela_stored_electron_per_m_s"] + totals["vela_stored_hole_per_m_s"]
    vela_rebuilt = totals["vela_rebuilt_electron_per_m_s"] + totals["vela_rebuilt_hole_per_m_s"]
    sent_alpha_vela_flux = (
        totals["sentaurus_alpha_vela_flux_electron_per_m_s"]
        + totals["sentaurus_alpha_vela_flux_hole_per_m_s"]
    )
    vela_alpha_sent_flux = (
        totals["vela_alpha_sentaurus_flux_electron_per_m_s"]
        + totals["vela_alpha_sentaurus_flux_hole_per_m_s"]
    )
    sent_alpha_sent_flux = (
        totals["sentaurus_alpha_sentaurus_flux_electron_per_m_s"]
        + totals["sentaurus_alpha_sentaurus_flux_hole_per_m_s"]
    )
    sent_alpha_sent_flux_magnitude = (
        totals["sentaurus_alpha_sentaurus_flux_magnitude_electron_per_m_s"]
        + totals["sentaurus_alpha_sentaurus_flux_magnitude_hole_per_m_s"]
    )
    vela_alpha_sent_flux_magnitude = (
        totals["vela_alpha_sentaurus_flux_magnitude_electron_per_m_s"]
        + totals["vela_alpha_sentaurus_flux_magnitude_hole_per_m_s"]
    )
    direct = totals["sentaurus_generation_vela_geometry_per_m_s"]
    summary = {
        "bias_V": args.bias,
        "matched_edges": len(selected),
        **totals,
        "vela_stored_total_per_m_s": vela_stored,
        "vela_rebuilt_total_per_m_s": vela_rebuilt,
        "vela_rebuild_over_stored": vela_rebuilt / vela_stored,
        "sentaurus_alpha_vela_flux_total_per_m_s": sent_alpha_vela_flux,
        "sentaurus_alpha_vela_flux_over_vela": sent_alpha_vela_flux / vela_stored,
        "vela_alpha_sentaurus_flux_total_per_m_s": vela_alpha_sent_flux,
        "vela_alpha_sentaurus_flux_over_vela": vela_alpha_sent_flux / vela_stored,
        "sentaurus_alpha_sentaurus_flux_total_per_m_s": sent_alpha_sent_flux,
        "sentaurus_alpha_sentaurus_flux_over_vela": sent_alpha_sent_flux / vela_stored,
        "sentaurus_alpha_sentaurus_flux_magnitude_total_per_m_s": (
            sent_alpha_sent_flux_magnitude
        ),
        "sentaurus_alpha_sentaurus_flux_magnitude_over_vela": (
            sent_alpha_sent_flux_magnitude / vela_stored
        ),
        "vela_alpha_sentaurus_flux_magnitude_total_per_m_s": (
            vela_alpha_sent_flux_magnitude
        ),
        "vela_alpha_sentaurus_flux_magnitude_over_vela": (
            vela_alpha_sent_flux_magnitude / vela_stored
        ),
        "sentaurus_generation_vela_geometry_total_per_m_s": direct,
        "sentaurus_generation_vela_geometry_over_vela": direct / vela_stored,
    }
    details.sort(key=lambda row: abs(row["vela_stored_source_per_m_s"]), reverse=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.out_dir / "edge_factorization.csv", details)
    write_rows(args.out_dir / "edge_factorization_top100.csv", details[:100])
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(args.out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
