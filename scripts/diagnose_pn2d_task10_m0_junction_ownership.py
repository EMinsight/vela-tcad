#!/usr/bin/env python3
"""Run controlled PN2D M0 junction-node doping variants.

The Task 10 dose-preserving mesh construction changed only the three nodes on
the x=1 um metallurgical junction.  This diagnostic keeps the mesh, physics,
solver, and voltage lattice fixed while varying the donor/acceptor ownership
of those nodes.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any


VARIANTS = {
    "old_double_species": (1.0, 1.0),
    "new_n_owned": (1.0, 0.0),
    "balanced_half": (0.5, 0.5),
    "neutral_zero": (0.0, 0.0),
    "p_owned": (0.0, 1.0),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def resolve_config_path(config_path: Path, raw: str) -> Path:
    path = Path(raw)
    if path.is_absolute():
        return path.resolve()
    return (config_path.parent / path).resolve()


def junction_nodes(mesh_path: Path, junction_x_um: float) -> list[int]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = [
        int(node["id"])
        for node in mesh["nodes"]
        if math.isclose(
            float(node["x"]),
            junction_x_um,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
    ]
    if not nodes:
        raise RuntimeError(f"no mesh nodes found at x={junction_x_um:g} um")
    return sorted(nodes)


def nodal_control_areas_um2(mesh_path: Path) -> dict[int, float]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    areas = {node_id: 0.0 for node_id in coordinates}
    cells = mesh.get("cells", mesh.get("triangles"))
    if cells is None:
        raise RuntimeError("mesh JSON has neither cells nor triangles")
    for cell in cells:
        node_ids = [int(value) for value in cell["node_ids"]]
        if len(node_ids) != 3:
            raise RuntimeError("dose audit requires triangular cells")
        (x0, y0), (x1, y1), (x2, y2) = [
            coordinates[node_id] for node_id in node_ids
        ]
        area = 0.5 * abs(
            (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        )
        for node_id in node_ids:
            areas[node_id] += area / 3.0
    return areas


def doping_metrics(
    path: Path,
    nodal_areas: dict[int, float],
    junction_ids: list[int],
) -> dict[str, Any]:
    rows = read_csv(path)
    values = {
        int(row["node_id"]): (
            float(row["donors_cm3"]),
            float(row["acceptors_cm3"]),
        )
        for row in rows
    }
    total_dose = sum(
        nodal_areas[node_id] * (donor + acceptor)
        for node_id, (donor, acceptor) in values.items()
    )
    signed_net_dose = sum(
        nodal_areas[node_id] * (donor - acceptor)
        for node_id, (donor, acceptor) in values.items()
    )
    absolute_net_dose = sum(
        nodal_areas[node_id] * abs(donor - acceptor)
        for node_id, (donor, acceptor) in values.items()
    )
    junction = [values[node_id] for node_id in junction_ids]
    return {
        "total_impurity_dose_cm3_um2": total_dose,
        "signed_net_dose_cm3_um2": signed_net_dose,
        "absolute_net_dose_cm3_um2": absolute_net_dose,
        "junction_donors_cm3": [item[0] for item in junction],
        "junction_acceptors_cm3": [item[1] for item in junction],
        "junction_net_doping_cm3": [item[0] - item[1] for item in junction],
        "junction_total_impurity_cm3": [
            item[0] + item[1] for item in junction
        ],
    }


def variant_doping(
    source: Path,
    destination: Path,
    node_ids: list[int],
    donor_fraction: float,
    acceptor_fraction: float,
    concentration_cm3: float,
) -> None:
    rows = read_csv(source)
    node_set = set(node_ids)
    seen: set[int] = set()
    for row in rows:
        node_id = int(row["node_id"])
        if node_id not in node_set:
            continue
        row["donors_cm3"] = f"{donor_fraction * concentration_cm3:.17g}"
        row["acceptors_cm3"] = f"{acceptor_fraction * concentration_cm3:.17g}"
        seen.add(node_id)
    if seen != node_set:
        raise RuntimeError("junction nodes are missing from the doping CSV")
    write_csv(
        destination,
        rows,
        ["node_id", "donors_cm3", "acceptors_cm3"],
    )


def configure_case(
    base: dict[str, Any],
    case_dir: Path,
    doping_path: Path,
    minimum_step_V: float,
) -> dict[str, Any]:
    config = copy.deepcopy(base)
    config["node_doping_file"] = str(doping_path.resolve())
    config["output_csv"] = str((case_dir / "iv.csv").resolve())
    sweep = config["sweep"]
    sweep["min_step"] = minimum_step_V
    sweep["stop_on_failure"] = True
    sweep["write_state_file"] = str((case_dir / "last_state.csv").resolve())
    sweep["write_state_every_point_prefix"] = str(
        (case_dir / "states" / "state").resolve()
    )
    sweep["write_vtk"] = False
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str((case_dir / "newton_history.csv").resolve()),
        "attempts_csv_file": str((case_dir / "newton_attempts.csv").resolve()),
        "iterations_csv_file": str(
            (case_dir / "newton_iterations.csv").resolve()
        ),
        "rejected_state_directory": str(
            (case_dir / "rejected_states").resolve()
        ),
    }
    diagnostics["bv_process_probe"] = {
        "enabled": True,
        "csv_file": str((case_dir / "process_probe.csv").resolve()),
    }
    return config


def float_at_bias(
    rows: list[dict[str, str]], bias: float, field: str
) -> float | None:
    for row in rows:
        if math.isclose(
            float(row["bias_V"]), bias, rel_tol=0.0, abs_tol=1.0e-9
        ):
            return float(row[field])
    return None


def source_at_bias(path: Path, bias: float) -> float | None:
    if not path.is_file():
        return None
    values = [
        float(row["source_integral"])
        for row in read_csv(path)
        if math.isclose(
            float(row["bias_V"]), bias, rel_tol=0.0, abs_tol=1.0e-9
        )
    ]
    return sum(values) if values else None


def bias_token(bias: float) -> str:
    sign = "m" if bias < 0.0 else ""
    return sign + f"{abs(bias):.6f}".replace(".", "p")


def junction_state_metrics(
    case_dir: Path, bias: float, junction_ids: list[int]
) -> dict[str, Any] | None:
    path = (
        case_dir
        / "states"
        / f"state_bias_{bias_token(bias)}.csv"
    )
    if not path.is_file():
        return None
    selected = [
        row for row in read_csv(path) if int(row["node_id"]) in junction_ids
    ]
    if len(selected) != len(junction_ids):
        raise RuntimeError("junction state is incomplete")

    def values(field: str) -> list[float]:
        return [float(row[field]) for row in selected]

    return {
        "state_csv": str(path),
        "state_sha256": sha256(path),
        "psi_V": values("psi"),
        "phin_V": values("phin"),
        "phip_V": values("phip"),
        "electrons_m3": values("electrons_m3"),
        "holes_m3": values("holes_m3"),
    }


def process_metrics(case_dir: Path, bias: float) -> dict[str, Any] | None:
    path = case_dir / "process_probe.csv"
    if not path.is_file():
        return None
    rows = [
        row
        for row in read_csv(path)
        if math.isclose(
            float(row["bias_V"]), bias, rel_tol=0.0, abs_tol=1.0e-9
        )
    ]
    if not rows:
        return None

    def magnitude(row: dict[str, str], x_field: str, y_field: str) -> float:
        return math.hypot(float(row[x_field]), float(row[y_field]))

    top_source = max(rows, key=lambda row: abs(float(row["source_integral"])))
    by_carrier: dict[str, dict[str, Any]] = {}
    for carrier in sorted({row["carrier"] for row in rows}):
        subset = [row for row in rows if row["carrier"] == carrier]
        by_carrier[carrier] = {
            "record_count": len(subset),
            "max_electric_field": max(
                magnitude(row, "electric_field_x", "electric_field_y")
                for row in subset
            ),
            "max_qf_gradient": max(
                magnitude(row, "qf_gradient_x", "qf_gradient_y")
                for row in subset
            ),
            "max_impact_field": max(
                abs(float(row["impact_field"])) for row in subset
            ),
            "max_alpha": max(abs(float(row["alpha"])) for row in subset),
            "max_selected_flux_magnitude": max(
                abs(float(row["selected_flux_magnitude"])) for row in subset
            ),
            "source_integral": sum(
                float(row["source_integral"]) for row in subset
            ),
        }
    return {
        "process_probe_csv": str(path),
        "process_probe_sha256": sha256(path),
        "record_count": len(rows),
        "source_integral_total": sum(
            float(row["source_integral"]) for row in rows
        ),
        "by_carrier": by_carrier,
        "top_abs_source_record": {
            key: top_source[key]
            for key in (
                "carrier",
                "cell_id",
                "local_edge",
                "edge_id",
                "node0",
                "node1",
                "density0",
                "density1",
                "midpoint_density",
                "impact_field",
                "alpha",
                "selected_flux_magnitude",
                "generation_rate",
                "source_integral",
                "active_branch_fingerprint",
            )
        },
    }


def summarize_case(
    name: str,
    case_dir: Path,
    doping_path: Path,
    junction_ids: list[int],
    nodal_areas: dict[int, float],
) -> dict[str, Any]:
    iv_path = case_dir / "iv.csv"
    attempts_path = case_dir / "newton_attempts.csv"
    iv = read_csv(iv_path) if iv_path.is_file() else []
    attempts = read_csv(attempts_path) if attempts_path.is_file() else []
    last = iv[-1] if iv else {}
    rejected = [row for row in attempts if row["status"] == "rejected"]
    return {
        "variant": name,
        "junction_nodes": junction_ids,
        "doping_csv": str(doping_path),
        "doping_sha256": sha256(doping_path),
        "returncode": None,
        "iv_rows": len(iv),
        "last_bias_V": float(last["bias_V"]) if last else None,
        "last_converged": last.get("converged") == "1",
        "last_failure_reason": last.get("failure_reason", ""),
        "attempt_count": len(attempts),
        "rejected_attempt_count": len(rejected),
        "current_total_A_per_um_at_m17V": float_at_bias(
            iv, -17.0, "current_total_A_per_um"
        ),
        "source_integral_at_m17V": source_at_bias(
            case_dir / "process_probe.csv", -17.0
        ),
        "doping_metrics": doping_metrics(
            doping_path, nodal_areas, junction_ids
        ),
        "junction_state_at_m17V": junction_state_metrics(
            case_dir, -17.0, junction_ids
        ),
        "process_metrics_at_m17V": process_metrics(case_dir, -17.0),
        "reached_m20V": (
            float_at_bias(iv, -20.0, "current_total_A_per_um") is not None
            and any(
                math.isclose(
                    float(row["bias_V"]),
                    -20.0,
                    rel_tol=0.0,
                    abs_tol=1.0e-9,
                )
                and row["converged"] == "1"
                for row in iv
            )
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--junction-x-um", type=float, default=1.0)
    parser.add_argument("--concentration-cm3", type=float, default=1.0e17)
    parser.add_argument("--minimum-step-V", type=float, default=1.0e-5)
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="Reuse complete case artifacts instead of rerunning the solver",
    )
    parser.add_argument(
        "--variants",
        default=",".join(VARIANTS),
        help="Comma-separated subset of controlled variants",
    )
    args = parser.parse_args()

    runner = args.runner.resolve()
    base_config_path = args.base_config.resolve()
    root = args.out_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    base = json.loads(base_config_path.read_text(encoding="utf-8"))
    mesh_path = resolve_config_path(base_config_path, base["mesh_file"])
    source_doping = resolve_config_path(
        base_config_path, base["node_doping_file"]
    )
    junction_ids = junction_nodes(mesh_path, args.junction_x_um)
    nodal_areas = nodal_control_areas_um2(mesh_path)
    requested = [item.strip() for item in args.variants.split(",") if item.strip()]
    unknown = [name for name in requested if name not in VARIANTS]
    if unknown:
        raise SystemExit(f"unknown variants: {', '.join(unknown)}")

    summaries = []
    for name in requested:
        case_dir = root / name
        case_dir.mkdir(parents=True, exist_ok=True)
        (case_dir / "states").mkdir(parents=True, exist_ok=True)
        donor_fraction, acceptor_fraction = VARIANTS[name]
        doping_path = case_dir / "doping.csv"
        variant_doping(
            source_doping,
            doping_path,
            junction_ids,
            donor_fraction,
            acceptor_fraction,
            args.concentration_cm3,
        )
        config = configure_case(
            base, case_dir, doping_path, args.minimum_step_V
        )
        config_path = case_dir / "simulation.json"
        config_path.write_text(
            json.dumps(config, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if args.reuse and (case_dir / "iv.csv").is_file():
            returncode = 0
        else:
            with (case_dir / "runner.log").open("w", encoding="utf-8") as log:
                completed = subprocess.run(
                    [str(runner), "--config", str(config_path)],
                    cwd=case_dir,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            returncode = completed.returncode
        summary = summarize_case(
            name, case_dir, doping_path, junction_ids, nodal_areas
        )
        summary["returncode"] = returncode
        summaries.append(summary)

    output = {
        "schema": "vela.pn2d_task10_m0_junction_ownership.v1",
        "runner": str(runner),
        "base_config": str(base_config_path),
        "mesh_file": str(mesh_path),
        "mesh_sha256": sha256(mesh_path),
        "source_doping": str(source_doping),
        "junction_x_um": args.junction_x_um,
        "junction_nodes": junction_ids,
        "minimum_step_V": args.minimum_step_V,
        "variants": summaries,
    }
    summary_path = root / "summary.json"
    summary_path.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
