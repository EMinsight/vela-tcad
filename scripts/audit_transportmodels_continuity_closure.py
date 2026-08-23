#!/usr/bin/env python3
"""Audit TransportModels contact-flux versus SRH continuity closure."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/reports"
    / "idvg_deep_off_precision_20260822"
)
OUTPUT = ROOT / "continuity_closure_audit"

CLOSURE_CASES = (
    (
        "DD",
        -1.0,
        "hard_gate_full_step",
        ROOT / "dd_m1p000000_hard_gate/curve.csv",
    ),
    (
        "DD",
        -0.68,
        "hard_gate_full_step",
        ROOT / "dd_m0p680000_hard_gate/curve.csv",
    ),
    (
        "DD",
        -0.68,
        "hard_gate_damping_0p5",
        ROOT
        / "newton_calibration/floor2e11_qf1e2_damp5e1"
        / "dd_m0p680000/curve.csv",
    ),
    (
        "DG",
        -0.52,
        "hard_gate_full_step",
        ROOT / "dg_m0p520000_hard_gate/curve.csv",
    ),
)

PHYSICAL_CASES = (
    ("DD", -0.68, ROOT / "dd_m0p680000_qfref"),
    ("DG", -0.52, ROOT / "dg_m0p520000_qfref"),
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], name: str) -> float:
    return float(row[name])


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def closure_rows() -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for model, bias, variant, path in CLOSURE_CASES:
        row = read_rows(path)[-1]
        electron_contact = number(row, "global_electron_contact_flux")
        hole_contact = number(row, "global_hole_contact_flux")
        electron_source = number(row, "global_electron_integrated_source")
        hole_source = number(row, "global_hole_integrated_source")
        result.append(
            {
                "model": model,
                "bias_V": bias,
                "variant": variant,
                "converged": int(row["converged"]),
                "carrier_row_violations": int(row["carrier_row_violations"]),
                "carrier_row_max_ratio": number(row, "carrier_row_max_ratio"),
                "electron_contact_flux_scaled": electron_contact,
                "electron_integrated_source_scaled": electron_source,
                "electron_mismatch_scaled": electron_contact - electron_source,
                "electron_closure_ratio": number(
                    row, "global_electron_continuity_closure_ratio"
                ),
                "hole_contact_flux_scaled": hole_contact,
                "hole_integrated_source_scaled": hole_source,
                "hole_mismatch_scaled": hole_contact - hole_source,
                "hole_closure_ratio": number(
                    row, "global_hole_continuity_closure_ratio"
                ),
                "pair_source_relative_difference": abs(
                    electron_source - hole_source
                )
                / max(abs(electron_source), abs(hole_source), 1.0e-300),
                "failure_reason": row["failure_reason"],
            }
        )
    return result


def physical_case_summary(
    model: str, bias: float, directory: Path
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    config = json.loads((directory / "config.json").read_text(encoding="utf-8"))
    mesh = json.loads(Path(config["mesh_file"]).read_text(encoding="utf-8"))
    doping_rows = read_rows(Path(config["node_doping_file"]))
    net_doping = {
        int(row["node_id"]): float(row["donors_cm3"]) - float(row["acceptors_cm3"])
        for row in doping_rows
    }
    contact_bias = {row["name"]: float(row["bias"]) for row in config["contacts"]}
    strongest_contact = ""
    strongest_doping = 0.0
    for contact in mesh["contacts"]:
        mean_net = sum(net_doping[int(node)] for node in contact["node_ids"]) / len(
            contact["node_ids"]
        )
        # Mirror NewtonSolver::configureQuasiFermiReferences: equal later
        # contacts do not replace the first strongest n-type contact.
        if mean_net > strongest_doping:
            strongest_doping = mean_net
            strongest_contact = contact["name"]

    terminal = read_rows(directory / "terminal_balance.csv")
    srh = read_rows(directory / "srh_balance.csv")[-1]
    edges = read_rows(directory / "contact_edges.csv")

    electron_sum = sum(number(row, "current_electron_A_per_um") for row in terminal)
    hole_sum = sum(number(row, "current_hole_A_per_um") for row in terminal)
    total_sum = sum(number(row, "current_total_A_per_um") for row in terminal)
    generation = number(srh, "srh_generation_current_A_per_um")

    contact_rows: list[dict[str, Any]] = []
    for terminal_row in terminal:
        contact = terminal_row["contact"]
        contact_edges = [row for row in edges if row["current_contact"] == contact]
        active_edges = [row for row in contact_edges if number(row, "mun") > 0.0]
        edge_electron_A_per_um = sum(
            number(row, "current_electron") for row in active_edges
        ) * 1.0e-6
        contact_rows.append(
            {
                "model": model,
                "bias_V": bias,
                "contact": contact,
                "active_electron_edges": len(active_edges),
                "zero_electron_current_edges": sum(
                    number(row, "current_electron") == 0.0 for row in active_edges
                ),
                "zero_reported_phin_drop_edges": sum(
                    number(row, "phin0") == number(row, "phin1")
                    for row in active_edges
                ),
                "edge_electron_current_A_per_um": edge_electron_A_per_um,
                "terminal_electron_current_A_per_um": number(
                    terminal_row, "current_electron_A_per_um"
                ),
                "edge_terminal_relative_difference": abs(
                    edge_electron_A_per_um
                    - number(terminal_row, "current_electron_A_per_um")
                )
                / max(
                    abs(number(terminal_row, "current_electron_A_per_um")),
                    1.0e-300,
                ),
            }
        )

    summary = {
        "model": model,
        "bias_V": bias,
        "electron_reference_contact": strongest_contact,
        "electron_reference_V": contact_bias[strongest_contact],
        "strongest_contact_net_doping_cm3": strongest_doping,
        "electron_terminal_sum_A_per_um": electron_sum,
        "hole_terminal_sum_A_per_um": hole_sum,
        "total_terminal_kcl_A_per_um": total_sum,
        "srh_generation_A_per_um": generation,
        "electron_generation_magnitude_relative_error": abs(
            abs(electron_sum) - generation
        )
        / max(generation, 1.0e-300),
        "hole_generation_magnitude_relative_error": abs(
            abs(hole_sum) - generation
        )
        / max(generation, 1.0e-300),
        "kcl_to_generation_ratio": abs(total_sum) / max(generation, 1.0e-300),
    }
    return summary, contact_rows


def markdown_report(
    closures: list[dict[str, Any]],
    physical: list[dict[str, Any]],
    contacts: list[dict[str, Any]],
) -> str:
    lines = [
        "# TransportModels 接触通量—SRH 连续性闭合审计（2026-08-22）",
        "",
        "## 审计结论",
        "",
        "SRH 体积分、符号和控制体权重不是当前主因。所有审计点的电子/空穴积分源完全相同，空穴接触通量能以 1e-10 或更好精度闭合该源；失败仅集中在电子接触通量。",
        "",
        "DD -0.68 V 的电子接触通量对 Newton 阻尼敏感，使用全步和 0.5 阻尼时甚至改变符号，而积分源保持稳定。这是极低电子电流的数值分辨率问题，不是物理源项变化。",
        "",
        "## 求解器内部闭合量",
        "",
        "| 模型 | Vg/V | 变体 | 局部违规 | 电子闭合比 | 空穴闭合比 | 电子接触通量 | 电子积分源 |",
        "|---|---:|---|---:|---:|---:|---:|---:|",
    ]
    for row in closures:
        lines.append(
            "| {model} | {bias_V:.2f} | {variant} | {carrier_row_violations} | "
            "{electron_closure_ratio:.3e} | {hole_closure_ratio:.3e} | "
            "{electron_contact_flux_scaled:.3e} | {electron_integrated_source_scaled:.3e} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "DD -0.68 V 的 0.5 阻尼变体已经没有局部载流子行违规，但电子闭合比仍为约 1.03；因此局部相对残差门槛不能保证极低总电流的全局绝对闭合。",
            "",
            "## 物理端口电流与 SRH 产生",
            "",
            "| 模型 | Vg/V | 电子 reference | ΣIe/(A/um) | ΣIh/(A/um) | SRH 产生/(A/um) | 电子误差 | 空穴误差 | |KCL|/SRH |",
            "|---|---:|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in physical:
        lines.append(
            "| {model} | {bias_V:.2f} | {electron_reference_contact} ({electron_reference_V:.1f} V) | {electron_terminal_sum_A_per_um:.3e} | "
            "{hole_terminal_sum_A_per_um:.3e} | {srh_generation_A_per_um:.3e} | "
            "{electron_generation_magnitude_relative_error:.3e} | "
            "{hole_generation_magnitude_relative_error:.3e} | "
            "{kcl_to_generation_ratio:.3e} |".format(**row)
        )
    lines.extend(
        [
            "",
            "空穴端口和 SRH 产生量吻合，说明 SRH 广义 Fermi 公式、积分符号和硅区控制体体积相互一致。电子端口和 KCL 则未达到深关断所需精度。",
            "",
            "## 接触边审计",
            "",
            "| 模型 | Vg/V | 接触 | 活跃电子边 | 零电子电流边 | 报告为零的 phin 差边 | Ie/(A/um) |",
            "|---|---:|---|---:|---:|---:|---:|",
        ]
    )
    for row in contacts:
        lines.append(
            "| {model} | {bias_V:.2f} | {contact} | {active_electron_edges} | "
            "{zero_electron_current_edges} | {zero_reported_phin_drop_edges} | "
            "{terminal_electron_current_A_per_um:.3e} |".format(**row)
        )
    lines.extend(
        [
            "",
            "DD -0.68 V 的源端 30 条活跃电子边全部给出零电流，而漏端 30 条边在物理 phin 输出差同样为零时仍通过内部增量产生非零电流。该不对称性与单一全局准费米 reference 以及另一端约 1.1 V 基值附近的 ULP/线性求解分辨率有关。",
            "",
            "## 实现核对",
            "",
            "- 全局闭合使用无边界行替换的物理连续性项；接触节点只累计 SG 通量，自由节点累计 SRH/碰撞电离源。该离散恒等式和现有单元测试的符号约定正确。",
            "- Newton 实际边界行仍为 Dirichlet identity；局部载流子行门槛不会检查被替换的接触行，因此需要独立全局闭合门槛。",
            "- 当前重启 CSV 只保存物理 phin/phip，不保存 phinIncrement/phipIncrement 和 reference；重新读取会丢失低于绝对电势 ULP 的内部增量，不能用于精确复放此问题。",
            "",
            "## 后续建议",
            "",
            "1. 先扩展重启/失败状态格式，持久化准费米 increment 和 reference，建立同一失败状态的逐边可重复审计。",
            "2. 比较单一全局 reference 与按接触盆地/连通分区 reference，避免源漏相差 1.1 V 时只有一端处于零附近。",
            "3. 在不改变物理模型的前提下，对电子接触通量及全局自由行求和增加 long-double/补偿求和 A/B；只有同一状态闭合显著改善时才进入生产路径。",
            "4. 保留现有全局硬门槛，不通过调大 stall floor 接受这些点。",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    closures = closure_rows()
    physical: list[dict[str, Any]] = []
    contacts: list[dict[str, Any]] = []
    for model, bias, directory in PHYSICAL_CASES:
        summary, contact_rows = physical_case_summary(model, bias, directory)
        physical.append(summary)
        contacts.extend(contact_rows)

    write_rows(OUTPUT / "solver_closure.csv", closures)
    write_rows(OUTPUT / "physical_balance.csv", physical)
    write_rows(OUTPUT / "contact_edge_balance.csv", contacts)
    payload = {
        "solver_closure": closures,
        "physical_balance": physical,
        "contact_edge_balance": contacts,
    }
    (OUTPUT / "audit.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    report = markdown_report(closures, physical, contacts)
    (OUTPUT / "report.md").write_text(report, encoding="utf-8")
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
