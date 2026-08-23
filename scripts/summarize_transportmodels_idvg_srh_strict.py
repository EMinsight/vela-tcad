#!/usr/bin/env python3
"""Summarize strict TransportModels Id-Vg SRH and sweep diagnostics.

Deep-off points are accepted only when both the logarithmic current comparison
and the port-conservation resolution requirement are meaningful.  A point for
which |Id| < margin * |KCL residual| is explicitly reported as numerically
unresolved rather than assigned a misleading relative-current error.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "idvg_srh_strict_2026-08-21"
)
REFERENCE_ROOT = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "generated/reference_curves"
)
REPORT_DIR = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/reports"
    / "idvg_srh_strict_20260821"
)
MARKDOWN = (
    REPO
    / "docs/validation/transportmodels_idvg_srh_strict_2026-08-21.md"
)
DEEP_OFF_BIASES = (-1.0, -0.84, -0.68, -0.52)
LOG_ERROR_LIMIT_DEX = 0.15
SUBSTRATE_SRH_CLOSURE_LIMIT = 0.01


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str) -> float:
    return float(row[key])


def nearest(rows: list[dict[str, str]], key: str, target: float) -> dict[str, str]:
    result = min(rows, key=lambda row: abs(number(row, key) - target))
    if abs(number(result, key) - target) > 1.0e-10:
        raise ValueError(f"No row at {key}={target} in available data")
    return result


def log_error(candidate: float, reference: float) -> float | None:
    if candidate == 0.0 or reference == 0.0:
        return None
    return abs(math.log10(abs(candidate)) - math.log10(abs(reference)))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def deep_off_rows(model: str) -> list[dict[str, Any]]:
    run = read_rows(RUN_ROOT / f"{model}_forward" / "curve.csv")
    balance = read_rows(RUN_ROOT / f"{model}_forward" / "srh_balance.csv")
    reference = read_rows(
        REFERENCE_ROOT / f"transportmodels_sentaurus2022_{model}_idvg_reference.csv"
    )
    result: list[dict[str, Any]] = []
    for bias in DEEP_OFF_BIASES:
        current_row = nearest(run, "bias_V", bias)
        balance_row = nearest(balance, "bias_V", bias)
        reference_row = nearest(reference, "bias_V", bias)
        candidate = number(current_row, "current_total_A_per_um")
        sentaurus = number(reference_row, "current_total")
        error = log_error(candidate, sentaurus)
        kcl_residual = number(
            balance_row, "four_terminal_kcl_residual_A_per_um"
        )
        id_to_kcl = (
            abs(candidate) / abs(kcl_residual)
            if kcl_residual != 0.0
            else (math.inf if candidate != 0.0 else 0.0)
        )
        resolution_ok = candidate != 0.0 and id_to_kcl >= 10.0
        closure_error = number(
            balance_row, "substrate_generation_magnitude_relative_error"
        )
        log_ok = error is not None and error <= LOG_ERROR_LIMIT_DEX
        closure_ok = closure_error <= SUBSTRATE_SRH_CLOSURE_LIMIT
        accepted = resolution_ok and log_ok and closure_ok
        result.append(
            {
                "model": model.upper(),
                "bias_V": bias,
                "sentaurus_Id_A_per_um": sentaurus,
                "vela_Id_A_per_um": candidate,
                "log10_abs_Id_error_dex": "" if error is None else error,
                "srh_generation_A_per_um": number(
                    balance_row, "srh_generation_current_A_per_um"
                ),
                "substrate_hole_A_per_um": number(
                    balance_row, "substrate_hole_current_A_per_um"
                ),
                "substrate_srh_closure_relative_error": closure_error,
                "four_terminal_kcl_residual_A_per_um": kcl_residual,
                "Id_to_KCL_ratio": id_to_kcl,
                "numerical_status": (
                    "resolved" if resolution_ok else "numerically_unresolved"
                ),
                "log_error_pass": log_ok,
                "substrate_srh_closure_pass": closure_ok,
                "deep_off_acceptance": "pass" if accepted else "not_accepted",
            }
        )
    return result


def dd_roundtrip_rows() -> list[dict[str, Any]]:
    rows = [
        row
        for row in read_rows(RUN_ROOT / "dd_roundtrip" / "curve.csv")
        if row.get("converged") == "1"
    ]
    groups: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[round(number(row, "bias_V"), 12)].append(row)
    result: list[dict[str, Any]] = []
    for bias, matches in sorted(groups.items(), reverse=True):
        if len(matches) < 2:
            continue
        forward = number(matches[0], "current_total_A_per_um")
        reverse = number(matches[-1], "current_total_A_per_um")
        relative = abs(reverse - forward) / max(abs(forward), 1.0e-300)
        result.append(
            {
                "bias_V": bias,
                "forward_Id_A_per_um": forward,
                "reverse_Id_A_per_um": reverse,
                "relative_difference": relative,
                "log10_abs_difference_dex": abs(
                    math.log10(max(abs(reverse), 1.0e-300))
                    - math.log10(max(abs(forward), 1.0e-300))
                ),
            }
        )
    return result


def fmt(value: Any, digits: int = 4) -> str:
    if value == "" or value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def main() -> int:
    deep = deep_off_rows("dd") + deep_off_rows("dg")
    roundtrip = dd_roundtrip_rows()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(REPORT_DIR / "deep_off_acceptance.csv", deep)
    write_csv(REPORT_DIR / "dd_roundtrip_consistency.csv", roundtrip)

    summary = {
        "criteria": {
            "deep_off_log_error_limit_dex": LOG_ERROR_LIMIT_DEX,
            "substrate_srh_closure_relative_error_limit": SUBSTRATE_SRH_CLOSURE_LIMIT,
            "resolution_rule": "abs(Id) >= 10 * abs(four-terminal KCL residual)",
        },
        "deep_off_points": deep,
        "dd_roundtrip_overlap": roundtrip,
        "dd_forward_completed": True,
        "dd_reverse_completed": False,
        "dd_reverse_last_converged_bias_V": 1.72,
        "dg_forward_completed": False,
        "dg_last_converged_bias_V": -0.4,
        "dg_failure_bias_V": -0.3987109375,
        "dg_failure_residual_norm": 1.43739005196077e-10,
    }
    (REPORT_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# TransportModels Id–Vg strict SRH and deep-off validation",
        "",
        "## Acceptance rules",
        "",
        "- A point is numerically resolved only when `|Id| >= 10 * |four-terminal KCL residual|`.",
        f"- Resolved deep-off points require log-current error <= {LOG_ERROR_LIMIT_DEX} dex.",
        f"- The silicon SRH-generation/substrate-hole closure error must be <= {SUBSTRATE_SRH_CLOSURE_LIMIT:.0%}.",
        "- Unresolved points are not assigned a passing relative-current comparison.",
        "",
        "## Deep-off results",
        "",
        "| Model | Vg (V) | Sentaurus Id (A/um) | Vela Id (A/um) | log error (dex) | SRH/substrate closure | Id/KCL | Status | Acceptance |",
        "|---|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in deep:
        lines.append(
            "| {model} | {bias} | {sent} | {vela} | {log} | {closure} | {ratio} | {status} | {acceptance} |".format(
                model=row["model"],
                bias=fmt(row["bias_V"]),
                sent=fmt(row["sentaurus_Id_A_per_um"]),
                vela=fmt(row["vela_Id_A_per_um"]),
                log=fmt(row["log10_abs_Id_error_dex"]),
                closure=fmt(row["substrate_srh_closure_relative_error"]),
                ratio=fmt(row["Id_to_KCL_ratio"]),
                status=row["numerical_status"],
                acceptance=row["deep_off_acceptance"],
            )
        )
    lines.extend(
        [
            "",
            "## DD direction check",
            "",
            "| Vg (V) | Forward Id (A/um) | Reverse Id (A/um) | Relative difference |",
            "|---:|---:|---:|---:|",
        ]
    )
    for row in roundtrip:
        lines.append(
            f"| {fmt(row['bias_V'])} | {fmt(row['forward_Id_A_per_um'])} | "
            f"{fmt(row['reverse_Id_A_per_um'])} | {fmt(row['relative_difference'])} |"
        )
    lines.extend(
        [
            "",
            "The DD forward sweep completed to 2.2 V. The same-process reverse sweep matched the overlapping 2.04, 1.88 and 1.72 V points, then stopped near 1.623 V at the strict residual floor.",
            "",
            "The DG sweep resolved the requested deep-off points through -0.52 V and was bridged to -0.40 V. It then stopped at -0.3987109375 V with residual 1.4374e-10 after adaptive step reduction; later DG biases are therefore not claimed as completed.",
            "",
        ]
    )
    MARKDOWN.parent.mkdir(parents=True, exist_ok=True)
    MARKDOWN.write_text("\n".join(lines), encoding="utf-8")
    print(REPORT_DIR / "summary.json")
    print(MARKDOWN)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
