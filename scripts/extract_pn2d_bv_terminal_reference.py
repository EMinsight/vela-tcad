#!/usr/bin/env python3
"""Extract PN2D Sentaurus BV terminal current from a PLT file."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Sequence


SCHEMA = "vela.pn2d_bv_terminal_reference.v1"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-plt", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--out-meta", type=Path, default=None)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument("--bias-column", default=None)
    parser.add_argument("--current-column", default=None)
    parser.add_argument("--expected-min-bias", type=float, default=-20.0)
    parser.add_argument("--sign", choices=("native", "flip"), default="native")
    return parser.parse_args(argv)


def fmt(value: float) -> str:
    if value == 0.0:
        return "0"
    return f"{value:.12g}"


def parse_quoted_list(text: str, key: str) -> list[str]:
    match = re.search(rf"{re.escape(key)}\s*=\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        return []
    return re.findall(r'"([^"]+)"', match.group(1))


def parse_values_block(text: str, column_count: int) -> list[list[float]]:
    match = re.search(r"Values\s*=\s*\[(.*?)\]", text, re.DOTALL)
    if not match:
        match = re.search(r"Data\s*\{(.*?)\}", text, re.DOTALL | re.IGNORECASE)
    if not match:
        raise ValueError("PLT contains neither Values=[...] nor Data {...} block")
    numbers = [
        float(value)
        for value in re.findall(
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
            match.group(1),
        )
    ]
    if column_count <= 0:
        raise ValueError("PLT dataset list is empty")
    if len(numbers) % column_count != 0:
        raise ValueError(
            f"PLT numeric value count {len(numbers)} is not divisible by dataset count {column_count}"
        )
    return [numbers[index:index + column_count] for index in range(0, len(numbers), column_count)]


def load_plt_rows(path: Path, bias_column: str, current_column: str, sign: str) -> list[dict[str, float]]:
    if not path.exists():
        raise FileNotFoundError(f"Sentaurus PLT not found: {path}")
    text = path.read_text(errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    if bias_column not in datasets:
        raise ValueError(f"PLT bias column not found: {bias_column}")
    if current_column not in datasets:
        raise ValueError(f"PLT current column not found: {current_column}")
    bias_index = datasets.index(bias_column)
    current_index = datasets.index(current_column)
    multiplier = -1.0 if sign == "flip" else 1.0
    rows: list[dict[str, float]] = []
    seen: set[float] = set()
    for raw in parse_values_block(text, len(datasets)):
        bias = raw[bias_index]
        current = multiplier * raw[current_index]
        if not (math.isfinite(bias) and math.isfinite(current)):
            continue
        key = round(bias, 12)
        if key in seen:
            continue
        seen.add(key)
        rows.append({"bias_V": bias, "current_total": current})
    return sorted(rows, key=lambda row: row["bias_V"], reverse=True)


def validate_rows(rows: list[dict[str, float]], expected_min_bias: float) -> None:
    if not rows:
        raise ValueError("Sentaurus terminal current reference has no finite rows")
    finite = [row for row in rows if math.isfinite(row["current_total"])]
    if len(finite) != len(rows):
        raise ValueError("Sentaurus terminal current reference contains non-finite current")
    if all(row["current_total"] == 0.0 for row in rows):
        raise ValueError("Sentaurus terminal current reference current_total column is all zero")
    min_bias = min(row["bias_V"] for row in rows)
    max_bias = max(row["bias_V"] for row in rows)
    if min_bias > expected_min_bias + 1.0e-9:
        raise ValueError(
            f"Sentaurus terminal current reference does not cover expected minimum bias "
            f"{expected_min_bias:g} V; range is {min_bias:g}..{max_bias:g} V"
        )


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["bias_V", "current_total", "current_total_unit", "current_total_A_per_um"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "bias_V": fmt(row["bias_V"]),
                "current_total": fmt(row["current_total"]),
                "current_total_unit": "A",
                "current_total_A_per_um": "",
            })


def write_meta(
    path: Path,
    *,
    sentaurus_plt: Path,
    rows: list[dict[str, float]],
    contact: str,
    bias_column: str,
    current_column: str,
    sign: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "source": str(sentaurus_plt),
        "contact": contact,
        "bias_column": bias_column,
        "current_column": current_column,
        "current_total_unit": "A",
        "current_total_A_per_um": None,
        "current_sign_convention": sign,
        "row_count": len(rows),
        "bias_range_V": [min(row["bias_V"] for row in rows), max(row["bias_V"] for row in rows)],
        "notes": [
            "Sentaurus PLT TotalCurrent is stored in native A.",
            "For 2D per-depth comparisons, convert explicitly outside this fixture; do not silently treat A as A/um.",
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    contact = args.contact
    bias_column = args.bias_column or f"{contact} OuterVoltage"
    current_column = args.current_column or f"{contact} TotalCurrent"
    try:
        rows = load_plt_rows(args.sentaurus_plt, bias_column, current_column, args.sign)
        validate_rows(rows, args.expected_min_bias)
        write_csv(args.out_csv, rows)
        meta_path = args.out_meta or args.out_csv.with_suffix(args.out_csv.suffix + ".meta.json")
        write_meta(
            meta_path,
            sentaurus_plt=args.sentaurus_plt,
            rows=rows,
            contact=contact,
            bias_column=bias_column,
            current_column=current_column,
            sign=args.sign,
        )
    except Exception as exc:  # noqa: BLE001 - CLI should report clean failures.
        print(str(exc), file=sys.stderr)
        return 2
    print(json.dumps({"out_csv": str(args.out_csv), "row_count": len(rows)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
