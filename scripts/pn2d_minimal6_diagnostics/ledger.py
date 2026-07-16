"""Unique, deterministically serialized Minimal6 quantity ledger."""

import csv
import json
import math
from pathlib import Path
import re

from .contracts import QuantityRecord


_FIELDS = [
    "run_id", "topology", "bias_V", "carrier", "support_kind", "support_id",
    "quantity", "source", "formula_version", "value", "unit", "sign_convention",
    "raw_source_path", "raw_source_sha256", "geometric_zero_reason",
]


def _canonical_support_id(value: str) -> tuple[tuple[int, object], ...]:
    return tuple(
        (0, int(part)) if part.isdigit() else (1, part.casefold())
        for part in re.split(r"(\d+)", value)
        if part
    )


def _sort_key(record: QuantityRecord) -> tuple[object, ...]:
    return (
        record.state.topology,
        float(record.state.bias_V),
        record.carrier,
        record.support_kind.value,
        _canonical_support_id(record.support_id),
        record.quantity,
        record.source.value,
        record.formula_version,
        record.state.run_id,
    )


class DiagnosticLedger:
    def __init__(self) -> None:
        self._records: dict[tuple[object, ...], QuantityRecord] = {}

    def add(self, record: QuantityRecord) -> None:
        if not isinstance(record, QuantityRecord):
            raise TypeError("diagnostic ledger accepts only QuantityRecord values")
        if record.key in self._records:
            raise ValueError("duplicate diagnostic ledger key")
        self._records[record.key] = record

    def records(self) -> list[QuantityRecord]:
        return sorted(self._records.values(), key=_sort_key)

    @staticmethod
    def _row(record: QuantityRecord) -> dict[str, object]:
        if not math.isfinite(float(record.state.bias_V)) or not math.isfinite(float(record.value)):
            raise ValueError("diagnostic ledger cannot serialize non-finite values")
        return {
            "run_id": record.state.run_id,
            "topology": record.state.topology,
            "bias_V": record.state.bias_V,
            "carrier": record.carrier,
            "support_kind": record.support_kind.value,
            "support_id": record.support_id,
            "quantity": record.quantity,
            "source": record.source.value,
            "formula_version": record.formula_version,
            "value": record.value,
            "unit": record.unit,
            "sign_convention": record.sign_convention,
            "raw_source_path": record.raw_source_path,
            "raw_source_sha256": record.raw_source_sha256,
            "geometric_zero_reason": record.geometric_zero_reason,
        }

    def write_json(self, path: Path) -> None:
        rows = [self._row(record) for record in self.records()]
        Path(path).write_text(
            json.dumps(rows, sort_keys=True, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )

    def write_csv(self, path: Path) -> None:
        rows = [self._row(record) for record in self.records()]
        with Path(path).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
