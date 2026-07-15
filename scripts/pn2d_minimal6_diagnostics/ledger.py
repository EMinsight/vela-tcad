import csv
import json
from pathlib import Path
from .contracts import QuantityRecord

_FIELDS = ["run_id", "topology", "bias_V", "carrier", "support_kind", "support_id", "quantity", "source", "formula_version", "value", "unit", "sign_convention", "raw_source_path", "raw_source_sha256", "geometric_zero_reason"]

class DiagnosticLedger:
    def __init__(self): self._records = {}
    def add(self, record: QuantityRecord):
        key = (record.state.run_id, record.state.topology, record.state.bias_V, record.carrier, record.support_kind, record.support_id, record.quantity, record.source, record.formula_version)
        if key in self._records: raise ValueError("duplicate diagnostic ledger key")
        self._records[key] = record
    def records(self): return [self._records[key] for key in sorted(self._records, key=str)]
    @staticmethod
    def _row(r: QuantityRecord) -> dict[str, object]:
        return {"run_id": r.state.run_id, "topology": r.state.topology, "bias_V": r.state.bias_V, "carrier": r.carrier, "support_kind": r.support_kind, "support_id": r.support_id, "quantity": r.quantity, "source": r.source, "formula_version": r.formula_version, "value": r.value, "unit": r.unit, "sign_convention": r.sign_convention, "raw_source_path": r.raw_source_path, "raw_source_sha256": r.raw_source_sha256, "geometric_zero_reason": r.geometric_zero_reason}
    def write_json(self, path: Path) -> None:
        Path(path).write_text(json.dumps([self._row(r) for r in self.records()], sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    def write_csv(self, path: Path) -> None:
        with Path(path).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=_FIELDS)
            writer.writeheader(); writer.writerows(self._row(r) for r in self.records())