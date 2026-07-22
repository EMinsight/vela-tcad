#!/usr/bin/env python3
"""Build the deterministic PN2D Minimal6 physics inverse-audit package."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_diagnostics.inverse_inputs import load_input_bundle
from scripts.pn2d_minimal6_diagnostics.inverse_report import (
    build_analysis_artifacts, write_report_manifest,
)


def _sentaurus_version(*roots: Path) -> str:
    versions = set()
    for root in roots:
        manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
        provenance = manifest.get("provenance", {})
        value = provenance.get("sentaurus_version", manifest.get("sentaurus_version"))
        if isinstance(value, str) and value:
            versions.add(value)
    if len(versions) > 1:
        raise ValueError("input manifests declare inconsistent Sentaurus versions")
    return next(iter(versions)) if versions else "not_declared_by_input_contract"


def build_report_package(*, vela_root: str | Path, sentaurus_root: str | Path,
                         supplemental_sentaurus_root: str | Path,
                         out_dir: str | Path, phase_base: str) -> dict:
    roots = {
        "vela_root": str(Path(vela_root).resolve()),
        "sentaurus_root": str(Path(sentaurus_root).resolve()),
        "supplemental_sentaurus_root": str(Path(supplemental_sentaurus_root).resolve()),
    }
    # All input validation completes before the output directory is created.
    bundle = load_input_bundle(roots["vela_root"], roots["sentaurus_root"],
                               roots["supplemental_sentaurus_root"])
    version = _sentaurus_version(Path(roots["sentaurus_root"]),
                                  Path(roots["supplemental_sentaurus_root"]))
    build_analysis_artifacts(bundle, out_dir, phase_base=phase_base,
                             input_roots=roots, sentaurus_version=version)
    write_report_manifest(out_dir, bundle, roots)
    from scripts.verify_pn2d_minimal6_physics_inverse_audit import verify_report
    return verify_report(out_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-root", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--supplemental-sentaurus-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--phase-base", default="a5524cf")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = build_report_package(
        vela_root=args.vela_root, sentaurus_root=args.sentaurus_root,
        supplemental_sentaurus_root=args.supplemental_sentaurus_root,
        out_dir=args.out_dir, phase_base=args.phase_base,
    )
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
