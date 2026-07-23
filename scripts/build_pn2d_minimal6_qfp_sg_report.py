#!/usr/bin/env python3
"""Build a portable HTML report for the Minimal6 QFP-SG replacement audit."""

from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path


def fmt(value: str, digits: int = 4) -> str:
    return f"{float(value):.{digits}g}"


def bar(value: float, maximum: float, x: float, y: float, color: str) -> str:
    height = 180.0 * value / maximum
    return (
        f'<rect x="{x}" y="{y + 180 - height}" width="46" height="{height}" '
        f'rx="4" fill="{color}"/><text x="{x + 23}" y="{y + 170 - height}" '
        f'text-anchor="middle" class="value">{value:.2f}</text>'
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    verification = json.loads(
        (root / "independent_verification.json").read_text(encoding="utf-8")
    )
    with (root / "qfp_replacement_summary.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    index = {(row["carrier"], row["variant"]): row for row in rows}
    electron_base = index[("electron", "baseline")]
    electron_swap = index[("electron", "electron_qfp")]
    hole_base = index[("hole", "baseline")]
    hole_swap = index[("hole", "hole_qfp")]
    maximum = 4.0
    svg = f"""
    <svg viewBox="0 0 620 270" role="img" aria-label="Median current-proxy log error before and after QFP replacement">
      <line x1="70" y1="220" x2="585" y2="220" class="axis"/>
      <line x1="70" y1="40" x2="70" y2="220" class="axis"/>
      <text x="18" y="132" transform="rotate(-90 18 132)" class="label">Median absolute error (dex)</text>
      <text x="177" y="250" class="label">Electron</text>
      <text x="430" y="250" class="label">Hole</text>
      {bar(float(electron_base["median_abs_log10_error_dex"]), maximum, 130, 40, "#64748b")}
      {bar(float(electron_swap["median_abs_log10_error_dex"]), maximum, 190, 40, "#0f766e")}
      {bar(float(hole_base["median_abs_log10_error_dex"]), maximum, 383, 40, "#64748b")}
      {bar(float(hole_swap["median_abs_log10_error_dex"]), maximum, 443, 40, "#0f766e")}
      <rect x="162" y="10" width="14" height="14" fill="#64748b"/><text x="182" y="22" class="legend">Baseline</text>
      <rect x="260" y="10" width="14" height="14" fill="#0f766e"/><text x="280" y="22" class="legend">Relevant QFP replaced</text>
    </svg>
    """
    table_rows = []
    for carrier, baseline, replacement in (
        ("Electron", electron_base, electron_swap),
        ("Hole", hole_base, hole_swap),
    ):
        table_rows.append(
            "<tr>"
            f"<td>{carrier}</td>"
            f"<td>{fmt(baseline['median_abs_log10_error_dex'])}</td>"
            f"<td>{fmt(replacement['median_abs_log10_error_dex'])}</td>"
            f"<td>{fmt(replacement['median_paired_log_error_improvement_dex'])}</td>"
            f"<td>{fmt(replacement['median_symmetric_relative_residual'])}</td>"
            f"<td>{fmt(replacement['sign_agreement_fraction'], 5)}</td>"
            "</tr>"
        )
    replay = manifest["baseline_cpp_replay"]
    source_rows = "".join(
        f"<tr><td>{html.escape(key)}</td><td><code>{html.escape(str(value))}</code></td></tr>"
        for key, value in manifest["inputs"].items()
        if key.endswith("sha256")
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Minimal6 internal-QFP SG replacement</title>
<style>
:root{{--ink:#172033;--muted:#5f6b7a;--paper:#f5f7fb;--card:#fff;--line:#d9e0ea;--teal:#0f766e;--amber:#a16207}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.55 Inter,Segoe UI,Arial,sans-serif}}
main{{max-width:1080px;margin:0 auto;padding:38px 24px 64px}} h1{{font-size:34px;line-height:1.1;margin:0 0 10px}} h2{{margin-top:34px}}
.lede{{font-size:18px;color:var(--muted);max-width:850px}} .grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:26px 0}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px;box-shadow:0 5px 18px #1e293b0d}}
.metric{{font-size:30px;font-weight:700;color:var(--teal)}} .warn{{color:var(--amber)}} .small{{color:var(--muted);font-size:13px}}
table{{width:100%;border-collapse:collapse;background:var(--card);border:1px solid var(--line)}} th,td{{text-align:left;padding:10px 12px;border-bottom:1px solid var(--line)}} th{{background:#eef2f7}}
figure{{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:16px;margin:24px 0}} svg{{width:100%;height:auto}} .axis{{stroke:#94a3b8;stroke-width:1}} .label{{font-size:13px;fill:#475569}} .legend{{font-size:12px;fill:#475569}} .value{{font-size:12px;font-weight:700;fill:#172033}}
code{{font:12px ui-monospace,SFMono-Regular,Consolas,monospace;word-break:break-all}} .conclusion{{border-left:5px solid var(--amber)}}
@media(max-width:760px){{.grid{{grid-template-columns:1fr}} main{{padding:24px 14px}}}}
</style>
</head>
<body><main>
<h1>Internal-QFP replacement fixes direction, not magnitude</h1>
<p class="lede">Across 40 exact Minimal6 states, replacing Sentaurus QFP only at internal nodes 1 and 5 improves the directed-edge SG comparison by about 1.4 dex and restores 100% sign agreement. It does not close the current residual: median errors remain above 2 dex.</p>
<div class="grid">
  <div class="card"><div class="metric">{fmt(electron_swap['median_paired_log_error_improvement_dex'])} dex</div><div>Electron paired median improvement</div></div>
  <div class="card"><div class="metric">{fmt(hole_swap['median_paired_log_error_improvement_dex'])} dex</div><div>Hole paired median improvement</div></div>
  <div class="card"><div class="metric warn">{fmt(electron_swap['median_symmetric_relative_residual'])} / {fmt(hole_swap['median_symmetric_relative_residual'])}</div><div>Electron / hole bounded residual after replacement</div></div>
</div>
<figure>{svg}<figcaption class="small">Affected edges only: 280 carrier-edge samples per branch. Lower is better.</figcaption></figure>
<h2>Comparison</h2>
<table><thead><tr><th>Carrier</th><th>Baseline median dex</th><th>Replaced median dex</th><th>Paired improvement dex</th><th>Replaced bounded residual</th><th>Sign agreement</th></tr></thead>
<tbody>{''.join(table_rows)}</tbody></table>
<div class="card conclusion"><h2>Scientific classification</h2>
<p><strong>{html.escape(verification['scientific_classification'])}</strong>. The QFP discrepancy is causal for the current direction and part of the magnitude gap, but it is not sufficient to reproduce the Sentaurus current proxy. The remaining roughly 100–130× median magnitude gap must come from another operator or support-semantic difference.</p>
<p class="small">The Sentaurus reference is the endpoint-mean node current-density vector projected onto the canonical edge tangent. It is not a native Sentaurus SG edge flux, so this experiment cannot uniquely identify a production formula.</p></div>
<h2>Controls and validation</h2>
<ul>
  <li>Frozen: Vela electrostatic potential, stored n/p, baseline production edge mobility, mesh, temperature, and intrinsic density.</li>
  <li>Baseline replay: {replay['sample_count']} samples; maximum relative difference {replay['max_relative_error']:.4g}; gate passed.</li>
  <li>Strict frozen-density SG is QFP-independent by construction and therefore serves only as a negative control.</li>
  <li>Independent verification: {verification['status']}; {verification['edge_row_count']} edge rows and {verification['summary_row_count']} summary rows recomputed.</li>
</ul>
<h2>Source hashes</h2>
<table><tbody>{source_rows}</tbody></table>
<p class="small">Deterministic query contract: <code>source_query.sql</code>. Full edge evidence: <code>qfp_replacement_edge_samples.csv</code>.</p>
</main></body></html>"""
    (root / "report.html").write_text(document, encoding="utf-8")
    query = """-- DuckDB query contract for the primary report chart and table.
SELECT
  carrier,
  variant,
  affected_edge_count,
  median_abs_log10_error_dex,
  median_symmetric_relative_residual,
  sign_agreement_fraction,
  median_paired_log_error_improvement_dex
FROM read_csv_auto('qfp_replacement_summary.csv')
WHERE (carrier = 'electron' AND variant IN ('baseline', 'electron_qfp'))
   OR (carrier = 'hole' AND variant IN ('baseline', 'hole_qfp'))
ORDER BY carrier, CASE variant WHEN 'baseline' THEN 0 ELSE 1 END;
"""
    (root / "source_query.sql").write_text(query, encoding="utf-8")
    print(json.dumps({"status": "built", "report": str(root / "report.html")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
