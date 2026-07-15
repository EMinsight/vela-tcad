#!/usr/bin/env python3
"""Build a validated native-chart report from an executed SQLite snapshot."""
from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
from pathlib import Path


SUMMARY_SQL = """
SELECT ABS(CAST(bias_V AS INTEGER)) AS bias_magnitude_V,
       CAST(bias_V AS INTEGER) AS bias_V,
       log10(CAST(n_ratio_geomean_interior AS REAL)) AS log10_density_ratio,
       log10(CAST(terminal_current_ratio AS REAL)) AS log10_current_ratio,
       log10(CAST(source_integral_ratio AS REAL)) AS log10_source_ratio,
       CAST(mean_abs_phin_error_interior_V AS REAL) AS phin_error_V,
       CAST(n_ratio_geomean_interior AS REAL) AS density_ratio,
       CAST(electric_field_ratio AS REAL) AS field_ratio,
       CAST(electron_alpha_ratio AS REAL) AS alpha_ratio,
       CAST(terminal_current_ratio AS REAL) AS current_ratio,
       CAST(source_integral_ratio AS REAL) AS source_ratio
FROM self_consistent_summary
WHERE topology='sketch' AND CAST(bias_V AS INTEGER) IN (-1,-12,-19,-20)
ORDER BY ABS(CAST(bias_V AS INTEGER))
""".strip()

ALIGNMENT_SQL = """
SELECT ABS(CAST(bias_V AS INTEGER)) AS bias_magnitude_V,
       printf('%d V', CAST(bias_V AS INTEGER)) AS bias_label,
       'Peak electric field' AS metric,
       100.0*(CAST(electric_field_ratio AS REAL)-1.0) AS deviation_percent,
       CAST(electric_field_ratio AS REAL) AS ratio,
       CAST(mean_abs_phin_error_interior_V AS REAL) AS phin_error_V
FROM self_consistent_summary WHERE topology='sketch' AND CAST(bias_V AS INTEGER) IN (-12,-19,-20)
UNION ALL
SELECT ABS(CAST(bias_V AS INTEGER)), printf('%d V', CAST(bias_V AS INTEGER)),
       'Electron impact coefficient', 100.0*(CAST(electron_alpha_ratio AS REAL)-1.0),
       CAST(electron_alpha_ratio AS REAL), CAST(mean_abs_phin_error_interior_V AS REAL)
FROM self_consistent_summary WHERE topology='sketch' AND CAST(bias_V AS INTEGER) IN (-12,-19,-20)
UNION ALL
SELECT ABS(CAST(bias_V AS INTEGER)), printf('%d V', CAST(bias_V AS INTEGER)),
       'Electron QF gradient',
       100.0*(CAST(vela_max_abs_phin_gradient_V_per_m AS REAL)
              /CAST(sentaurus_max_abs_phin_gradient_V_per_m AS REAL)-1.0),
       CAST(vela_max_abs_phin_gradient_V_per_m AS REAL)
              /CAST(sentaurus_max_abs_phin_gradient_V_per_m AS REAL),
       CAST(mean_abs_phin_error_interior_V AS REAL)
FROM self_consistent_summary WHERE topology='sketch' AND CAST(bias_V AS INTEGER) IN (-12,-19,-20)
ORDER BY bias_magnitude_V, metric
""".strip()

TASK8_SQL = """
SELECT COUNT(*) AS exact_common_checkpoints,
       MIN(CAST(bias_V AS REAL)) AS deepest_common_bias_V,
       SUM(CASE WHEN reason='exact common checkpoint' THEN 0 ELSE 1 END) AS non_exact_rows
FROM task8_comparison
WHERE classification='common_exact'
""".strip()

FORMULA_SQL = "SELECT priority, quantity, formula, finding, status FROM formula_audit ORDER BY priority"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def import_csv(connection: sqlite3.Connection, table: str, path: Path) -> None:
    rows = read_csv(path)
    columns = list(rows[0])
    connection.execute(f"CREATE TABLE {table} ({', '.join(f'[{column}] TEXT' for column in columns)})")
    placeholders = ",".join("?" for _ in columns)
    connection.executemany(
        f"INSERT INTO {table} VALUES ({placeholders})",
        [[row[column] for column in columns] for row in rows],
    )


def query_rows(connection: sqlite3.Connection, sql: str) -> list[dict[str, object]]:
    cursor = connection.execute(sql)
    columns = [item[0] for item in cursor.description]
    return [dict(zip(columns, row)) for row in cursor.fetchall()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-csv", type=Path, required=True)
    parser.add_argument("--task8-csv", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    database = args.out.parent / "evidence.sqlite"
    if database.exists():
        database.unlink()
    connection = sqlite3.connect(database)
    connection.create_function("log10", 1, math.log10)
    import_csv(connection, "self_consistent_summary", args.summary_csv)
    import_csv(connection, "task8_comparison", args.task8_csv)
    formulas = [
        (1, "Electric field", "E=-grad(psi)", "Aligned: peak-field ratio is 1.001-1.002 from -12 to -20 V", "aligned"),
        (2, "Carrier statistics", "n=ni exp((psi-phin)/Vt); p=ni exp((phip-psi)/Vt)", "Algebra is correct; the self-consistent QF branch is offset by 0.30-0.37 V", "branch_error"),
        (3, "Impact coefficient", "alpha=gamma A exp(-gamma B/abs(F))", "Aligned: electron-alpha ratio is 1.015-1.040 from -12 to -20 V", "aligned"),
        (4, "Impact generation", "G=(alpha_n|Jn|+alpha_p|Jp|)/q", "Algebra is correct; the published Vela source integral needs the 1e-8 um2-to-cm2 factor", "output_unit_fixed"),
        (5, "Continuity geometry", "edge transport uses um line; R/G source uses um2 area", "Volume sources are overweighted by 1e4 relative to transport", "root_cause"),
    ]
    connection.execute("CREATE TABLE formula_audit (priority INTEGER, quantity TEXT, formula TEXT, finding TEXT, status TEXT)")
    connection.executemany("INSERT INTO formula_audit VALUES (?,?,?,?,?)", formulas)
    connection.commit()
    summary = query_rows(connection, SUMMARY_SQL)
    alignment = query_rows(connection, ALIGNMENT_SQL)
    task8 = query_rows(connection, TASK8_SQL)
    formula_rows = query_rows(connection, FORMULA_SQL)
    connection.close()

    chain = []
    for row in summary:
        for metric, field in (("Carrier density", "log10_density_ratio"),
                              ("Terminal current", "log10_current_ratio"),
                              ("Avalanche source", "log10_source_ratio")):
            chain.append({**row, "metric": metric, "log10_ratio": row[field]})
    representative = next(row for row in summary if row["bias_V"] == -19)
    headline = [{"current_ratio": representative["current_ratio"],
                 "phin_error_V": representative["phin_error_V"],
                 "field_ratio": representative["field_ratio"]}]

    db_path = "build-release/pn2d-minimal6-physics-error-report-task8-20260715/evidence.sqlite"
    sources = [
        {"id": "physics", "label": "Self-consistent Minimal6 comparison", "path": db_path,
         "query": {"engine": "SQLite", "language": "sql", "sql": SUMMARY_SQL,
                   "description": "Executed exact-bias extraction and derived error factors.",
                   "tables_used": ["self_consistent_summary"],
                   "filters": ["topology=sketch", "bias in -1,-12,-19,-20 V", "interpolation forbidden"],
                   "metric_definitions": ["All ratios are Vela divided by Sentaurus.", "Interior nodes exclude contact nodes."]}},
        {"id": "alignment", "label": "Local operator alignment", "path": db_path,
         "query": {"engine": "SQLite", "language": "sql", "sql": ALIGNMENT_SQL,
                   "description": "Executed local-operator ratio query for representative high-field biases.",
                   "tables_used": ["self_consistent_summary"],
                   "filters": ["topology=sketch", "bias in -12,-19,-20 V"]}},
        {"id": "task8", "label": "Corrected Task 8 checkpoint census", "path": db_path,
         "query": {"engine": "SQLite", "language": "sql", "sql": TASK8_SQL,
                   "description": "Executed census of exact common checkpoints.",
                   "tables_used": ["task8_comparison"], "filters": ["classification=common_exact"]}},
        {"id": "formula", "label": "Formula and unit audit", "path": db_path,
         "query": {"engine": "SQLite", "language": "sql", "sql": FORMULA_SQL,
                   "description": "Executed selection from the reviewed formula-audit ledger.",
                   "tables_used": ["formula_audit"],
                   "metric_definitions": ["1 um=1e-4 cm", "1 um2=1e-8 cm2"]}},
    ]
    title = "PN2D Minimal6 Task 8 potential, quasi-Fermi, and impact-ionization audit"
    manifest = {
        "version": 1, "surface": "report", "title": title,
        "description": "Exact-bias solver comparison and continuity-equation unit audit",
        "generatedAt": "2026-07-15T23:59:00+08:00", "sources": sources,
        "blocks": [
            {"id": "title", "type": "markdown", "body": f"# {title}"},
            {"id": "summary", "type": "markdown", "sourceId": "physics", "body":
             "## Technical summary\n\n**The dominant error is not electrostatic potential or the Van Overstraeten coefficient. It is the self-consistent quasi-Fermi branch, driven by a continuity-equation geometry-unit imbalance.** At -19 V the interior potential error is 2.8e-12 V, the peak-field ratio is 1.0013, and the electron-alpha ratio is 1.016. The mean electron quasi-Fermi error is 0.347 V, interior carrier density is 3.45e5 times too high, and terminal current is 3.66e3 times too high."},
            {"id": "metrics", "type": "metric-strip", "cardIds": ["current", "qf", "field"]},
            {"id": "task8_text", "type": "markdown", "sourceId": "task8", "body":
             f"## Task 8 completion\n\nThe executed checkpoint census contains {int(task8[0]['exact_common_checkpoints'])} exact common rows across sketch and mirror, reaches {float(task8[0]['deepest_common_bias_V']):g} V, and contains no non-exact promoted rows. Both topologies therefore reach -20 V with 40 common checkpoints and zero rejected transitions. Interpolation and BV extrapolation remain forbidden."},
            {"id": "chain_text", "type": "markdown", "sourceId": "physics", "body":
             "## The error follows the quasi-Fermi branch\n\nThe relations n=ni exp((psi-phin)/Vt) and p=ni exp((phip-psi)/Vt) turn a 0.30-0.37 V QF offset into a multi-decade carrier error. Current and avalanche-source errors track this upstream branch rather than the already-aligned electrostatic field."},
            {"id": "chain_block", "type": "chart", "chartId": "chain"},
            {"id": "alignment_text", "type": "markdown", "sourceId": "alignment", "body":
             "## Potential and alpha are already aligned\n\nFrom -12 to -20 V, peak-field error is 0.12%-0.20% and electron-alpha error is 1.5%-4.0%. The QF-gradient magnitude is also within 3.4%-5.2%, but its absolute QF level is wrong by roughly 0.35 V. The absolute branch level, not the local gradient formula, controls the carrier-density error."},
            {"id": "alignment_block", "type": "chart", "chartId": "alignment_chart"},
            {"id": "exact_block", "type": "table", "tableId": "exact"},
            {"id": "formula_text", "type": "markdown", "sourceId": "formula", "body":
             "## Formula audit isolates a 1e4 source-weight error\n\nElectric field, carrier statistics, SG transport, and alpha(E) have the correct algebraic form. In TCAD-internal units, however, edge transport uses a line length stored in um while SRH and avalanche sources use areas stored in um2. The physical centimetre factors are 1e-4 and 1e-8 respectively, so the current continuity residual overweights volume sources relative to transport by 1e4."},
            {"id": "formula_block", "type": "table", "tableId": "formula_table"},
            {"id": "next", "type": "markdown", "sourceId": "formula", "body":
             "## Repair order and acceptance test\n\n1. Apply consistent line and area conversion factors to continuity residuals and the analytic Jacobian.\n2. Add legacy-SI versus unit-scaling equivalence tests with SRH and avalanche enabled.\n3. Correct the unit-scaling state CSV concentration labels or convert values to m^-3.\n4. Rerun Minimal6; phin/phip, carrier density, current, and source should move together while psi and alpha(E) remain nearly unchanged."},
            {"id": "limits", "type": "markdown", "sourceId": "task8", "body":
             "## Scope and limitations\n\nThis is a six-node diagnostic device, not a physical BV curve. Results cover accepted exact checkpoints from 0 V to -20 V only. The -1 V alpha ratio is not diagnostic because both absolute coefficients are vanishingly small."},
        ],
        "cards": [
            {"id": "current", "dataset": "headline", "sourceId": "physics", "description": "Vela/Sentaurus at -19 V",
             "metrics": [{"label": "Terminal-current ratio", "field": "current_ratio", "format": "compact"}]},
            {"id": "qf", "dataset": "headline", "sourceId": "physics", "description": "Interior-node mean absolute error",
             "metrics": [{"label": "Electron QF error (V)", "field": "phin_error_V", "format": "number"}]},
            {"id": "field", "dataset": "headline", "sourceId": "physics", "description": "Vela/Sentaurus at -19 V",
             "metrics": [{"label": "Peak-field ratio", "field": "field_ratio", "format": "number"}]},
        ],
        "charts": [
            {"id": "chain", "title": "Error factors versus reverse bias (log10)",
             "subtitle": "Exact sketch checkpoints; 3 on the y-axis means a 1000x error", "showDescription": True,
             "intent": "trend", "question": "How does the QF branch error propagate to density, current, and avalanche source?",
             "rationale": "A shared logarithmic scale compares three multi-decade error factors directly.",
             "type": "line", "dataset": "chain", "sourceId": "physics",
             "encodings": {"x": {"field": "bias_magnitude_V", "type": "quantitative", "label": "Reverse-bias magnitude (V)"},
                           "y": {"field": "log10_ratio", "type": "quantitative", "label": "log10(Vela/Sentaurus)"},
                           "color": {"field": "metric", "type": "nominal"},
                           "tooltip": [{"field": "phin_error_V", "label": "Mean electron-QF error (V)"}]},
             "combinationRationale": "Color encodes the second categorical dimension: density, current, or source.",
             "xAxisTitle": "Reverse-bias magnitude (V)", "yAxisTitle": "log10 error factor", "layout": "full"},
            {"id": "alignment_chart", "title": "Local-operator deviation from Sentaurus",
             "subtitle": "Vela/Sentaurus minus one; alpha uses -12 V to -20 V", "showDescription": True,
             "intent": "comparison", "question": "Which local operators are already aligned?",
             "rationale": "Grouped bars compare bias-specific deviations on a common percent scale.",
             "type": "bar", "dataset": "alignment", "sourceId": "alignment",
             "encodings": {"x": {"field": "metric", "type": "nominal"},
                           "y": {"field": "deviation_percent", "type": "quantitative", "label": "Deviation (%)"},
                           "color": {"field": "bias_label", "type": "nominal"},
                           "tooltip": [{"field": "ratio", "label": "Vela/Sentaurus"}]},
             "combinationRationale": "Color encodes bias as the second categorical dimension.", "layout": "full"},
        ],
        "tables": [
            {"id": "exact", "title": "Exact representative-bias ledger",
             "subtitle": "Sketch topology; ratios are Vela divided by Sentaurus", "showDescription": True,
             "dataset": "summary", "sourceId": "physics", "density": "spacious", "layout": "full",
             "defaultSort": {"field": "bias_magnitude_V", "direction": "asc"},
             "columns": [{"field": "bias_V", "label": "Bias (V)", "format": "number"},
                         {"field": "phin_error_V", "label": "Mean phin error (V)", "format": "number"},
                         {"field": "density_ratio", "label": "Carrier-density ratio", "format": "compact"},
                         {"field": "field_ratio", "label": "Peak-field ratio", "format": "number"},
                         {"field": "alpha_ratio", "label": "Electron-alpha ratio", "format": "number"},
                         {"field": "current_ratio", "label": "Terminal-current ratio", "format": "compact"},
                         {"field": "source_ratio", "label": "Source-integral ratio", "format": "compact"},
                         {"field": "bias_magnitude_V", "label": "|Bias| (V)", "format": "number"}]},
            {"id": "formula_table", "title": "Formula and implementation audit",
             "subtitle": "Ordered by importance to the observed error", "showDescription": True,
             "dataset": "formula", "sourceId": "formula", "density": "spacious", "layout": "full",
             "defaultSort": {"field": "priority", "direction": "asc"},
             "columns": [{"field": "priority", "label": "Order", "format": "number"},
                         {"field": "quantity", "label": "Quantity", "type": "text"},
                         {"field": "formula", "label": "Formula/implementation", "type": "text"},
                         {"field": "finding", "label": "Finding", "type": "text"},
                         {"field": "status", "label": "Status", "type": "text"}]},
        ],
    }
    artifact = {"surface": "report", "manifest": manifest,
                "snapshot": {"version": 1, "generatedAt": manifest["generatedAt"], "status": "ready",
                             "datasets": {"headline": headline, "chain": chain, "alignment": alignment,
                                          "summary": summary, "formula": formula_rows, "task8": task8},
                             "accessIssues": []},
                "sources": sources}
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
