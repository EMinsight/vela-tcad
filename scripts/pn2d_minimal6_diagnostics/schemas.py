import math

DISCLAIMER = "minimal6 diagnostic sweep; not a physical BV curve"
_FORMULA_FIELDS = {"input_provenance", "audit_provenance", "state_matrix", "row_counts", "waterfall_paths", "interactions", "dominance_rules", "sentaurus_internal_semantics_residual", "vela_parameter_agreement", "artifact_hashes", "records"}
_BV_FIELDS = {"solver_configurations", "accepted_transitions", "failed_transitions", "checkpoints", "terminal_currents", "maximum_fields", "source_integrals", "convergence_metadata", "curve_artifact_hashes", "records"}
_BV_COMPARISON_FIELDS = {"deepest_common_bias_V", "missing_tails", "topology_sensitivity", "fixed_state_recheck", "failure_transitions", "artifact_hashes", "input_artifacts", "closure"}
def _finite(value):
    if isinstance(value, float) and not math.isfinite(value): raise ValueError("non-finite values are forbidden")
    if isinstance(value, dict):
        for item in value.values(): _finite(item)
    elif isinstance(value, list):
        for item in value: _finite(item)
def _validate(report, schema, fields):
    if report.get("schema") != schema: raise ValueError("invalid report schema")
    if report.get("diagnostic_disclaimer") != DISCLAIMER: raise ValueError("missing diagnostic disclaimer")
    missing = sorted(fields - set(report))
    if missing: raise ValueError(f"missing required report fields: {missing}")
    _finite(report)
def validate_formula_difference_v1(report:dict)->None:
    _validate(report, "vela.pn2d_minimal6_formula_difference.v1", _FORMULA_FIELDS)
def validate_bv_comparison_v1(report:dict)->None:
    _validate(report, "vela.pn2d_minimal6_bv_comparison.v1", _BV_FIELDS)
    missing = sorted(_BV_COMPARISON_FIELDS - set(report))
    if missing: raise ValueError(f"missing required comparison fields: {missing}")
    if report.get("interpolation") != "forbidden": raise ValueError("comparison must forbid interpolation")
    if report.get("closure", {}).get("status") != "closed": raise ValueError("comparison closure is not closed")