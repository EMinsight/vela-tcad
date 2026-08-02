# PN2D Source-Volume-Factor Adaptive Bisection Design

## Goal

Narrow the stable/unstable boundary of the diagnostic
`impact_ionization.source_volume_factor` axis and determine whether the highest
stable scalar factor can recover the Sentaurus breakdown-knee marker. This is a
controlled physical follow-up experiment, not a production-default change.

## Fixed Baseline

Use the existing `qflim0p05_high_source_scan` configuration and artifacts as
the immutable baseline. Every candidate keeps these settings identical:

- PN2D Sentaurus-2018 mesh and doping inputs;
- reverse sweep from `0 V` to `-20 V`;
- `quasi_fermi_update_limit_V=0.05`;
- secant predictor for `psi`, `phin`, and `phip`;
- terminal-current consistency branch gate;
- carrier-density branch gate with electron p95 limit `2.0 dex`;
- `driving_force=quasi_fermi_gradient`;
- the same Newton, line-search, adaptive-step, mobility, recombination, and
  impact-ionization settings.

Only `impact_ionization.source_volume_factor` and candidate-specific output
paths may differ.

The initial bracket is:

```text
stable lower bound:   0.921875
unstable upper bound: 0.9375
```

## Adaptive Algorithm

Run candidates sequentially. For each iteration, choose the exact arithmetic
midpoint of the current stable and unstable bounds.

A candidate is stable only if it reaches an accepted converged row at exactly
`-20 V`. A run that exits nonzero, stops before `-20 V`, produces an incomplete
CSV, or reports a branch/Newton failure is unstable. An unstable run must not
be reclassified by relaxing continuation or branch-guard settings.

After a stable result, replace the lower bound with the candidate. After an
unstable result, replace the upper bound. Stop when the bracket width is at
most `0.001953125` (one eighth of the initial width). This requires at most
three new solver runs.

## Candidate Evidence

Each candidate gets a separate ignored artifact directory containing:

- the complete generated simulation JSON;
- solver CSV and last-state CSV;
- captured command, exit code, stdout, and stderr;
- a machine-readable run summary;
- knee-shape comparison against the existing Sentaurus reference;
- Newton failure diagnostics when the run is unstable.

The summary records deepest accepted bias, converged row count, terminal
current, maximum accepted p95 carrier-density jump, failure bias/reason,
maximum absolute log10 current error over `-20..-10 V`, and the first 1 V
growth-ratio crossings above `1.5` and `2.0`.

## Failure Classification

For every unstable candidate, inspect the first failed transition and classify
it using existing diagnostics:

- `continuation_branch_selection` when predictor choice, branch acceptance,
  retry trajectory, or carrier-density jump changes discontinuously;
- `newton_step_pathology` when the residual is small but the proposed Newton
  step is disproportionately large or line search rejects all reductions;
- `source_volume_magnitude` when the branch trajectory remains shared and only
  current/source magnitude changes smoothly;
- `insufficient_diagnostics` when required failure evidence is absent.

Carrier positivity, finiteness, block residuals, Newton step norm, line-search
attempts, and the largest Poisson-residual nodes are retained in the failure
record when available.

## Final Decision Rules

The experiment reports the highest stable factor and lowest unstable factor at
the target bracket width.

- If a stable candidate produces a `>2.0` knee marker with improved Sentaurus
  current error, scalar source ownership remains a viable diagnostic axis.
- If the highest stable candidate still lacks the `>2.0` marker, or improvement
  is only smooth magnitude rescaling, scalar source ownership is insufficient
  to recover the knee on its own.
- If instability is classified as continuation or Newton-step behavior, the
  next experiment must target that mechanism rather than relaxing the physical
  acceptance gates.

No result changes the default `source_volume_factor`, establishes a physical
BV curve, or validates a production breakdown voltage.

## Verification And Documentation

Before publication:

- validate every candidate config differs from the baseline only in the factor
  and output paths;
- verify the bracket update history is internally consistent;
- rerun the existing knee-shape and factor-branch diagnostic tools;
- run focused regression for affected analysis scripts if code changes are
  needed;
- update `docs/validation/pn2d_bv_validation.md` with exact evidence and
  limitations;
- keep generated solver outputs out of Git.
