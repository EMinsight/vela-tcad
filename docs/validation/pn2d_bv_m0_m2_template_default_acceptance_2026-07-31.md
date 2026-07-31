# PN2D BV M0/M2 prospective template-default acceptance

Date: 2026-07-31

Status: complete; candidate rejected; template default rolled back.

## Decision

The prospective machine decision is:

```text
pn2d_bv_template_default_not_accepted
```

`production_default_change_authorized=false` and `authorized_surface=none`.
The PN2D BV template default was therefore restored atomically to
`legacy_cell_reconstructed`.  The SG/Laux implementation, the named opt-in
profile, diagnostics, and tests remain available.  The global C++ defaults
were never changed.

The machine-readable evidence is under:

```text
build-release/pn2d-bv-template-default-prospective-20260731/
```

In particular:

- `acceptance/acceptance.json` is the unified decision;
- `acceptance/M0_acceptance.json` is the M0 level decision;
- `acceptance/M2_acceptance.json` is the requested M2 machine-readable
  acceptance;
- `M0/contract-domain-parity.json` and
  `M2/contract-domain-parity.json` contain the same-grid curve evidence;
- each level has independent `run-a` and `run-b` execution and state
  manifests.

Generated simulation outputs remain ignored and are not committed.

## Prospective contract

The contract was frozen before the new runs:

```text
docs/validation/contracts/pn2d_bv_m0_m2_template_default_acceptance_v1.json
SHA-256 53c4647dfe2916384ff861f2bb85e989eeabdd44b2330689754dadc38638f402
```

It requires both M0 and M2 to pass the same-grid Sentaurus-golden gates,
all off/IIC/on branches to complete their exact lattice, duplicate IV,
physics configuration, process-probe, and state evidence to be deterministic,
and global continuity closure to pass.  Cross-mesh convergence is explicitly
observation-only.

The evaluated SG/Laux candidate template hash was:

```text
b39120af630ba5baceb93516e63ce3603c51c005cbf6d37890809aec1f8f43f9
```

After the rejected decision, the committed template surface was rolled back
to:

```text
avalanche_current_support_profile=legacy_cell_reconstructed
current_approximation=cell_reconstructed
source_mapping_mode=triangle_gss_gradqf_truncated
cell_reconstructed_midpoint_density=gss_logistic
```

The post-rollback template hash is:

```text
c525e0b78d3142bdf37086dbb70086ead54a48f3a276000c25286294bc9b372b
```

## Same-grid results

| Metric | M0 | M2 | Frozen limit | Result |
|---|---:|---:|---:|---|
| effective curve median error (dex) | 0.002030 | 0.046355 | 0.05 | pass / pass |
| effective curve P95 error (dex) | 0.003708 | 0.075860 | 0.10 | pass / pass |
| effective curve maximum error (dex) | 0.003749 | 0.078827 | 0.15 | pass / pass |
| effective gain median error (dex) | 0.002186 | 0.041009 | 0.05 | pass / pass |
| effective gain maximum error (dex) | 0.003702 | 0.073480 | 0.10 | pass / pass |
| knee median error (dex) | 0.003045 | 0.057983 | 0.05 | pass / **fail** |
| knee maximum error (dex) | 0.003749 | 0.078827 | 0.10 | pass / pass |
| absolute `V_break` error (V) | 0.001 | 0.011 | 0.10 | pass / pass |
| absolute `V_slope` error (V) | 0.000819 | unavailable | 0.10 | pass / **fail** |
| adjacent-slope RMSE (dex/V) | 0.002659 | 0.055059 | 0.20 | pass / pass |

Both M2 curves have no `1 dex/V` slope crossing inside the frozen voltage
window, so `V_slope` is `null`.  The contract predeclared that estimator as
required; it was not waived after the result.

## Exact-lattice and determinism evidence

For both levels and both runs:

- avalanche-off completes 29/29 requested points through -20 V;
- avalanche-on completes 29/29 requested points through -20 V;
- off/on IV, physics configuration, process-probe, and available state hashes
  are identical between run A and run B.

IIC/postprocess is incomplete:

| Level | Requested point | Last failed/observed bias | Failure |
|---|---:|---:|---|
| M0 A/B | -6 V | -5.528062500160932 V | `carrier_row_convergence_line_search_rejected` |
| M2 A/B | -3 V | -2.104377500115871 V | `carrier_row_convergence_line_search_rejected` |

The failure bias and class repeat, but the partial IIC IV file hashes differ.
Consequently `complete_exact_lattice=false` and
`duplicate_determinism=false`.  IIC is postprocess-only and does not feed the
avalanche source into continuity; this failure does not by itself identify a
self-consistent SG/Laux source-feedback defect.  It still blocks the
predeclared diagnostic contract.

## Closure

The self-consistent avalanche-on branch passes the machine-readable global
continuity closure at every contract bias:

| Level | max electron ratio | max hole ratio | Limit |
|---|---:|---:|---:|
| M0 | 7.13696e-6 | 7.12422e-6 | 0.01 |
| M2 | 4.70138e-7 | 2.15322e-6 | 0.01 |

The evaluator also binds both process-probe files, checks their A/B hash
identity, required columns, and exact contract-bias coverage.

## Evidence-chain hardening

The final evaluator fails closed on:

- a parity curve that is not the recorded Vela branch output or the explicitly
  supplied Sentaurus aggregate;
- a render manifest that cannot regenerate the recorded base configuration;
- a branch config or IV whose actual hash differs from `execution.json`;
- a state manifest that is not passed, does not cover the complete execution
  lattice, or contains a snapshot hash mismatch;
- a process probe that lacks required fields or contract-bias coverage;
- non-finite machine metrics or non-standard JSON numeric values.

Failed executions still emit a typed partial state manifest so their available
hash evidence is retained.

## Verification

- focused Python contract/template suite: 28/28 passed;
- Release CTest after rollback: 506/506 passed;
- fresh default render: legacy profile and legacy three-field combination;
- explicit SG/Laux render: atomic three-field opt-in;
- explicit legacy render: atomic rollback;
- mixed and omitted three-field configurations: fail closed;
- global C++ omitted-field defaults: unchanged and covered by C++ assertions.

## Independent reviews

The second scientific review returns `REJECT_DEFAULT_CHANGE`.  It accepts the
strong M0 and useful M2 same-grid evidence but does not waive the M2 knee or
IIC failures.

The second code review returns `REJECT_DEFAULT_CHANGE`.  It requires the
template rollback and identified fail-closed evidence-chain checks.  Those
checks and their negative tests were completed before final verification.

## Next authorized work

SG/Laux remains opt-in.  A future default proposal requires a new,
pre-run-frozen contract.  Before that proposal:

1. diagnose the IIC/postprocess interaction with carrier-row convergence
   without changing avalanche physics;
2. predeclare how a slope crossing outside the voltage window is handled,
   or extend the prospective voltage window before running;
3. rerun independent M0/M2 default-candidate evidence under that new contract.

No threshold, bias point, or old acceptance artifact may be rewritten to
convert this rejected proposal into a pass.
