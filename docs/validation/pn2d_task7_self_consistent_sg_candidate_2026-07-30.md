# PN2D Task 7 self-consistent SG/Laux candidate qualification

Date: 2026-07-30

Outcome: `tradeoff_without_parity`.

## Candidate and scope

The controlled opt-in candidate replaces the archived triangle current/source
support with the complete directed SG current vector and matching element-box
geometry:

- `current_approximation=element_edge_sg_gss_laux`;
- `source_mapping_mode=element_vertex_box_measure`; and
- inactive compatibility field
  `cell_reconstructed_midpoint_density=bernoulli`.

The third assignment only removes the capability-guard-incompatible
`gss_logistic` selection owned by the retired triangle source. Mesh, doping,
mobility, Van Overstraeten coefficients, high-field drive, contacts,
continuation schedule, tolerances, and production defaults remain unchanged.

## Exact-lattice and determinism result

Run A completed avalanche-off, IIC/postprocess-only, and self-consistent
avalanche-on at all 29 exact points from `0` through `-20 V`. Run B
independently completed the same 29-point self-consistent avalanche-on branch.

The two avalanche-on IV files are byte-identical:

`de0d40cedfc7ca4f19f4c284e395877c345a4215e76bbbda9bd9d9a245f40de6`

The normalized process manifest contains 68,034 field records, 783 aggregate
records, and 87 exact-target Newton-attempt records. Its SHA-256 is:

`dd4da6d6ebca779010bca4fadc7377fd39ab9cb05bf9ecb0edfb56b32b96a43c`

The first launch returned Windows `0xC0000135` before writing an IV row because
the UCRT64 DLL directory was not on `PATH`. The recorded incomplete execution
was replaced by an in-place `--resume` run after applying the repository's
UCRT64 environment contract; no solver or convergence failure occurred.

## Curve and knee score

| Metric | Sealed baseline | SG/Laux candidate |
|---|---:|---:|
| knee log-current RMSE (dex) | `11.4007360393` | `0.0146387652` |
| RMSE improvement | n/a | `99.8716%` |
| `V_break` error (V) | `0.2320` | `0.0210` |
| `V_slope` error (V) | unavailable | `0.01643` |
| adjacent-slope RMSE (dex/V) | n/a | `0.04161` |
| knee median/max current error (dex) | n/a | `0.01708/0.02202` |
| maximum global-error worsening (dex) | n/a | `0.000797` |

The candidate also gives global-lattice median/P95 log-current errors of
`0.00478/0.01708 dex`. The global maximum remains `0.42072 dex` only at
`0 V`, where both currents are below the scientific numerical floor; it is
not a candidate regression.

## Process-chain and Task 6 feedback

The candidate process-chain analysis remains fully observed:

- 203/203 source-terminal closure rows pass;
- no stage observation is missing;
- first departure remains `state` on adjacent `-19.7/-19.8 V`; and
- the typed process outcome remains `density_qfp_feedback_cause`.

The state agreement improves materially:

| Internal metric | Sealed baseline | SG/Laux candidate |
|---|---:|---:|
| QFP RMSE (V) | `0.380441524` | `9.98598e-5` |
| density log-RMSE (dex) | `5.756812957` | `0.001673079` |
| Task 6 initial QFP RMSE (V) | `0.428636879` | `0.000114547` |

Duplicate Task 6 feedback outputs are byte-identical at both adjacent biases.
The feedback diagnosis changes from `continuation_only_cause` to
`electron_qfp_feedback_cause`: electron-QFP-only replacement is causal at both
biases, while the prior full-coupled reversal is no longer the dominant
signature. Task 6 still reports `task8_authorized=false`.

## Remaining failed gate

Ten of the eleven Task 7 authorization gates pass. The only failed gate is
`nonmonotonic_interval_removed`.

The candidate removes the baseline high-current reverse intervals
`-18→-18.5→-19→-19.25 V`, but introduces three low-current reverse intervals:

- `-3→-4 V`;
- `-4→-5 V`; and
- `-6→-7 V`.

All three occur near `3e-17 A/um`. This is a large scientific improvement over
the baseline high-field pathology, but the predeclared gate requires removal
rather than relocation and therefore fails closed. No current floor was added
after observing the result.

The machine-readable scorecard is:

`build-release/pn2d-task7-sg-laux-selfconsistent-scorecard-20260730/acceptance.json`

Its SHA-256 is:

`9f236ebed8bfdb370051fa7f54d2da3eabdd43743e6ec528f3bc3a541a5f8f17`

Task 8 and all production-default changes remain unauthorized.

## Verification

- Task 7 Python regression: 4/4 passed;
- Python syntax compilation: passed;
- element-edge SG/Laux test: 134 assertions in 8 cases passed;
- process-observability CTest: passed; and
- duplicate exact-lattice and Task 6 output hashes: passed.
