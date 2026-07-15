# PN2D Minimal6 Task 1-8 handoff (2026-07-15)

## Objective and current conclusion

The active objective is to complete every remaining item in
`docs/superpowers/plans/2026-07-14-pn2d-minimal6-formula-difference-and-sweep.md`.
Tasks 1-7 have implementation and focused-regression evidence. Task 8 has real
fixed-state evidence and both real nonlinear sweep outcomes, but the final
comparison regeneration, validation-document update, full regression rerun, and
final validation commit remain to be completed.

All curves and sweep reports remain subject to this statement:

> minimal6 diagnostic sweep; not a physical BV curve

No production physics formula, production default, or solver gate was relaxed.

## Repository and branch state

- Repository: `D:\code-repo\vela-tcad`
- Implementation worktree:
  `D:\code-repo\vela-tcad\.worktrees\pn2d-minimal6-formula-sweep`
- Implementation branch: `codex-pn2d-minimal6-formula-sweep`
- Requested continuation branch: `codex-pn2d-minimal6-operator-audit`
- Baseline before this handoff commit: `752f91e Verify PN2D sweep comparison artifacts`
- Use `git log -1 --oneline` to obtain the handoff commit hash after checkout.

The requested continuation branch is attached to the main checkout. The handoff
commit is created in the implementation worktree first and then the main checkout
branch is fast-forwarded to that commit.

## Task status

| Task | Status | Evidence / limitation |
| --- | --- | --- |
| 1 | Complete | Six exact extended-field states, replay provenance, and a separately generated 355-member recovery seal. |
| 2 | Complete | Ledger contracts, units, deterministic serialization, and both public report schemas. |
| 3 | Complete | Independent geometry, physics, alpha, intrinsic-density, integration, and support-conversion routines. |
| 4 | Complete with explicit stop conclusion | Fixed-state report validates but returns `insufficient_data`; no causal ranking was fabricated. |
| 5 | Complete | Five PNG/PDF figure pairs and reviewed figure manifest. |
| 6 | Complete | Vela and Sentaurus drivers, real accepted endpoints, and preserved first failures. |
| 7 | Implemented; regenerate after this commit | Existing comparison was generated before the Vela source-family correction and must be replaced. |
| 8 | In progress | Real sweeps are complete; final comparison, docs, full verification, and final validation commit remain. |

## Authoritative fixed-state and recovery evidence

Authoritative state root:

`D:\code-repo\vela-tcad\build-release\reference_tcad\pn2d_sentaurus2018_minimal6\state_exports\minimal6_states_live_20260713_v2`

- Run ID: `minimal6_states_live_20260713_v2`
- Manifest SHA-256:
  `b44ad95d5df6d57383ba3d5b292818568e358d67f0fc0424ee72f95b673e8aaa`
- Matrix: sketch/mirror x 0/-12/-19 V
- Recovery seal:
  `D:\code-repo\vela-tcad\build-release\reference_tcad\pn2d_sentaurus2018_minimal6\recovery_validation\minimal6_states_live_20260713_v2\recovery_validation.json`
- Seal member count: 355
- Seal SHA-256:
  `9466ee2db317fda4707254403fa33c2e1ebee666f8fdc72b776fa4a8ab689ec3`
- Seal verification result: zero member-hash mismatches

Fixed-state formula report:

`D:\code-repo\vela-tcad\.worktrees\pn2d-minimal6-formula-sweep\build-release\pn2d-minimal6-formula-diff-task8-20260715`

It contains six records, exact 36/54/24 node/edge/triangle identities, five
reviewed figures, and the explicit `insufficient_data` conclusion.

## Real nonlinear sweep outcomes

### Vela

Manifest:

`D:\code-repo\vela-tcad\.worktrees\pn2d-minimal6-formula-sweep\build-release\pn2d-minimal6-vela-task8-final-r2-20260715\sweep_manifest.json`

- Accepted exact checkpoints: sketch -1 V and mirror -1 V.
- At -1 V, the production SG total source and the sum of per-edge reconstructed
  sources close exactly in the recorded output.
- First retained failure for both topologies: transition -1 V to -2 V.
- Failure class: `nonfinite_residual`; solver exit code 1.
- Physics and solver settings were not changed or retried.

### Sentaurus

Remote computation completed for sketch and mirror through -20 V. The remote
controller accepted every target from -2 V through -20 V.

Downloaded raw root:

`D:\code-repo\vela-tcad\.worktrees\pn2d-minimal6-formula-sweep\build-release\pn2d-minimal6-sentaurus-remote-task8-v2-20260715`

Imported sweep root:

`D:\code-repo\vela-tcad\.worktrees\pn2d-minimal6-formula-sweep\build-release\pn2d-minimal6-sentaurus-task8-final-20260715`

- Accepted exact checkpoints: 40 (sketch/mirror x -1..-20 V).
- First failure: none.
- Two attempted 0 V rows in this generated manifest were rejected before the
  field-selector fix because the extended TDR contains scalar and vector
  `ElectricField` datasets. Do not cite those two rejected rows as solver
  failures. The code in this handoff commit fixes the selector by requiring one
  full `(name, region, components, unit)` contract match.

The 40 accepted -1..-20 V rows are valid and can be used directly for the final
comparison. A clean 0 V sweep package may be regenerated later if desired.

## Code corrections in the handoff commit

1. Vela source families are no longer rejected for lacking a Sentaurus-native
   field. The driver records the production Vela SG total as the Vela native
   source and independently sums `edge_source_integral` as the reconstructed
   source. It requires closure at relative tolerance `1e-12`.
2. Sentaurus field selection now filters by normalized name, region, component
   count, and unit before enforcing uniqueness. A scalar/vector namesake no
   longer causes a false duplicate rejection, while two contract-compatible
   duplicates still fail closed.
3. Focused diagnostic-sweep regressions pass: 21 tests, zero failures.

## Resume commands at home

Open PowerShell and run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
Set-Location "D:\code-repo\vela-tcad"
git branch --show-current
git log -1 --oneline
git status --short
```

The branch should be `codex-pn2d-minimal6-operator-audit` and the tracked tree
