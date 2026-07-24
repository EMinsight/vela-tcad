# PN2D Minimal6 Task 10 scientific review request

Date: 2026-07-24

Review status: requested; deterministic independent verifiers pass.

## Requested verdict

Review whether the following conclusions are supported without fitted
parameters:

1. Sentaurus native low-field mobility is bias invariant on the six
   `mirror/sketch x -1/-10/-20 V` controls.
2. Native element electric field is the strongest tested high-field mobility
   drive, but the config-only self-consistent candidate improves QFP only
   modestly and misses all parity targets.
3. The current production impact-source gap is dominated by the triangle
   current proxy/support, not by the Van Overstraeten alpha coefficient.
4. No production mobility, SG, impact, SRH, Poisson, or continuation formula
   change is justified by the current evidence.

## Scientific checks already closed

| Check | Result |
|---|---:|
| exact physical states | 40, no interpolation |
| imported-state density maximum error | `4.426181e-6 dex` |
| electric-field mobility replay median, electron/hole | `0.000833/0.000443 dex` |
| box-current replay sign agreement | 100 percent |
| Sentaurus nodal alpha-current source closure median | `5.25100e-5 dex` |
| Sentaurus node-to-triangle projection median | `2.58704e-10 dex` |
| triangle local-volume geometry closure | `0.0` relative |
| triangle source-to-node mapping closure | `1.50363e-16` relative |
| source-unit isolated SRH/impact ratio error | `3.30872e-24` |
| analytic/finite-difference Jacobian maximum | `2.98950e-9` |

## Required reviewer challenges

- Confirm every field unit before comparison:
  `V/cm -> V/m`, `cm^-1 -> m^-1`, `A/cm2`, and mesh `m/um`.
- Confirm electron/hole conventional-current orientation and distinguish
  native Sentaurus node/element observations from projected box/local-edge
  values.
- Confirm the 160 zero-volume local edges are geometric zeros and never
  included as finite dex samples.
- Confirm fixed-imported-state evidence is not used as proof of a
  self-consistent branch.
- Confirm the source-unit result is independent of the mobility experiment.
- Reject any attempt to fit alpha, saturation velocity, field scale, or
  geometry to compensate for current-support error.

## Evidence roots

- `pn2d-minimal6-sentaurus-highfield-drive-20260724-a/b`
- `pn2d-minimal6-highfield-box-current-20260724-a/b`
- `pn2d-minimal6-current-factorial-20260724-a/b`
- `pn2d-minimal6-phase-e-highfield-followup-20260724-a/b`
- `pn2d-minimal6-electric-field-phase-f-triangle-source-20260724-a/b`
- `pn2d-minimal6-impact-factorization-final-20260724-a/b`
- `pn2d-source-unit-isolated-audit-20260724-a/b`

The review should stop with `insufficient_data` rather than infer a proprietary
Sentaurus directed-edge operator from projected node or element fields.
