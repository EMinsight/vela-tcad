# Slot-LDMOS Vela IALMob ablation status (2026-08-20)

## Controlled comparison

The Vela pair shares the sealed Slot-LDMOS mesh, Stage 05 converged state,
Van Overstraeten self-consistent avalanche, frozen avalanche-source Jacobian,
Masetti doping dependence, quasi-Fermi-gradient high-field saturation,
external resistance, contacts, tolerances, and direct-Newton handoff.  The only
physics delta is:

- off: `mobility.model = masetti_field`
- on: `mobility.model = masetti_field_lombardi` with
  `surface.surface_interface = [Silicon_1, Oxide_1]`

The selected regions share 255 mesh edges.  Deck generation rejects a missing
or empty interface and audits the normalized pair for changes outside mobility
and isolated output paths.

## 60 V outer-load-line result

| Metric | IALMob off | IALMob on | On - off |
|---|---:|---:|---:|
| Inner voltage (V) | 0.733397745 | 0.806855273 | +0.073457527 |
| Drain current (A/um) | 5.926660216e-11 | 5.919314450e-11 | -0.123944% |
| Integrated avalanche qG (A/um) | 6.201005213e-11 | 6.208662222e-11 | +0.123480% |
| Maximum electric field (V/cm) | 2.829837753e4 | 3.113324859e4 | +10.017787% |
| Maximum avalanche generation (cm^-3 s^-1) | 5.415392742e23 | 5.833983893e23 | +7.729655% |
| Nodal avalanche peak (x, y), um | (0.00466258, 3.49729) | (0.00466258, 3.49729) | no shift |

Both terminal solutions satisfy the 60 V, 1e12 ohm*um load line and global
continuity closure.  The on case required six scalar boundary evaluations; all
six device solves converged.

## Implementation findings

The first full-mesh Lombardi run exposed two post-processing defects:

1. Interface-distance preprocessing admitted every two-cell edge before the
   per-cell selector, causing both excessive cost and internal same-region
   edges to enter the nearest-interface search.  Geometry candidates are now
   filtered to the configured physical region pair.
2. VTK/release diagnostics used the unresolved input mobility config rather
   than a config carrying the surface geometry cache.  The output path now
   resolves that cache once and uses it for current-density avalanche and
   mobility fields.  Nodal mobility evaluation also consumes the cached normal
   field and distance.

After these changes, checkpoint-resumed IALMob post-processing completed in
about six seconds instead of continuing for more than ten minutes.

## BVDS status

The 60 V probe reaches only about 5.9e-11 A/um, so it cannot define the
1e-7 A/um BVDS criterion.  A paired fixed-inner-voltage continuation was
prepared and attempted.  The off case converged at 0.85 V, but the transition
to the next target invoked repeated quasi-Fermi branch recovery after changing
from the external-resistor state to a fixed-voltage boundary.  This pathway
was stopped because it measures boundary-reference recovery cost rather than
IALMob physics.

Consequently the IALMob BVDS shift remains unmeasured.  The next valid route is
to continue both cases under the external-resistor boundary, with an actual
inner-voltage step cap of at most 1 V near the previously observed 15 V solver
limit, then interpolate the 1e-7 A/um crossing.  No BVDS number should be
reported until both curves cross that criterion.

## External-resistor continuation execution

The strict external-resistor pair was prepared with the shared
`1e12 ohm*um` resistor, an inner-voltage request cap of `1 V`, and an opt-in
adaptive device continuation that halves a failed device step.  The off case
resumed the preserved boundary checkpoints and added one converged evaluation:

| Case | Inner voltage (V) | Drain current (A/um) | Target outer voltage (V) | Load-line residual (V) | Newton iterations |
|---|---:|---:|---:|---:|---:|
| IALMob off | 13.233397745426 | 5.357545312919e-10 | 1000 | -451.012070963 | 65 |

The following requested steps all failed with `line_search_non_decrease` at
Newton iteration 2, while retaining the 13.233397745426 V checkpoint:

| Requested step (V) | Failed inner voltage (V) |
|---:|---:|
| 0.5 | 13.733397745426 |
| 0.25 | 13.483397745426 |
| 0.125 | 13.358397745426 |

This is a voltage-controlled branch limit, not evidence that the
`1e-7 A/um` current threshold has been reached.  The off current remains about
187 times below that criterion, so the on curve was not launched past 60 V and
no IALMob BVDS interpolation was performed.

## Fold-crossing ablations

Pseudo-arclength was seeded by the converged off states at 12.983397745426 V
and 13.233397745426 V.  The following corrector ablations were executed:

| Corrector | Result | Evidence |
|---|---|---|
| Frozen avalanche-source Jacobian, strict line search | Fail | At `Delta s=1e-4`, iteration 3 rejected all updates; residual 1.374941573e-4 and best trial 1.374959138e-4. |
| Finite-difference avalanche-source Jacobian only in arclength | Fail | The augmented Jacobian solve failed before the first corrector update. |
| Frozen Jacobian, relative non-decrease allowance 1e-4 | Fail | After 80 iterations, residual drifted to 1.383309348e-4 instead of approaching the unchanged 1e-8 tolerance. |

The relaxed line-search candidate is therefore rejected.  The production
physics remains self-consistent avalanche with the frozen source Jacobian.
The next solver gate is either a genuinely coupled series-resistor formulation
that adds the inner drain voltage and load-line equation to the Newton unknowns,
or a nonsingular, descent-producing avalanche-source Jacobian suitable for the
arclength bordered system.  Until one of those gates closes, both external-load
curves cannot be advanced to the `1e-7 A/um` crossing.
