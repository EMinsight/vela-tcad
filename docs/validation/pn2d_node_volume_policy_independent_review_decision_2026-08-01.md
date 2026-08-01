# PN2D node-volume policy independent review decision

Date: 2026-08-01

## Outcome

Both required independent reviews are complete:

| Review | Verdict | Meaning |
|---|---|---|
| scientific | `APPROVE_WITH_CONDITIONS` | supports a qualified non-obtuse PN2D M0/M2 template proposal |
| code | `APPROVE_WITH_CONDITIONS` | supports implementing an atomic patch, but does not approve the current worktree |

Combined typed outcome:

```text
independent_reviews_complete_atomic_patch_and_rereview_required
```

Production default change authorized now: **no**.

## Why the reviews do not yet authorize the switch

The scientific evidence is strong for the actual M0/M2 grid family: direct
Sentaurus box measures match Vela mixed-Voronoi to floating-point roundoff,
off/IIC/on and forward controls pass, and both meshes are non-obtuse. However,
the result is not a general Sentaurus MixAverage equivalence claim for arbitrary
meshes.

The current repository contains no actual default-value patch. The BV template
still defaults to the legacy avalanche profile, does not render a node-volume
policy, and the acceptance aggregator is not bound to the default render. The
code reviewer therefore had no final patch to approve.

## Authorized next step

The next authorized task is narrowly limited to:

1. create an isolated atomic PN2D BV template patch;
2. bind SG/Laux and mixed-Voronoi through one profile selector;
3. retain an explicit legacy+barycentric rollback;
4. add fail-closed half-migration, render-binding, mutation, IV-isolation, and
   compatibility tests;
5. limit or qualify the template for non-obtuse PN2D meshes;
6. rerun the frozen acceptance against the actual default render;
7. conduct a fresh code review of the patch and a scientific scope check of
   the new bound evidence.

No global C++ parser default, PN2D IV template, solver parameter, physical
model coefficient, or acceptance threshold is authorized to change.

## Review records

- `docs/validation/pn2d_node_volume_policy_independent_scientific_review_2026-08-01.md`
- `docs/validation/pn2d_node_volume_policy_independent_code_review_2026-08-01.md`
