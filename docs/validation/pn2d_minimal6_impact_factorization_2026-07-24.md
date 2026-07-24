# PN2D Minimal6 impact-source factorization

Date: 2026-07-24

Status: complete with typed outcome `current_support_dominant`.

## Scope

The audit uses the final electric-field-mobility self-consistent states but
does not change the impact configuration. For every exact
`mirror/sketch x -1..-20 V` state it compares on one explicit triangle
local-edge support:

- Sentaurus native nodal `ImpactIonization`;
- Sentaurus nodal `alpha * |J| / q`;
- the same Sentaurus fields projected to Vela triangle local edges;
- the production Vela triangle-GSS source;
- a Sentaurus-alpha-only hybrid with Vela current proxy; and
- a Sentaurus-current-only hybrid with Vela alpha.

Sentaurus local-edge current remains a declared endpoint-to-triangle
projection, not a native directed-edge observation.

## Deterministic evidence

Roots:

- `build-release/pn2d-minimal6-impact-factorization-final-20260724-a`
- `build-release/pn2d-minimal6-impact-factorization-final-20260724-b`

The independent verifier checks every local formula, state aggregation, hash,
geometric zero, source-to-node mapping, and A/B output equality.

| Metric | Count | Median dex | P95 dex | Maximum dex |
|---|---:|---:|---:|---:|
| Sentaurus nodal alpha-current vs native source | 40 | `5.25100e-5` | `3.22301e-4` | `4.23421e-4` |
| node-to-triangle support projection | 40 | `2.58704e-10` | `2.15771e-8` | `1.46178e-7` |
| fully projected Sentaurus source | 40 | `5.25103e-5` | `3.22280e-4` | `4.23275e-4` |
| Vela production triangle source | 40 | `10.764920` | `11.904647` | `11.950516` |
| Sentaurus alpha + Vela current proxy | 40 | `12.360983` | `12.384311` | `12.385500` |
| Vela alpha + Sentaurus current projection | 40 | `0.435579` | `0.498671` | `0.498926` |
| local alpha | 640 | `0.404687` | `35.400546` | `244.844303` |
| local current proxy | 480 | `5.333305` | `12.838626` | `12.850894` |

The very large alpha tail occurs where both coefficients are exponentially
small. It does not overturn the paired source experiment: replacing only the
current proxy removes more than ten decades of integrated-source error,
whereas replacing only alpha does not close the source.

## Geometry and mapping

The 40-state lattice contains 160 zero-volume triangle local edges, or 320
carrier-local rows. Every zero-volume source remains exactly zero.

| Closure | Maximum error |
|---|---:|
| sum of three local truncated volumes vs triangle area | `0.0` relative |
| local source split to two endpoint nodes vs local total | `1.50363e-16` relative |
| Sentaurus nodal alpha-current reconstruction | `4.23421e-4 dex` |
| Sentaurus node-to-triangle projection | `1.46178e-7 dex` |

This confirms the user's right-triangle interpretation for the present mesh:
the circumcenter lies on the hypotenuse, so that diagonal has zero truncated
partial volume; the other two local edges partition the triangle area.

## Production decision

The current production source mapping is explicitly named and validated as
`triangle_gss_gradqf_truncated`. It is a QFP-gradient-only design; it does not
promise an `electric_field` alpha branch. Adding electric-field alpha on this
support would require a separately named mapping mode and tests.

That design must not be implemented yet. Task 9 requires current support and
sign to close before alpha parity is claimed. The production triangle current
proxy differs from the Sentaurus projected current by a median `5.33 dex` and
dominates the source difference. No impact coefficient, field scale, geometry
factor, or alpha law is fitted or changed.

The next production experiment should replace the triangle
`mobility * midpoint_density * edge_QFP_difference` proxy with a diagnostic
box-current support that is independently closed to terminal current, while
holding alpha and geometry fixed. Only after that gate passes should a
separate electric-field triangle-alpha mode be considered.
