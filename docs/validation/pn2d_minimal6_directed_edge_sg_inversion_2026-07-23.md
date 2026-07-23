# PN2D Minimal6 directed-edge SG inversion audit

Date: 2026-07-23

## Question

Determine whether the remaining current mismatch is primarily caused by:

1. mapping Sentaurus node current vectors to Vela directed edges;
2. the mobility definition;
3. electron/hole QFP sign conventions;
4. endpoint carrier densities; or
5. a different discrete current operator.

No production formula in `include/` or `src/` was modified.

## Input and support contract

- Exact state lattice: 40 states, comprising `mirror` and `sketch` at reverse
  biases -1 V through -20 V.
- Canonical mesh support: 6 nodes, 4 triangles, and 9 sorted directed edges.
- Sentaurus TDR audit: all 40 states contain one two-component
  `eCurrentDensity` and one two-component `hCurrentDensity` dataset with six
  values and `global_vertex_order` mapping.
- Native Sentaurus directed-edge current is therefore unavailable. The
  observed current is a node vector, not a nine-edge flux.
- The P1 line mean along an edge is exactly the endpoint mean. It is not an
  independent current-support reconstruction.
- Two non-fitted reconstructions were compared:
  `endpoint_mean_tangent` and `adjacent_cell_mean_tangent`.
- `qfp_aligned_endpoint_magnitude_control` removes the observed vector
  direction from the comparison. It is a localization control and is
  ineligible for formula acceptance.

The 40-state audit produced 720 support rows, 34,560 staged SG samples, and
4,320 effective-mobility inversion samples.

## Staged replacement result

The table below uses the endpoint-mean tangent reconstruction and all
nonzero-current edges. Error is
`abs(log10(abs(candidate/reference)))`.

| Carrier | SG form | Branch | Median error (dex) | p95 error (dex) | Sign agreement |
|---|---|---|---:|---:|---:|
| electron | QFP | Vela all | 3.500744 | 6.297851 | 0.864286 |
| electron | QFP | Sentaurus QFP only | 2.101840 | 2.751742 | 1.000000 |
| electron | QFP | Sentaurus QFP + mobility | 0.190124 | 0.946733 | 1.000000 |
| electron | QFP | Sentaurus all | 0.191079 | 0.732083 | 1.000000 |
| hole | QFP | Vela all | 3.524351 | 6.308279 | 0.864286 |
| hole | QFP | Sentaurus QFP only | 2.071857 | 2.810195 | 1.000000 |
| hole | QFP | Sentaurus QFP + mobility | 0.185204 | 1.100887 | 1.000000 |
| hole | QFP | Sentaurus all | 0.186200 | 0.886237 | 1.000000 |
| electron | density | Vela all | 2.497817 | 8.750375 | 0.864286 |
| electron | density | Sentaurus density only | 1.809082 | 2.315188 | 1.000000 |
| electron | density | Sentaurus density + mobility | 0.127914 | 0.510179 | 1.000000 |
| electron | density | Sentaurus all | 0.127680 | 0.513124 | 1.000000 |
| hole | density | Vela all | 2.474210 | 8.808828 | 0.864286 |
| hole | density | Sentaurus density only | 1.779098 | 2.373641 | 1.000000 |
| hole | density | Sentaurus density + mobility | 0.127747 | 0.664333 | 1.000000 |
| hole | density | Sentaurus all | 0.127537 | 0.667278 | 1.000000 |

The sequential paired median reductions are:

| Carrier | SG form | State-variable step (dex) | Mobility step (dex) | Potential step (dex) |
|---|---|---:|---:|---:|
| electron | QFP | 1.422250 | 1.837320 | 7.63e-11 |
| hole | QFP | 1.476220 | 1.793500 | -7.63e-11 |
| electron | density | 0.857421 | 1.676290 | 6.59e-11 |
| hole | density | 0.924478 | 1.636320 | -6.68e-11 |

The medians are paired per sample and are not algebraically additive after
aggregation. The important ordering is stable: replacing the carrier state
fixes the sign, and replacing mobility removes about another 1.6-1.8 dex.
Replacing electrostatic potential after those steps has no material effect.

## Support-mapping sensitivity

For the full Sentaurus replacement:

| Carrier | SG form | Endpoint mean median (dex) | Adjacent-cell median (dex) | Difference (dex) |
|---|---|---:|---:|---:|
| electron | QFP | 0.191079 | 0.227826 | 0.036747 |
| hole | QFP | 0.186200 | 0.224873 | 0.038673 |
| electron | density | 0.127680 | 0.179744 | 0.052064 |
| hole | density | 0.127537 | 0.179553 | 0.052016 |

The choice between the two deterministic node-to-edge mappings changes the
pooled median by only about 0.04-0.05 dex. It cannot explain the previous
approximately 2 dex residual.

Localization by edge class shows that the remaining discrepancy is larger on
boundary-to-interior edges:

| Carrier | SG form | Internal edge median (dex) | Boundary-to-interior median (dex) |
|---|---|---:|---:|
| electron | QFP | 0.165078 | 0.398243 |
| hole | QFP | 0.165222 | 0.435103 |
| electron | density | 0.127680 | 0.254503 |
| hole | density | 0.127537 | 0.295181 |

The adjacent-cell QFP reconstruction gives an even smaller internal-edge
median, about 0.113 dex for both carriers, but a larger pooled median. This
pattern is consistent with boundary/current-support semantics contributing
to the remaining sub-dex mismatch.

## Effective mobility inversion

For each edge, the full Sentaurus state was evaluated at unit mobility. The
positive mobility required to match each reconstructed current was then
computed exactly from SG linearity in mobility. Zero operators and
sign-incompatible samples remain typed and are never converted to negative
physical mobility.

On endpoint-mean tangent support, 280 of 360 samples per carrier and SG form
have a nonzero operator. None is sign-incompatible; the other 80 are zero
operators on zero-driving-force edges.

| Carrier | SG form | Median required/Sentaurus mobility (dex) | Median required/Vela production mobility (dex) |
|---|---|---:|---:|
| electron | QFP | +0.191079 | +2.102530 |
| hole | QFP | +0.186200 | +2.072547 |
| electron | density | -0.101679 | +1.809772 |
| hole | density | -0.106559 | +1.779788 |

Thus the mobility required by the reconstructed Sentaurus current is within
about 0.1-0.2 dex of the exported Sentaurus mobility, but is about 1.8-2.1 dex
above the effective Vela production edge mobility on the pooled nonzero
support. This quantitatively identifies mobility/operator-state interaction
as the dominant source of the previous approximately 2 dex magnitude
residual.

QFP-SG and density-SG do not yield identical required mobility because the
exported Sentaurus QFP, density, potential, and fixed `ni=1e16 m^-3` are not
an exactly self-consistent endpoint state under Vela's equilibrium
transformation. Their 0.06-dex difference is evidence for a remaining state
definition/operator convention, not evidence for a sign error.

## Scientific conclusion

1. The earlier near-90-degree pooled direction error was not a mobility
   problem; with carrier conventions separated, both carriers have the same
   accepted negative-QFP-gradient sign.
2. QFP or density replacement alone corrects current direction but does not
   close magnitude.
3. Mobility replacement after the carrier-state replacement is the dominant
   closure step and reduces the median mismatch to 0.13-0.19 dex.
4. Electrostatic-potential replacement is immaterial on this state lattice.
5. Endpoint versus adjacent-cell support mapping contributes only about
   0.04-0.05 dex to the pooled median, although boundary-to-interior edges
   retain a larger sub-dex mismatch.
6. The available evidence supports investigating why Vela's production edge
   mobility is 1.8-2.1 dex below the mobility required by the Sentaurus
   current/state pair. It does not support changing the QFP sign or SG
   formula.
7. A production formula change remains unauthorized because no native
   Sentaurus directed-edge flux is observable.

## Remote probe and remaining gate

A read-only retry to the authorized Sentaurus host on 2026-07-23 timed out
while connecting to port 22. No remote deck or result was changed.

When the host is reachable, the next gate is a single-state O-2018.06-SP2
probe for an element/edge current-output qualifier. A 40-state regeneration
is justified only if that probe produces a dataset with native cell or
nine-edge support. If it again produces six global-node values, the native
edge-flux question remains unobservable and no additional 40-state export is
scientifically useful.

## Evidence

- Generated root:
  `build-release/pn2d-minimal6-edge-flux-inversion-20260723-a`
- Main manifest: `manifest.json`
- Raw replacement lattice: `sg_replacement_samples.csv`
- Replacement summary: `sg_replacement_summary.csv`
- Sequential paired contributions: `replacement_contributions.csv`
- Required mobility samples: `mobility_inversion_samples.csv`
- Required mobility summary: `mobility_inversion_summary.csv`
- Independent verification: `independent_verification.json`

The standalone verifier passed with zero failures and independently checked
all output hashes, raw residuals, grouped medians, paired contributions,
mobility classifications, and required-mobility algebra.
