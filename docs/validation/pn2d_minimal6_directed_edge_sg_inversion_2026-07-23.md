# PN2D Minimal6 directed-edge SG inversion audit

> Quantitative correction: the node-to-edge tables below used the pre-fix
> Vela restart-unit and fixed-`ni` inputs. They are superseded by
> `pn2d_minimal6_state_unit_and_transport_reaudit_2026-07-23.md`. The native
> Sentaurus element observations remain valid.

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

The native-element regeneration below supersedes the node-vector
reconstruction for conclusions about current direction. The staged
node-to-edge replacement remains useful as a conditional sensitivity audit.

1. The earlier near-90-degree pooled direction error was not a mobility or
   QFP-sign problem. On identical native element support, current and
   negative QFP gradient are aligned to much less than one degree.
2. QFP or density replacement alone corrects current direction but does not
   close magnitude.
3. Mobility replacement after the carrier-state replacement is the dominant
   closure step inside the conditional node-current reconstruction.
4. Electrostatic-potential replacement is immaterial on this state lattice.
5. Endpoint versus adjacent-cell support mapping contributes only about
   0.04-0.05 dex inside that reconstruction, although boundary-to-interior
   edges retain a larger sub-dex mismatch.
6. The apparent 1.8-2.1 dex Vela mobility gap is conditional on a
   reconstructed node-current reference and must not be interpreted as an
   independently established mobility defect.
7. A production formula change remains unauthorized because no native
   Sentaurus directed-edge flux is observable.

## Native element regeneration

After the authorized Sentaurus host became reachable, a single-state
O-2018.06-SP2 probe established the following native element fields:

- `ElectricField/Element/Vector`;
- `eGradQuasiFermi/Element/Vector` and
  `hGradQuasiFermi/Element/Vector`;
- `eMobility/Element` and `hMobility/Element`; and
- `eCurrentDensity/Element/Vector` and
  `hCurrentDensity/Element/Vector`.

Each field has exactly four values in `region_cell_order`, one per triangle.
The requested element potential, carrier densities, and QFP scalars are not
available. Native directed-edge flux is also still unavailable.

All 40 states were regenerated: 20/20 mirror and 20/20 sketch passed. The
downloaded archive SHA-256 is
`34eb46ae228a0477482209ef540356539e80eff9f862e6f20ec1324658a75583`.
All 40 TDR files independently matched the remote SHA-256 ledger.

## Same-element transport result

The Sentaurus manual defines `eGradQuasiFermi` and `hGradQuasiFermi` as
negative QFP gradients. Under the active isothermal drift-diffusion model and
Einstein relation, both carriers therefore obey the observable direction
contract

`J = q * mobility * density * GradQuasiFermi`.

The current, mobility, and QFP-gradient fields were compared on identical
native element support. Carrier density cannot be checked on that same
support, so an effective density was inferred by least-squares projection and
compared only to diagnostic node-to-cell density controls.

| Carrier | Valid | Current/gradient angle median / p95 (deg) | Orthogonal residual median / p95 | Effective vs arithmetic node-density gap median / p95 (dex) | Effective vs geometric node-density gap median / p95 (dex) |
|---|---:|---:|---:|---:|---:|
| electron | 160/160 | 0.008824 / 0.351724 | 0.000154 / 0.006139 | 6.340327 / 7.253257 | 4.968272 / 7.397695 |
| hole | 160/160 | 0.011670 / 0.350602 | 0.000204 / 0.006119 | 5.923051 / 6.927266 | 4.029447 / 6.903342 |

There are no sign-incompatible or degenerate samples. This directly rejects
mobility and carrier-sign convention as explanations for the earlier
near-90-degree direction mismatch. The magnitude does not close against any
node-density interpolation: the element current's magnitude support cannot be
combined with node density as if they were the same discrete operator.

Two output roots were generated independently and are byte-identical. A
standalone verifier recomputed all 320 carrier samples without importing the
main analysis implementation. It found zero density or orthogonal-residual
reconstruction error and a maximum angle reconstruction difference of
`1.454e-12` degrees.

## Remaining gate

The next scientifically discriminating input is a native directed-edge
carrier flux or an independently documented conversion from Sentaurus's
element current representation to the SG edge operator. Without it, neither
the remaining magnitude gap nor a production formula replacement is
identifiable.

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
- Native transport source and export:
  `build-release/pn2d-minimal6-transport-elements-20260723-b`
- Same-element deterministic audits:
  `build-release/pn2d-minimal6-transport-element-closure-20260723-b/audit-a`
  and `audit-b`
- Same-element independent verification:
  `build-release/pn2d-minimal6-transport-element-closure-20260723-b/independent-verification-a.json`
  and `independent-verification-b.json`

The original standalone verifier passed with zero failures and independently
checked all output hashes, raw residuals, grouped medians, paired
contributions, mobility classifications, and required-mobility algebra. The
new same-element verifier independently passed all 320 carrier samples.
