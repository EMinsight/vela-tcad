# PN2D Minimal6 state-unit and transport re-audit

Date: 2026-07-23

Status: authoritative quantitative correction for the 2026-07-23 Minimal6
QFP, mobility, and directed-edge reports.

## Corrections

The prior inverse inputs treated Vela restart columns named `electrons_m3`
and `holes_m3` as SI even though unit-scaling mode wrote internal cm^-3. The
sealer then divided by 1e6, creating an artificial six-decade density gap.
The re-audit fixes three contracts:

- restart CSV I/O converts between internal concentration and physical m^-3;
- the `old_slotboom` reference concentration is converted through the active
  unit system; and
- Minimal6 uses the Sentaurus-2018 model-derived silicon value
  `ni = 1.4638914958767616e10 cm^-3`.

The fixed-state audit now loads the deck `materials_file`. Offline QFP-SG
replay recovers frozen node effective-ni from `(psi, QFP, n/p)`, including
BGN. Production SG and mobility formulas were not changed.

## Regeneration and sealing

The corrected Vela sweep accepted 40/40 states with zero failed transitions.
At -1 V, the restart CSV stores contact majority density near `1e23 m^-3` and
minority density `2.7411172687e9 m^-3`. Vela, Sentaurus, and supplemental
sealed roots each contain 40 exact states and valid seals. Sentaurus inputs
were independently reimported with sealed importer SHA-256
`8238ab6db4bae2167a8bd7c385edf44a5016568238afe16326c6076fb8752f86`.

## Corrected node comparison

Both topologies and reverse biases -1 through -20 V are pooled. Potential
entries are absolute differences; densities are absolute log10 ratios.

| Quantity/support | n | Median | p95 | Maximum |
|---|---:|---:|---:|---:|
| Electrostatic potential, all nodes (V) | 240 | 1.22326e-7 | 1.22326e-7 | 1.22326e-7 |
| Electrostatic potential, nodes 1/5 (V) | 80 | 3.75211e-12 | 1.46025e-11 | 2.24802e-11 |
| Electron QFP, nodes 1/5 (V) | 80 | 0.331106 | 0.347977 | 0.367601 |
| Hole QFP, nodes 1/5 (V) | 80 | 0.329839 | 0.347838 | 0.367407 |
| Electron density, nodes 1/5 (dex) | 80 | 5.56234 | 5.84576 | 6.17544 |
| Hole density, nodes 1/5 (dex) | 80 | 5.54106 | 5.84343 | 6.17217 |
| Electron majority contact (dex) | 80 | 4.34e-16 | 4.34e-16 | 4.34e-16 |
| Hole majority contact (dex) | 80 | 5.71e-15 | 2.70e-14 | 3.21e-14 |
| Electron minority contact (dex) | 80 | 6.48117e-6 | 6.48117e-6 | 6.48117e-6 |
| Hole minority contact (dex) | 80 | 6.48117e-6 | 6.48117e-6 | 6.48117e-6 |

The global six-decade failure is removed: contact states agree. The remaining
density gap is localized to internal nodes 1 and 5 and accompanies an
approximately 0.33 V internal QFP gap. Electron- and hole-derived Vela
effective-ni agree within `1.93e-16 dex`; the BGN range is `1.65563e16` to
`1.96229e16 m^-3`.

## QFP-SG replacement

Offline SG closes against C++ on all 720 baseline carrier-edge samples, with
maximum relative error `3.35492e-16` against a `5e-11` gate.

| Carrier | Branch | Median/p95 error (dex) | Sign agreement | Paired improvement |
|---|---|---:|---:|---:|
| Electron | Vela baseline | 3.66595 / 6.47450 | 86.43% | 0 |
| Electron | Sentaurus eQFP at nodes 1/5 | 1.93702 / 2.36927 | 100% | 1.75184 dex |
| Hole | Vela baseline | 3.68958 / 6.48495 | 86.43% | 0 |
| Hole | Sentaurus hQFP at nodes 1/5 | 1.90704 / 2.42772 | 100% | 1.80582 dex |

QFP replacement fixes sign and removes about 1.75-1.81 dex, but does not
close magnitude alone.

## Conditional directed-edge inversion

Sentaurus exposes six node current vectors, not nine native directed-edge
fluxes. The following uses deterministic node-to-edge reconstructions.

| Support | Carrier/form | Vela median | Full replacement median |
|---|---|---:|---:|
| Endpoint mean | Electron QFP SG | 3.66595 | 0.127685 |
| Endpoint mean | Hole QFP SG | 3.68958 | 0.127541 |
| Endpoint mean | Electron density SG | 3.66623 | 0.127680 |
| Endpoint mean | Hole density SG | 3.68984 | 0.127537 |
| Adjacent-cell mean | Electron QFP SG | 3.64491 | 0.179749 |
| Adjacent-cell mean | Hole QFP SG | 3.66941 | 0.179558 |

Endpoint QFP and mobility paired reductions are 1.75-1.81 dex and 1.75-1.80
dex, respectively. Replacing electrostatic potential changes the median by
only about `2.2e-8 dex`. Required mobility is within `0.1275-0.1277 dex` of
the exported Sentaurus mobility under the endpoint reconstruction.

## Native element result and decision

On identical four-element Sentaurus support, current and negative QFP
gradient remain collinear: electron angle median/p95 is
`0.008824/0.351724 deg`; hole is `0.011670/0.350602 deg`.

No native element carrier density is exported. Density inferred from
`J = q mu n GradQF` differs from arithmetic node-to-cell controls by medians
of 6.340 dex (electron) and 5.923 dex (hole); this is a support diagnostic,
not an observed element-density mismatch.

The evidence establishes that potential already agrees, the remaining state
gap is internal, QFP corrects sign and much of magnitude, and Sentaurus
mobility plus state replacement reduces the conditional reconstruction to
0.13-0.18 dex. A production formula change is still unauthorized because no
native Sentaurus directed-edge flux or native element density is observable.

## Evidence

- `build-release/pn2d-minimal6-task8-vela-unitfix-20260723-b`
- `build-release/pn2d-minimal6-inverse-inputs-unitfix-20260723-b`
- `build-release/pn2d-minimal6-physics-inverse-audit-unitfix-20260723-b`
- `build-release/pn2d-minimal6-mobility-diagnosis-unitfix-20260723-b`
- `build-release/pn2d-minimal6-qfp-sg-replacement-unitfix-20260723-b`
- `build-release/pn2d-minimal6-edge-flux-inversion-unitfix-20260723-b`
- `build-release/pn2d-minimal6-transport-element-closure-unitfix-20260723-b`

Validation passed: 68 focused Python tests, all independent artifact
verifiers, full Release build, and all 465 CTest tests.
