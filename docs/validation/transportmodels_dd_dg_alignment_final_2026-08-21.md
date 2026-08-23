# TransportModels DD/DG alignment final report

Status: **completed with production regression passing**.

## Final production result

The qualified production configuration is:

- OldSlotboom BGN with `fermi_statistics_correction=true`;
- Fermi-Dirac carrier statistics;
- density-gradient `direct_band_edge` compatibility coupling;
- all-material `sentaurus_box` DG operator with neutral continuous interface;
- Masetti + Lombardi surface mobility and `transport_cell_vector`
  quasi-Fermi-gradient high-field drive.

The post-P2 regression completed all 42 points.

| Curve/region | Metric | Result | Limit | Status |
|---|---|---:|---:|---|
| Id-Vg transition | maximum absolute log error | 0.04414 dex | 0.15 dex | pass |
| Id-Vg on | maximum relative error | 2.9405% | 10% | pass |
| Id-Vd full curve | maximum relative error | 2.8814% | 5% | pass |
| Id-Vd Vd=2 V | endpoint relative error | 1.3911% | 3% | pass |

The three deepest-off Id-Vg points remain dominated by the different numerical
current floors (Sentaurus approximately 1e-15 A/um and Vela approximately
1e-21 A/um); this produces a large log error but is outside the transport
acceptance regions.

## Ordered task outcomes

1. Same-settings DD rerun: 21/21 points converged. The old DD baseline setup
   was not the main cause of the DG discrepancy.
2. Sentaurus default/DirectQC by default/NoFermi 2x2: DirectQC changes Id by
   only 0.001-0.005 dex; Fermi-BGN changes threshold-region Id by up to
   approximately 0.205 dex.
3. Vela Fermi-BGN A/B: the Vela correction reproduces the Sentaurus direction
   and magnitude and is the change responsible for the final curve pass.
4. Five-bias spatial audit: source-end and channel-mid Qn/n profiles pass; the
   remaining local hotspot is the drain-end Si/SiO2 corner.
5. Sentaurus exponential DG transport coupling: parser, density relation,
   generalized-Einstein separation, contact current, and tests were implemented
   as an opt-in experimental mode. It is not production-qualified: the
   self-consistent outer loop diverges and Frozen-Q gives only 4/5 converged
   points with excessive low-Vg current changes.
6. Enormal/high-field audit: VTK now exports the exact per-cell surface-normal
   field and an area-weighted nodal reconstruction. Existing
   `transport_cell_vector` high-field drive is retained because the evidence
   does not justify a contact-specific rule change.
7. P3: all 23 DG manufactured/oracle tests pass (97 assertions). WKB remains a
   separate semiconductor-only branch because it is mutually exclusive with
   the qualified all-material domain. `sentaurus_box + neutral_continuous`
   remains the best self-consistent operator/interface contract. Global mesh
   refinement is deferred; a later study should remesh only the drain-end
   corner in both Sentaurus and Vela.

## Verification

- Density-gradient unit tests: 23 cases / 97 assertions passed.
- Carrier-statistics tests: 6 cases / 188 assertions passed.
- DG configuration/density/Jacobian tests: 5 selected cases / 68 assertions passed.
- VTK mobility/normal-field integration test: 1 case / 11 assertions passed.
- Sentaurus spatial oracle, Vela spatial oracle, spatial comparison, and final
  regression hash checks all passed.

## Primary artifacts

- `transportmodels_dg_post_p2_regression_v4_2026-08-21.{json,md}`
- `transportmodels_dg_transport_coupling_ab_2026-08-21.{json,md}`
- `transportmodels_idvg_spatial_oracle_2026-08-21.{json,md}`
- `transportmodels_dg_p3_assessment_2026-08-21.md`
- `dg_post_p2_regression_v4_2026-08-21/dg_idvg_idvd_comparison.{png,svg}`
