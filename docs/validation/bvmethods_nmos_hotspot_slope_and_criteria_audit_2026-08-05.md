# BVmethods NMOS hotspot slope and criterion audit (2026-08-05)

## Corrected headline

The previously quoted `6.3774942778 V` value is not the
`BreakAtIonIntegral(3 1.)` path threshold. It is the linear zero crossing of
`Iava-Id` extracted from the original sparse ABA-coupled plot between
5.939777636 V and 7.107808235 V.

The O-2018.06-SP2 manual states that
`BreakAtIonIntegral(number value)` stops a quasistationary sweep when the
specified number of largest, decreasingly ordered path ionization integrals
exceed the value. In the BVmethods deck, `BreakAtIonIntegral(3 1.)` terminates
the original path calculation at 10.448266731 V. It does not define the
6.377494278 V current-source interpolation.

The dense fixed-bias Sentaurus replay resolves the curvature hidden by the
sparse adaptive output and gives an `Iava-Id=0` crossing of 6.734425890 V.
At 6.4 V the sparse straight line predicts `Iava-Id=+4.0965e-8 A/um`, whereas
the direct fixed-bias calculation gives `-5.8913e-7 A/um` and
`Iava/Id=0.937893`. The 0.357 V change from 6.377494 V to 6.734426 V is
therefore a sparse interpolation error, not a difference between path and
volume-integral breakdown physics.

## Electron hotspot amplitude and slope

The comparison uses the direct peak of `eImpactIonization` in the Sentaurus
silicon node field and the direct peak of
`electron_alpha_m_inv * electron_flux_proxy` in the Vela SG edge ledger. The
two peak locations agree within 0.000147 um at every common state from 6.4 to
7.0 V.

| bias (V) | peak G Vela/Sentaurus | peak alpha Vela/Sentaurus | peak Jn Vela/Sentaurus | integrated electron source Vela/Sentaurus |
|---:|---:|---:|---:|---:|
| 6.4 | 0.917741 | 1.009512 | 0.935737 | 0.981941 |
| 6.5 | 0.917185 | 1.009422 | 0.935160 | 0.979853 |
| 6.6 | 0.916696 | 1.009334 | 0.934649 | 0.977583 |
| 6.7 | 0.916271 | 1.009249 | 0.934200 | 0.975210 |
| 6.8 | 0.915908 | 1.009166 | 0.933813 | 0.972807 |
| 6.9 | 0.915605 | 1.009086 | 0.933484 | 0.970435 |
| 7.0 | 0.915359 | 1.009008 | 0.933212 | 0.968144 |

The 6.4-to-7.0 V log slopes are:

| quantity | Sentaurus (dex/V) | Vela (dex/V) | Vela - Sentaurus (dex/V) |
|---|---:|---:|---:|
| peak electron generation | 0.386632 | 0.384751 | -0.001881 |
| peak electron alpha | 0.017790 | 0.017429 | -0.000362 |
| peak electron current density | 0.369279 | 0.367323 | -0.001957 |
| integrated electron source | 0.490555 | 0.480313 | -0.010242 |

Therefore, the hotspot amplitude grows at essentially the correct rate. Its
approximately 8.3% fixed amplitude deficit is almost entirely the electron
current-density deficit; the local electron ionization coefficient is already
about 0.9% high and must not be refitted.

The larger integrated-source slope difference is spatial. Defining an
effective generation support as `integrated electron source / peak electron
generation`, the Vela/Sentaurus support ratio decreases from 1.069954 at
6.4 V to 1.057666 at 7.0 V. The remaining bias-dependent source discrepancy is
therefore in the widening/redistribution of the high-generation region, not in
the peak coefficient law.

## Path criterion versus current-source criterion

The exact Sentaurus `MeanIonIntegral` field contains eight distinct positive
plateaus in the 6.32-to-7.0 V states. The top values are:

| bias (V) | rank 1 | rank 2 | rank 3 | Iava/Id |
|---:|---:|---:|---:|---:|
| 6.4 | 1.032743 | 0.882804 | 0.001903 | 0.937893 |
| 7.0 | 1.146872 | 0.975093 | 0.001759 | 1.051626 |

This illustrates the different observables:

- a single strongest path can exceed one while the device-integrated
  avalanche current is still below the drain current;
- `BreakAtIonIntegral(3 1.)` does not fire because the third distinct path is
  still far below one;
- the global `Iava/Id` criterion crosses at 6.734426 V in the dense replay;
- the original path sweep continues until a third high-integral path is
  present and exceeds one, terminating at 10.448267 V.

The path field is a ranked local topology observable. The integrated source is
a device-wide amplitude observable. They should be reported as separate
validation contracts and must not share one breakdown-voltage label.

## Next root-cause target

The next current-source closure task should hold the peak model fixed and
compare the transverse/longitudinal growth of the electron-generation support:

1. area and source fraction above 10%, 30%, 50%, and 80% of each simulator's
   peak generation;
2. same-coordinate electron current-density profiles around the drain-body
   high-field junction;
3. accumulated electron source versus distance from the common peak;
4. only after those close, re-evaluate the 6.734426 V Sentaurus versus
   6.920074 V Vela current-source crossing.

## Reproducible artifacts

- analysis implementation:
  `scripts/audit_bvmethods_nmos_hotspot_slope_and_criteria.py`;
- machine summary:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/qf_vector_hotspot_slope_criteria_audit_20260805/summary.json`;
- detailed outputs in the same directory:
  `sentaurus_criteria_vs_bias.csv`, `vela_hotspot_vs_bias.csv`,
  `hotspot_same_bias_compare.csv`, and `hotspot_slope_compare.csv`;
- manual source:
  `/usr/synopsys/sentaurus/O_2018.06-SP2/tcad/O-2018.06-SP2/manuals/PDFManual/data/sdevice_ug.pdf`, approximate breakdown analysis and
  `BreakAtIonIntegral` reference entries.
